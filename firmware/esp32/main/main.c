/*
 * ECG monitor: AD8232 -> ADC1 @ 250 Hz -> filter / QRS / delineation -> UART.
 *
 * Sampling path: GPTimer alarm ISR (250 Hz) notifies a high-priority task
 * that reads the ADC and lead-off pins and queues the sample. Nothing in this
 * path waits on anything except the timer. DSP and UART output run in a
 * separate task fed by the queue, so a slow host can never stall sampling.
 */
#include <stdbool.h>
#include <stdint.h>

#include "freertos/FreeRTOS.h"
#include "freertos/queue.h"
#include "freertos/task.h"
#include "driver/gpio.h"
#include "driver/gptimer.h"
#include "driver/uart.h"
#include "esp_adc/adc_cali.h"
#include "esp_adc/adc_cali_scheme.h"
#include "esp_adc/adc_oneshot.h"
#include "esp_attr.h"
#include "esp_err.h"
#include "esp_idf_version.h"
#include "esp_log.h"
#include "sdkconfig.h"

#include "board.h"
#include "ecg_config.h"
#include "ecg_monitor.h"
#include "ecg_report.h"

#if ESP_IDF_VERSION >= ESP_IDF_VERSION_VAL(5, 2, 0)
#define ECG_ADC_ATTEN       ADC_ATTEN_DB_12
#else
#define ECG_ADC_ATTEN       ADC_ATTEN_DB_11
#endif

#if CONFIG_FREERTOS_UNICORE
#define SAMPLER_CORE        0
#else
#define SAMPLER_CORE        1
#endif
#define PROCESS_CORE        0

#define LEAD_DEBOUNCE       ECG_MS(100)
#define GAP_FILL_MAX        ECG_MS(100)     /* longer gaps restart the analysis */
#define STATS_PERIOD        (10 * ECG_FS_HZ)
#define HAVE_LEAD_OFF       (ECG_LO_P_GPIO >= 0 && ECG_LO_N_GPIO >= 0)

typedef struct {
    uint32_t seq;
    int16_t raw;            /* ADC counts, -1 on read error */
    uint8_t lead_off;
} raw_sample_t;

static const char *TAG = "ecg";

static TaskHandle_t s_sampler;
static QueueHandle_t s_queue;
static adc_oneshot_unit_handle_t s_adc;
static adc_channel_t s_adc_chan;
static adc_cali_handle_t s_cali;
static bool s_cali_ok;

static volatile uint32_t s_missed_ticks;    /* sampler woke late */
static volatile uint32_t s_queue_full;      /* processing fell behind */
static uint32_t s_tx_dropped;               /* host link saturated */

static ecg_monitor_t s_mon;

/* ------------------------------------------------------------ sampling */

static bool IRAM_ATTR on_sample_tick(gptimer_handle_t timer,
                                     const gptimer_alarm_event_data_t *edata,
                                     void *ctx)
{
    BaseType_t woken = pdFALSE;

    (void)timer;
    (void)edata;
    (void)ctx;
    vTaskNotifyGiveFromISR(s_sampler, &woken);
    return woken == pdTRUE;
}

static void sampler_task(void *arg)
{
    uint32_t seq = 0;

    (void)arg;
    for (;;) {
        uint32_t ticks = ulTaskNotifyTake(pdTRUE, portMAX_DELAY);
        raw_sample_t s;
        int raw = -1;

        if (adc_oneshot_read(s_adc, s_adc_chan, &raw) != ESP_OK)
            raw = -1;
        s.raw = (int16_t)raw;
#if HAVE_LEAD_OFF
        s.lead_off = (uint8_t)(gpio_get_level(ECG_LO_P_GPIO) | gpio_get_level(ECG_LO_N_GPIO));
#else
        s.lead_off = 0;
#endif
        if (ticks > 1)
            s_missed_ticks += ticks - 1;

        /* replicate on missed ticks so the sample clock stays exactly 250 Hz */
        while (ticks--) {
            s.seq = seq++;
            if (xQueueSend(s_queue, &s, 0) != pdTRUE)
                s_queue_full++;
        }
    }
}

/* ---------------------------------------------------------- processing */

static void uart_send(const char *line, size_t len, bool droppable)
{
    size_t room = 0;

    if (len == 0)
        return;
    /* sample lines are dropped rather than queued when the link is saturated */
    if (droppable &&
        uart_get_tx_buffer_free_size(ECG_UART_NUM, &room) == ESP_OK && room < len) {
        s_tx_dropped++;
        return;
    }
    uart_write_bytes(ECG_UART_NUM, line, len);
}

static void send_event(uint32_t n, const char *code, int32_t value)
{
    char line[ECG_LINE_MAX];

    uart_send(line, ecg_fmt_event(line, sizeof(line), n, code, value), false);
}

static float raw_to_mv(int raw)
{
    int pin_mv;

    if (s_cali_ok && adc_cali_raw_to_voltage(s_cali, raw, &pin_mv) == ESP_OK)
        return (float)pin_mv / ECG_FRONTEND_GAIN;
    return (float)raw * (ECG_ADC_FULL_SCALE_MV / 4095.0f) / ECG_FRONTEND_GAIN;
}

static int32_t mv_to_uv(float mv)
{
    float uv = mv * 1000.0f;
    return (int32_t)(uv < 0.0f ? uv - 0.5f : uv + 0.5f);
}

static void feed(float mv)
{
    char line[ECG_LINE_MAX];
    ecg_beat_t beat;
    uint32_t n = s_mon.n;
    float y;

    if (ecg_monitor_process(&s_mon, mv, &y, &beat))
        uart_send(line, ecg_fmt_beat(line, sizeof(line), &beat), false);
    uart_send(line, ecg_fmt_sample(line, sizeof(line), n, mv_to_uv(y)), true);
}

static void report_counters(uint32_t n)
{
    static uint32_t missed, full, tx;
    uint32_t m = s_missed_ticks, f = s_queue_full, t = s_tx_dropped;

    if (m != missed)
        send_event(n, "OVR", (int32_t)(m - missed));
    if (f != full)
        send_event(n, "DROP", (int32_t)(f - full));
    if (t != tx)
        send_event(n, "TXD", (int32_t)(t - tx));
    missed = m;
    full = f;
    tx = t;
}

static void process_task(void *arg)
{
    char line[ECG_LINE_MAX];
    raw_sample_t s;
    uint32_t expect = 0, since_stats = 0;
    bool leads_off = false;
    int lo_count = 0;
    float last_mv = 0.0f;

    (void)arg;
    ecg_monitor_init(&s_mon, 0);
    uart_send(line, ecg_fmt_info(line, sizeof(line), ecg_filter_delay(&s_mon.filt)), false);

    for (;;) {
        int32_t gap;

        xQueueReceive(s_queue, &s, portMAX_DELAY);

        if (++since_stats >= STATS_PERIOD) {
            report_counters(s.seq);
            since_stats = 0;
        }

        /* debounced lead-off; reconnect restarts filters and thresholds */
        if ((bool)s.lead_off != leads_off) {
            if (++lo_count >= LEAD_DEBOUNCE) {
                leads_off = s.lead_off;
                lo_count = 0;
                send_event(s.seq, leads_off ? "LOFF" : "LON", ECG_NA);
                if (!leads_off) {
                    ecg_monitor_init(&s_mon, s.seq);
                    expect = s.seq;
                }
            }
        } else {
            lo_count = 0;
        }
        if (leads_off || s.raw < 0) {
            expect = s.seq + 1;
            continue;
        }

        /* bridge short queue overflows by holding the last value */
        gap = (int32_t)(s.seq - expect);
        if (gap > GAP_FILL_MAX) {
            ecg_monitor_init(&s_mon, s.seq);
        } else {
            while (gap-- > 0)
                feed(last_mv);
        }
        expect = s.seq + 1;

        last_mv = raw_to_mv(s.raw);
        feed(last_mv);
    }
}

/* --------------------------------------------------------------- setup */

static bool cali_setup(adc_unit_t unit, adc_channel_t chan)
{
    esp_err_t err = ESP_ERR_NOT_SUPPORTED;

#if ADC_CALI_SCHEME_CURVE_FITTING_SUPPORTED
    adc_cali_curve_fitting_config_t cfg = {
        .unit_id = unit,
        .chan = chan,
        .atten = ECG_ADC_ATTEN,
        .bitwidth = ADC_BITWIDTH_DEFAULT,
    };
    err = adc_cali_create_scheme_curve_fitting(&cfg, &s_cali);
#elif ADC_CALI_SCHEME_LINE_FITTING_SUPPORTED
    adc_cali_line_fitting_config_t cfg = {
        .unit_id = unit,
        .atten = ECG_ADC_ATTEN,
        .bitwidth = ADC_BITWIDTH_DEFAULT,
#if CONFIG_IDF_TARGET_ESP32
        .default_vref = 1100,
#endif
    };
    (void)chan;
    err = adc_cali_create_scheme_line_fitting(&cfg, &s_cali);
#else
    (void)unit;
    (void)chan;
#endif
    if (err != ESP_OK)
        ESP_LOGW(TAG, "ADC calibration unavailable (%s), using nominal scale",
                 esp_err_to_name(err));
    return err == ESP_OK;
}

static void adc_setup(void)
{
    adc_unit_t unit;
    adc_oneshot_unit_init_cfg_t ucfg = { 0 };
    adc_oneshot_chan_cfg_t ccfg = {
        .atten = ECG_ADC_ATTEN,
        .bitwidth = ADC_BITWIDTH_DEFAULT,
    };

    ESP_ERROR_CHECK(adc_oneshot_io_to_channel(ECG_ADC_GPIO, &unit, &s_adc_chan));
    if (unit != ADC_UNIT_1)
        ESP_LOGW(TAG, "GPIO%d is on ADC2; readings fail while Wi-Fi is active", ECG_ADC_GPIO);
    ucfg.unit_id = unit;
    ucfg.ulp_mode = ADC_ULP_MODE_DISABLE;
    ESP_ERROR_CHECK(adc_oneshot_new_unit(&ucfg, &s_adc));
    ESP_ERROR_CHECK(adc_oneshot_config_channel(s_adc, s_adc_chan, &ccfg));
    s_cali_ok = cali_setup(unit, s_adc_chan);
}

static void lead_off_setup(void)
{
#if HAVE_LEAD_OFF
    gpio_config_t io = {
        .pin_bit_mask = (1ULL << ECG_LO_P_GPIO) | (1ULL << ECG_LO_N_GPIO),
        .mode = GPIO_MODE_INPUT,
        .pull_up_en = GPIO_PULLUP_DISABLE,      /* AD8232 drives LO+/LO- push-pull */
        .pull_down_en = GPIO_PULLDOWN_DISABLE,
        .intr_type = GPIO_INTR_DISABLE,
    };
    ESP_ERROR_CHECK(gpio_config(&io));
#endif
}

static void uart_setup(void)
{
    uart_config_t cfg = {
        .baud_rate = ECG_UART_BAUD,
        .data_bits = UART_DATA_8_BITS,
        .parity = UART_PARITY_DISABLE,
        .stop_bits = UART_STOP_BITS_1,
        .flow_ctrl = UART_HW_FLOWCTRL_DISABLE,
        .source_clk = UART_SCLK_DEFAULT,
    };

    ESP_ERROR_CHECK(uart_driver_install(ECG_UART_NUM, 256, ECG_UART_TX_BUF, 0, NULL, 0));
    ESP_ERROR_CHECK(uart_param_config(ECG_UART_NUM, &cfg));
}

static void timer_setup(void)
{
    gptimer_handle_t timer;
    gptimer_config_t tcfg = {
        .clk_src = GPTIMER_CLK_SRC_DEFAULT,
        .direction = GPTIMER_COUNT_UP,
        .resolution_hz = 1000000,
    };
    gptimer_alarm_config_t acfg = {
        .alarm_count = 1000000 / ECG_FS_HZ,
        .reload_count = 0,
        .flags.auto_reload_on_alarm = true,
    };
    gptimer_event_callbacks_t cbs = {
        .on_alarm = on_sample_tick,
    };

    ESP_ERROR_CHECK(gptimer_new_timer(&tcfg, &timer));
    ESP_ERROR_CHECK(gptimer_register_event_callbacks(timer, &cbs, NULL));
    ESP_ERROR_CHECK(gptimer_set_alarm_action(timer, &acfg));
    ESP_ERROR_CHECK(gptimer_enable(timer));
    ESP_ERROR_CHECK(gptimer_start(timer));
}

void app_main(void)
{
    uart_setup();
    adc_setup();
    lead_off_setup();

    /* console shares the data UART; keep it quiet once setup has reported */
    esp_log_level_set("*", ESP_LOG_ERROR);

    s_queue = xQueueCreate(ECG_QUEUE_LEN, sizeof(raw_sample_t));
    configASSERT(s_queue);
    xTaskCreatePinnedToCore(process_task, "ecg_dsp", 4096, NULL, 5, NULL, PROCESS_CORE);
    xTaskCreatePinnedToCore(sampler_task, "ecg_adc", 3072, NULL,
                            configMAX_PRIORITIES - 2, &s_sampler, SAMPLER_CORE);
    configASSERT(s_sampler);

    timer_setup();
}
