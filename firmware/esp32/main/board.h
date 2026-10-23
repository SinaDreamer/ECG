#ifndef BOARD_H
#define BOARD_H

/*
 * AD8232 single-lead front end (e.g. SparkFun SEN-12650) on an ESP32 DevKitC.
 * The ECG input must be on ADC1; ADC2 is unavailable while Wi-Fi is active.
 * Set ECG_LO_P_GPIO / ECG_LO_N_GPIO to -1 if the lead-off pins are not wired.
 */
#define ECG_ADC_GPIO            36          /* AD8232 OUTPUT -> ADC1_CH0 (SENSOR_VP) */
#define ECG_LO_P_GPIO           34          /* AD8232 LO+, high when RA/LA lead is off */
#define ECG_LO_N_GPIO           35          /* AD8232 LO- */

#define ECG_FRONTEND_GAIN       1100.0f     /* INA 100 x output stage 11 */
#define ECG_ADC_FULL_SCALE_MV   3100.0f     /* nominal, used only without eFuse calibration */

#define ECG_UART_NUM            UART_NUM_0  /* USB bridge on DevKitC */
#define ECG_UART_BAUD           115200      /* stream uses ~5.6 kB/s, ~49% of the link */
#define ECG_UART_TX_BUF         4096

#define ECG_QUEUE_LEN           64          /* 256 ms of raw samples between tasks */

#endif
