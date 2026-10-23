#ifndef ECG_MONITOR_H
#define ECG_MONITOR_H

#include <stdint.h>
#include "ecg_config.h"
#include "ecg_filter.h"
#include "qrs_detector.h"

#define ECG_NA          INT32_MIN   /* measurement not available */
#define ECG_PEND_LEN    8
#define ECG_SEG_PRE     ECG_MS(450)
#define ECG_SEG_MAX     (ECG_SEG_PRE + ECG_MS(ECG_T_SEARCH_MAX_MS) + 8)

typedef struct {
    uint32_t r_index;   /* sample index of the R fiducial (filtered stream) */
    int32_t hr_x10;     /* heart rate, bpm x 10 */
    int32_t rr_ms;
    int32_t qrs_ms;
    int32_t pr_ms;
    int32_t qt_ms;
    int32_t r_uv;       /* R amplitude above baseline, uV */
    int32_t st_uv;      /* ST level at J + 60 ms (J + 40 ms above 120 bpm), uV */
} ecg_beat_t;

typedef struct {
    uint32_t r;
    float r_frac;       /* sub-sample R position, -0.5..0.5 */
    float rr;           /* samples, 0 if unknown */
    uint32_t t_lim;     /* next QRS bounds the T search */
    int limited;
} ecg_pending_t;

typedef struct {
    ecg_filter_t filt;
    qrs_detector_t qrs;
    float buf[ECG_BUF_LEN];
    float seg[ECG_SEG_MAX];
    uint32_t n;         /* index of the next sample */
    uint32_t first;     /* first valid index after reset */

    ecg_pending_t pend[ECG_PEND_LEN];
    int ph, pc;

    int have_prev_r;
    uint32_t prev_r;
    float prev_r_frac;
    int have_prev_tend;
    uint32_t prev_tend;
    float polarity;
} ecg_monitor_t;

/* start_index: sample index assigned to the first sample after (re)init. */
void ecg_monitor_init(ecg_monitor_t *m, uint32_t start_index);

/*
 * Process one raw sample (electrode-referred mV). Writes the filtered sample
 * to *y_mv. Returns 1 and fills *beat when a beat has been fully measured.
 * Beats are reported once the T wave has been observed, i.e. roughly
 * 0.7 * RR + 0.3 s after the R peak (bounded by ~0.9 s).
 */
int ecg_monitor_process(ecg_monitor_t *m, float x_mv, float *y_mv, ecg_beat_t *beat);

#endif
