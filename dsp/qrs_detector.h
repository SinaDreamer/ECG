#ifndef QRS_DETECTOR_H
#define QRS_DETECTOR_H

#include <stdint.h>
#include "biquad.h"
#include "ecg_config.h"

#define QRS_MWI_LEN     ECG_MS(QRS_MWI_MS)
#define QRS_SLOPE_RING  128     /* power of two, > MWI + peak confirm delay */
#define QRS_RR_AVG_N    8

/*
 * Pan-Tompkins QRS detector: 5-15 Hz bandpass, five-point derivative,
 * squaring, 150 ms moving-window integration, dual adaptive thresholds,
 * 200 ms refractory period, T-wave slope discrimination and searchback.
 */
typedef struct {
    biquad_t hp, lp;
    float hist[5];
    float mwi_buf[QRS_MWI_LEN];
    float mwi_sum;
    int mwi_pos;
    float slope[QRS_SLOPE_RING];
    uint32_t n;

    /* MWI peak finder */
    int armed;
    float pk_max;
    uint32_t pk_idx;
    float prev_mwi;

    /* adaptive thresholds */
    float spki, npki, thr1, thr2;
    int learning;
    float learn_max, learn_sum;
    uint32_t learn_count;
    uint32_t quiet;

    /* beat history */
    int have_qrs;
    uint32_t last_qrs;
    float last_slope;
    float rr[QRS_RR_AVG_N];
    int rr_count, rr_pos;
    float rr_avg;

    /* searchback candidate */
    int sb_valid;
    float sb_pk, sb_slope;
    uint32_t sb_idx;
} qrs_detector_t;

void qrs_init(qrs_detector_t *d, uint32_t start_index);

/*
 * Feed one filtered ECG sample. Returns the number of QRS complexes detected
 * (0..2) and writes their MWI peak sample indices to out[].
 */
int qrs_process(qrs_detector_t *d, float x, uint32_t out[2]);

/* Worst-case delay from an MWI peak to its detection being reported. */
#define QRS_DETECT_LATENCY  (QRS_MWI_LEN + ECG_MS(QRS_PEAK_CONFIRM_MS))

#endif
