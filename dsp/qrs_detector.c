#include "qrs_detector.h"
#include <math.h>
#include <string.h>

#define REFRACTORY  ECG_MS(QRS_REFRACTORY_MS)
#define TWAVE_WIN   ECG_MS(QRS_TWAVE_MS)
#define LEARN_LEN   ECG_MS(QRS_LEARN_MS)
#define CONFIRM     ECG_MS(QRS_PEAK_CONFIRM_MS)
#define SPKI_DECAY  ECG_MS(QRS_SPKI_DECAY_MS)

void qrs_init(qrs_detector_t *d, uint32_t start_index)
{
    memset(d, 0, sizeof(*d));
    biquad_highpass(&d->hp, (float)ECG_FS_HZ, QRS_BP_LO_HZ, 0.70710678f);
    biquad_lowpass(&d->lp, (float)ECG_FS_HZ, QRS_BP_HI_HZ, 0.70710678f);
    d->n = start_index;
    d->learning = 1;
}

static void update_thresholds(qrs_detector_t *d)
{
    d->thr1 = d->npki + 0.25f * (d->spki - d->npki);
    d->thr2 = 0.5f * d->thr1;
}

/* Largest |derivative| over the integration window ending at idx. */
static float window_slope(const qrs_detector_t *d, uint32_t idx)
{
    float m = 0.0f;
    int k;

    for (k = 0; k < QRS_MWI_LEN; ++k) {
        float s = d->slope[(idx - (uint32_t)k) & (QRS_SLOPE_RING - 1)];
        if (s > m)
            m = s;
    }
    return m;
}

static void accept_qrs(qrs_detector_t *d, uint32_t idx, float slope)
{
    if (d->have_qrs) {
        float rr = (float)(int32_t)(idx - d->last_qrs);
        float sum = 0.0f;
        int k;

        d->rr[d->rr_pos] = rr;
        d->rr_pos = (d->rr_pos + 1) % QRS_RR_AVG_N;
        if (d->rr_count < QRS_RR_AVG_N)
            d->rr_count++;
        for (k = 0; k < d->rr_count; ++k)
            sum += d->rr[k];
        d->rr_avg = sum / (float)d->rr_count;
    }
    d->have_qrs = 1;
    d->last_qrs = idx;
    d->last_slope = slope;
    d->sb_valid = 0;
    d->quiet = 0;
}

static int is_twave(const qrs_detector_t *d, uint32_t idx, float slope)
{
    return d->have_qrs &&
           (int32_t)(idx - d->last_qrs) < TWAVE_WIN &&
           slope < 0.5f * d->last_slope;
}

/* Returns 1 if the MWI peak is a QRS complex. */
static int classify(qrs_detector_t *d, float pk, uint32_t idx)
{
    float slope;

    if (d->have_qrs && (int32_t)(idx - d->last_qrs) < REFRACTORY)
        return 0;

    slope = window_slope(d, idx);

    if (pk > d->thr1 && !is_twave(d, idx, slope)) {
        d->spki = 0.125f * pk + 0.875f * d->spki;
        update_thresholds(d);
        accept_qrs(d, idx, slope);
        return 1;
    }

    if (pk > d->thr2 && !is_twave(d, idx, slope) &&
        (!d->sb_valid || pk > d->sb_pk)) {
        d->sb_valid = 1;
        d->sb_pk = pk;
        d->sb_idx = idx;
        d->sb_slope = slope;
    }
    d->npki = 0.125f * pk + 0.875f * d->npki;
    update_thresholds(d);
    return 0;
}

int qrs_process(qrs_detector_t *d, float x, uint32_t out[2])
{
    uint32_t n = d->n++;
    float *h = d->hist;
    float y, dv, sq, mwi;
    int nout = 0;
    int k;

    /* 5-15 Hz bandpass */
    y = biquad_run(&d->lp, biquad_run(&d->hp, x));

    /* five-point derivative: (2x[n] + x[n-1] - x[n-3] - 2x[n-4]) / 8T */
    h[4] = h[3];
    h[3] = h[2];
    h[2] = h[1];
    h[1] = h[0];
    h[0] = y;
    dv = (2.0f * h[0] + h[1] - h[3] - 2.0f * h[4]) * ((float)ECG_FS_HZ / 8.0f);
    d->slope[n & (QRS_SLOPE_RING - 1)] = fabsf(dv);

    /* squaring and moving-window integration */
    sq = dv * dv;
    d->mwi_sum += sq - d->mwi_buf[d->mwi_pos];
    d->mwi_buf[d->mwi_pos] = sq;
    if (++d->mwi_pos == QRS_MWI_LEN) {
        /* resum once per window to stop float drift */
        d->mwi_pos = 0;
        d->mwi_sum = 0.0f;
        for (k = 0; k < QRS_MWI_LEN; ++k)
            d->mwi_sum += d->mwi_buf[k];
    }
    mwi = d->mwi_sum > 0.0f ? d->mwi_sum / (float)QRS_MWI_LEN : 0.0f;

    if (d->learning) {
        if (mwi > d->learn_max)
            d->learn_max = mwi;
        d->learn_sum += mwi;
        if (++d->learn_count >= LEARN_LEN) {
            d->spki = d->learn_max / 3.0f;
            d->npki = 0.5f * d->learn_sum / (float)d->learn_count;
            update_thresholds(d);
            d->learning = 0;
        }
        d->prev_mwi = mwi;
        return 0;
    }

    /* searchback: no beat for 166% of the mean RR -> take best candidate */
    if (d->have_qrs && d->rr_avg > 0.0f && d->sb_valid &&
        !(d->armed && d->pk_max > d->thr1)) {
        int32_t since = (int32_t)(n - (uint32_t)CONFIRM - d->last_qrs);
        if ((float)since > QRS_SEARCHBACK_RR * d->rr_avg) {
            d->spki = 0.25f * d->sb_pk + 0.75f * d->spki;
            update_thresholds(d);
            accept_qrs(d, d->sb_idx, d->sb_slope);
            out[nout++] = d->last_qrs;
        }
    }

    /* signal amplitude dropped (electrode moved, gain change): relax SPKI */
    if (++d->quiet >= SPKI_DECAY) {
        d->spki *= 0.5f;
        update_thresholds(d);
        d->quiet = 0;
    }

    /* MWI peak: local max, confirmed after it halves or CONFIRM elapses */
    if (!d->armed) {
        if (mwi > d->prev_mwi) {
            d->armed = 1;
            d->pk_max = mwi;
            d->pk_idx = n;
        }
    } else if (mwi > d->pk_max) {
        d->pk_max = mwi;
        d->pk_idx = n;
    } else if (mwi < 0.5f * d->pk_max || (int32_t)(n - d->pk_idx) >= CONFIRM) {
        d->armed = 0;
        if (classify(d, d->pk_max, d->pk_idx))
            out[nout++] = d->pk_idx;
    }
    d->prev_mwi = mwi;

    return nout;
}
