#ifndef ECG_FILTER_H
#define ECG_FILTER_H

#include "biquad.h"
#include "ecg_config.h"

/*
 * Linear-phase high-pass: y[n] = x[n - (L-1)] - MA_L(MA_L(x))[n].
 * Two cascaded length-L boxcars give a non-negative sinc^2 low-pass, so the
 * high-pass passband never exceeds unity (ripple < 0.42 dB above 2 fc).
 */
typedef struct {
    float b1[ECG_HPF_MA_MAX];
    float b2[ECG_HPF_MA_MAX];
    float sum1, sum2;
    int pos, len;
} ma_hpf_t;

/* HPF 0.5 Hz -> notch 60 Hz -> LPF 40 Hz (4th-order Butterworth). */
typedef struct {
#if ECG_HPF_LINEAR_PHASE
    ma_hpf_t hp;
#else
    biquad_t hp;
#endif
    biquad_t notch;
    biquad_t lp[2];
    int primed;
} ecg_filter_t;

void ecg_filter_init(ecg_filter_t *f);
float ecg_filter_run(ecg_filter_t *f, float x);

/* Constant group delay of the high-pass stage, in samples. */
int ecg_filter_delay(const ecg_filter_t *f);

#endif
