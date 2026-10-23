#include "ecg_filter.h"

/* Section Qs of a 4th-order Butterworth: 1 / (2 cos(pi/8)), 1 / (2 cos(3pi/8)) */
#define BW4_Q1      0.54119610f
#define BW4_Q2      1.30656296f
#define BW2_Q       0.70710678f

/* 1 - sinc^2(k) = 1/sqrt(2)  ->  L = k fs / fc */
#define MA_HPF_K    0.572f

static void ma_hpf_init(ma_hpf_t *h)
{
    int len = (int)(MA_HPF_K * (float)ECG_FS_HZ / ECG_HPF_HZ + 0.5f);

    if (len > ECG_HPF_MA_MAX)
        len = ECG_HPF_MA_MAX;
    if (len < 2)
        len = 2;
    h->len = len;
    h->pos = 0;
    h->sum1 = h->sum2 = 0.0f;
}

static void ma_hpf_prime(ma_hpf_t *h, float x0)
{
    int k;

    for (k = 0; k < h->len; ++k) {
        h->b1[k] = x0;
        h->b2[k] = x0;
    }
    h->sum1 = h->sum2 = x0 * (float)h->len;
}

static float ma_hpf_run(ma_hpf_t *h, float x)
{
    float inv = 1.0f / (float)h->len;
    float m1, m2, delayed;
    int p = h->pos, k;

    h->sum1 += x - h->b1[p];
    h->b1[p] = x;
    m1 = h->sum1 * inv;

    h->sum2 += m1 - h->b2[p];
    h->b2[p] = m1;
    m2 = h->sum2 * inv;

    if (++p == h->len) {
        /* resum once per window to stop float drift */
        p = 0;
        h->sum1 = h->sum2 = 0.0f;
        for (k = 0; k < h->len; ++k) {
            h->sum1 += h->b1[k];
            h->sum2 += h->b2[k];
        }
    }
    h->pos = p;

    /* oldest sample in b1 is x[n - (L-1)], matching the MA^2 group delay */
    delayed = h->b1[p];
    return delayed - m2;
}

void ecg_filter_init(ecg_filter_t *f)
{
#if ECG_HPF_LINEAR_PHASE
    ma_hpf_init(&f->hp);
#else
    biquad_highpass(&f->hp, (float)ECG_FS_HZ, ECG_HPF_HZ, BW2_Q);
#endif
    biquad_notch(&f->notch, (float)ECG_FS_HZ, ECG_NOTCH_HZ, ECG_NOTCH_Q);
    biquad_lowpass(&f->lp[0], (float)ECG_FS_HZ, ECG_LPF_HZ, BW4_Q1);
    biquad_lowpass(&f->lp[1], (float)ECG_FS_HZ, ECG_LPF_HZ, BW4_Q2);
    f->primed = 0;
}

float ecg_filter_run(ecg_filter_t *f, float x)
{
    float y;

    if (!f->primed) {
        /* front-end DC offset would otherwise ring the HPF for seconds */
#if ECG_HPF_LINEAR_PHASE
        ma_hpf_prime(&f->hp, x);
#else
        biquad_prime(&f->hp, x);
#endif
        f->primed = 1;
    }
#if ECG_HPF_LINEAR_PHASE
    y = ma_hpf_run(&f->hp, x);
#else
    y = biquad_run(&f->hp, x);
#endif
    y = biquad_run(&f->notch, y);
    y = biquad_run(&f->lp[0], y);
    y = biquad_run(&f->lp[1], y);
    return y;
}

int ecg_filter_delay(const ecg_filter_t *f)
{
#if ECG_HPF_LINEAR_PHASE
    return f->hp.len - 1;
#else
    (void)f;
    return 0;
#endif
}
