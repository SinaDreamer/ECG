#include "biquad.h"
#include <math.h>

#define BQ_PI 3.14159265358979f

/* Bilinear-transform designs with frequency pre-warping (RBJ cookbook). */

static void bq_set(biquad_t *f, float b0, float b1, float b2,
                   float a0, float a1, float a2)
{
    f->b0 = b0 / a0;
    f->b1 = b1 / a0;
    f->b2 = b2 / a0;
    f->a1 = a1 / a0;
    f->a2 = a2 / a0;
    f->s1 = 0.0f;
    f->s2 = 0.0f;
}

void biquad_lowpass(biquad_t *f, float fs, float fc, float q)
{
    float w0 = 2.0f * BQ_PI * fc / fs;
    float c = cosf(w0);
    float alpha = sinf(w0) / (2.0f * q);

    bq_set(f, 0.5f * (1.0f - c), 1.0f - c, 0.5f * (1.0f - c),
           1.0f + alpha, -2.0f * c, 1.0f - alpha);
}

void biquad_highpass(biquad_t *f, float fs, float fc, float q)
{
    float w0 = 2.0f * BQ_PI * fc / fs;
    float c = cosf(w0);
    float alpha = sinf(w0) / (2.0f * q);

    bq_set(f, 0.5f * (1.0f + c), -(1.0f + c), 0.5f * (1.0f + c),
           1.0f + alpha, -2.0f * c, 1.0f - alpha);
}

void biquad_notch(biquad_t *f, float fs, float f0, float q)
{
    float w0 = 2.0f * BQ_PI * f0 / fs;
    float c = cosf(w0);
    float alpha = sinf(w0) / (2.0f * q);

    bq_set(f, 1.0f, -2.0f * c, 1.0f,
           1.0f + alpha, -2.0f * c, 1.0f - alpha);
}

void biquad_prime(biquad_t *f, float x0)
{
    float dc = (f->b0 + f->b1 + f->b2) / (1.0f + f->a1 + f->a2);
    float y = dc * x0;

    f->s2 = f->b2 * x0 - f->a2 * y;
    f->s1 = f->b1 * x0 - f->a1 * y + f->s2;
}
