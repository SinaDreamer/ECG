#ifndef BIQUAD_H
#define BIQUAD_H

/* Second-order IIR section, transposed direct form II, a0 normalised to 1. */
typedef struct {
    float b0, b1, b2, a1, a2;
    float s1, s2;
} biquad_t;

void biquad_lowpass(biquad_t *f, float fs, float fc, float q);
void biquad_highpass(biquad_t *f, float fs, float fc, float q);
void biquad_notch(biquad_t *f, float fs, float f0, float q);

/* Load the state as if x0 had been applied forever (no start-up transient). */
void biquad_prime(biquad_t *f, float x0);

static inline float biquad_run(biquad_t *f, float x)
{
    float y = f->b0 * x + f->s1;
    f->s1 = f->b1 * x - f->a1 * y + f->s2;
    f->s2 = f->b2 * x - f->a2 * y;
    return y;
}

#endif
