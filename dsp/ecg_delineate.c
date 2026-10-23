#include "ecg_delineate.h"
#include "ecg_config.h"
#include <math.h>

#ifndef FLAT_WIN_MS
#define FLAT_WIN_MS 20
#endif
#define FLAT_WIN    ECG_MS(FLAT_WIN_MS)  /* flat window ending the QRS; longer than a Q/S trough */
#define BASE_WIN    ECG_MS(20)
#define NOISE_MAX   96
#ifndef NOISE_K
#define NOISE_K     2.0f        /* threshold floor, x median |2nd difference| */
#endif

static int imax(int a, int b) { return a > b ? a : b; }
static int imin(int a, int b) { return a < b ? a : b; }

/* central difference, mV/sample */
static float d1(const float *x, int i)
{
    return 0.5f * (x[i + 1] - x[i - 1]);
}

/* wider, smoother difference for the low-slope P and T waves */
static float d2(const float *x, int i)
{
    return (x[i + 2] + x[i + 1] - x[i - 1] - x[i - 2]) * (1.0f / 6.0f);
}

/*
 * Median |second difference| over [a, b]. P and T waves have almost no
 * curvature at the sample scale, so this tracks broadband EMG noise only.
 */
static float noise_level(const float *x, int a, int b)
{
    float v[NOISE_MAX];
    int n = 0, i, j;

    for (i = a; i <= b && n < NOISE_MAX; ++i) {
        float s = fabsf(x[i + 1] - 2.0f * x[i] + x[i - 1]);
        for (j = n++; j > 0 && v[j - 1] > s; --j)
            v[j] = v[j - 1];
        v[j] = s;
    }
    return n ? v[n / 2] : 0.0f;
}

static float win_slope(const float *x, int a)
{
    float s = 0.0f;
    int k;

    for (k = 0; k < FLAT_WIN; ++k)
        s += fabsf(d1(x, a + k));
    return s * (1.0f / (float)FLAT_WIN);
}

static void find_qrs(const float *x, int r, int lo, int hi, int *on, int *off)
{
    float su = 0.0f, sd = 0.0f, thr, nf;
    int iu = r, id = r, i;

    /* steepest up- and down-stroke around the dominant deflection */
    for (i = imax(lo, r - ECG_MS(60)); i < r; ++i)
        if (fabsf(d1(x, i)) > su) {
            su = fabsf(d1(x, i));
            iu = i;
        }
    for (i = r + 1; i <= imin(hi, r + ECG_MS(60)); ++i)
        if (fabsf(d1(x, i)) > sd) {
            sd = fabsf(d1(x, i));
            id = i;
        }
    thr = ECG_QRS_SLOPE_FRAC * (su > sd ? su : sd);
    nf = NOISE_K * noise_level(x, imax(2, r - ECG_MS(450)), imax(2, r - ECG_MS(150)));
    if (nf > thr)
        thr = nf;

    /*
     * Walk outwards until the mean slope over a 20 ms window is flat, then
     * step to the first flat sample next to the QRS. The window length
     * carries the search across Q and S troughs; averaging tolerates EMG.
     */
    *on = lo;
    for (i = iu; i - FLAT_WIN + 1 >= lo; --i)
        if (win_slope(x, i - FLAT_WIN + 1) < thr) {
            while (i > lo && fabsf(d1(x, i)) >= thr)
                --i;
            *on = i;
            break;
        }
    *off = hi;
    for (i = id; i + FLAT_WIN - 1 <= hi; ++i)
        if (win_slope(x, i) < thr) {
            while (i < hi && fabsf(d1(x, i)) >= thr)
                ++i;
            *off = i;
            break;
        }
}

/* Mean of the flattest 20 ms window in the PR segment. */
static float find_baseline(const float *x, int on)
{
    float best = 1e30f, base = x[on];
    int s, k;

    for (s = imax(0, on - ECG_MS(60)); s + BASE_WIN <= on; ++s) {
        float mn = x[s], mx = x[s], sum = 0.0f;
        for (k = s; k < s + BASE_WIN; ++k) {
            if (x[k] < mn) mn = x[k];
            if (x[k] > mx) mx = x[k];
            sum += x[k];
        }
        if (mx - mn < best) {
            best = mx - mn;
            base = sum / (float)BASE_WIN;
        }
    }
    return base;
}

/* T end: tangent at the steepest point of the descending limb meets baseline. */
static float find_t_end(const float *x, int ts, int te, float base)
{
    float ta = 0.0f, best = 0.0f, sg, t;
    int tp = ts, is = -1, i;

    if (te - ts < ECG_MS(60))
        return -1.0f;
    for (i = ts; i <= te - ECG_MS(40); ++i)
        if (fabsf(x[i] - base) > fabsf(ta)) {
            ta = x[i] - base;
            tp = i;
        }
    if (fabsf(ta) < ECG_T_MIN_MV)
        return -1.0f;

    sg = ta > 0.0f ? 1.0f : -1.0f;
    for (i = tp + 1; i <= imin(te, tp + ECG_MS(200)); ++i)
        if (-sg * d2(x, i) > best) {
            best = -sg * d2(x, i);
            is = i;
        }
    if (is < 0)
        return -1.0f;

    t = (float)is + (base - x[is]) / d2(x, is);
    if (t < (float)is)
        t = (float)is;
    if (t > (float)te)
        t = (float)te;
    return t;
}

/* P onset: tangent at the steepest point of the ascending limb meets baseline. */
static float find_p_on(const float *x, int ps, int pe, float base)
{
    float pa = 0.0f, best = 0.0f, sg, t;
    int pp = ps, is = -1, i;

    if (pe - ps < ECG_MS(60))
        return -1.0f;
    for (i = ps; i <= pe; ++i)
        if (fabsf(x[i] - base) > fabsf(pa)) {
            pa = x[i] - base;
            pp = i;
        }
    /* must be a real wave, not the tail of T or the start of the QRS */
    if (fabsf(pa) < ECG_P_MIN_MV || pp <= ps + 1 || pp >= pe)
        return -1.0f;

    sg = pa > 0.0f ? 1.0f : -1.0f;
    for (i = imax(ps, pp - ECG_MS(100)); i < pp; ++i)
        if (sg * d2(x, i) > best) {
            best = sg * d2(x, i);
            is = i;
        }
    if (is < 0)
        return -1.0f;

    t = (float)is + (base - x[is]) / d2(x, is);
    if (t > (float)is)
        t = (float)is;
    if (t < (float)ps)
        t = (float)ps;
    return t;
}

void ecg_delineate(const float *x, int len, int r, int p_min, int t_max,
                   int st_off, ecg_delin_t *out)
{
    int lo = imax(2, r - ECG_MS(120));
    int hi = imin(imin(len - 3, t_max), r + ECG_MS(140));
    int i, j;

    find_qrs(x, r, lo, hi, &out->qrs_on, &out->qrs_off);
    out->baseline = find_baseline(x, out->qrs_on);

    out->r_amp = x[out->qrs_on];
    for (i = out->qrs_on; i <= out->qrs_off; ++i)
        if (x[i] > out->r_amp)
            out->r_amp = x[i];
    out->r_amp -= out->baseline;

    j = out->qrs_off + st_off;
    out->st_valid = j + 1 <= imin(t_max, len - 1);
    out->st = out->st_valid
        ? (x[j - 1] + x[j] + x[j + 1]) * (1.0f / 3.0f) - out->baseline
        : 0.0f;

    out->t_end = find_t_end(x, out->qrs_off + ECG_MS(60),
                            imin(t_max, len - 3), out->baseline);
    out->p_on = find_p_on(x, imax(imax(p_min, 2), out->qrs_on - ECG_MS(ECG_PR_SEARCH_MS)),
                          out->qrs_on - ECG_MS(30), out->baseline);
}
