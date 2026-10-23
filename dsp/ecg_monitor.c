#include "ecg_monitor.h"
#include "ecg_delineate.h"
#include <string.h>

#define BUF_MASK        (ECG_BUF_LEN - 1)
#define R_SEARCH        (QRS_MWI_LEN + ECG_MS(40))
#define R_DEDUPE        ECG_MS(160)
#define FINAL_LATENCY   (QRS_DETECT_LATENCY + ECG_MS(40))
#define MS_PER_SAMPLE   (1000.0f / (float)ECG_FS_HZ)

static int32_t rnd(float v)
{
    return (int32_t)(v < 0.0f ? v - 0.5f : v + 0.5f);
}

void ecg_monitor_init(ecg_monitor_t *m, uint32_t start_index)
{
    memset(m, 0, sizeof(*m));
    ecg_filter_init(&m->filt);
    qrs_init(&m->qrs, start_index);
    m->n = start_index;
    m->first = start_index;
}

static float at(const ecg_monitor_t *m, uint32_t i)
{
    return m->buf[i & BUF_MASK];
}

/* Oldest index still held in the history buffer (n = newest). */
static uint32_t oldest(const ecg_monitor_t *m, uint32_t n)
{
    uint32_t avail = n - m->first + 1;
    if (avail > ECG_BUF_LEN)
        avail = ECG_BUF_LEN;
    return n - avail + 1;
}

static void on_detection(ecg_monitor_t *m, uint32_t p, uint32_t n)
{
    uint32_t lo = p - R_SEARCH, i, ip = p, in = p, r;
    uint32_t old = oldest(m, n);
    float vp, vn, frac = 0.0f;
    ecg_pending_t *e;

    if ((int32_t)(lo - old) < 1)
        lo = old + 1;

    /* dominant deflection; polarity is tracked so RS complexes don't flip */
    vp = vn = at(m, p);
    for (i = lo; (int32_t)(i - p) <= 0; ++i) {
        float v = at(m, i);
        if (v > vp) { vp = v; ip = i; }
        if (v < vn) { vn = v; in = i; }
    }
    m->polarity = 0.8f * m->polarity + 0.2f * (vp + vn);
    r = m->polarity >= 0.0f ? ip : in;

    if (m->have_prev_r && (int32_t)(r - m->prev_r) < R_DEDUPE)
        return;

    /* parabolic interpolation of the peak for sub-sample RR */
    if ((int32_t)(n - r) >= 1 && (int32_t)(r - old) >= 1) {
        float y0 = at(m, r - 1), y1 = at(m, r), y2 = at(m, r + 1);
        float den = y0 - 2.0f * y1 + y2;
        if (den != 0.0f)
            frac = 0.5f * (y0 - y2) / den;
        if (frac > 0.5f) frac = 0.5f;
        if (frac < -0.5f) frac = -0.5f;
    }

    if (m->pc > 0) {
        e = &m->pend[(m->ph + m->pc - 1) % ECG_PEND_LEN];
        e->limited = 1;
        e->t_lim = r - ECG_MS(80);
    }
    if (m->pc == ECG_PEND_LEN) {
        m->ph = (m->ph + 1) % ECG_PEND_LEN;
        m->pc--;
    }
    e = &m->pend[(m->ph + m->pc) % ECG_PEND_LEN];
    m->pc++;
    e->r = r;
    e->r_frac = frac;
    e->rr = m->have_prev_r
        ? (float)(int32_t)(r - m->prev_r) + frac - m->prev_r_frac
        : 0.0f;
    e->limited = 0;

    m->have_prev_r = 1;
    m->prev_r = r;
    m->prev_r_frac = frac;
}

static int32_t ms_or_na(float samples, int32_t lo_ms, int32_t hi_ms)
{
    int32_t ms = rnd(samples * MS_PER_SAMPLE);
    return (ms >= lo_ms && ms <= hi_ms) ? ms : ECG_NA;
}

static int finalize(ecg_monitor_t *m, uint32_t n, ecg_beat_t *beat)
{
    ecg_pending_t *e;
    ecg_delin_t dl;
    int32_t span;
    uint32_t t_end, s0, s1, old, i;
    int len, p_min, st_off;

    if (m->pc == 0)
        return 0;
    e = &m->pend[m->ph];

    span = ECG_MS(ECG_T_SEARCH_MAX_MS);
    if (e->rr > 0.0f && ECG_T_SEARCH_RR_FRAC * e->rr < (float)span)
        span = (int32_t)(ECG_T_SEARCH_RR_FRAC * e->rr);
    t_end = e->r + (uint32_t)span;
    if (e->limited && (int32_t)(e->t_lim - t_end) < 0)
        t_end = e->t_lim;

    if (!(e->limited && (int32_t)(n - t_end) >= 3) &&
        (int32_t)(n - (e->r + (uint32_t)span)) < FINAL_LATENCY)
        return 0;

    m->ph = (m->ph + 1) % ECG_PEND_LEN;
    m->pc--;

    /* copy the beat's segment out of the ring buffer */
    old = oldest(m, n);
    s0 = e->r - ECG_SEG_PRE;
    if ((int32_t)(s0 - old) < 0)
        s0 = old;
    s1 = t_end + 3;
    if ((int32_t)(s1 - n) > 0)
        s1 = n;
    len = (int)(int32_t)(s1 - s0) + 1;
    if (len > ECG_SEG_MAX)
        len = ECG_SEG_MAX;
    if ((int32_t)(e->r - s0) < 2 || (int32_t)(t_end - s0) + 3 > len)
        return 0;
    for (i = 0; i < (uint32_t)len; ++i)
        m->seg[i] = at(m, s0 + i);

    p_min = 2;
    if (m->have_prev_tend) {
        int32_t d = (int32_t)(m->prev_tend - s0) + ECG_MS(20);
        if (d > p_min)
            p_min = (int)d;
    }

    st_off = (e->rr > 0.0f && e->rr < (float)ECG_MS(500))
        ? ECG_MS(ECG_ST_OFFSET_FAST_MS) : ECG_MS(ECG_ST_OFFSET_MS);
    ecg_delineate(m->seg, len, (int)(e->r - s0), p_min, (int)(t_end - s0),
                  st_off, &dl);

    beat->r_index = e->r;
    if (e->rr > 0.0f) {
        float rr_ms = e->rr * MS_PER_SAMPLE;
        beat->rr_ms = rnd(rr_ms);
        beat->hr_x10 = rnd(600000.0f / rr_ms);
    } else {
        beat->rr_ms = ECG_NA;
        beat->hr_x10 = ECG_NA;
    }
    beat->qrs_ms = ms_or_na((float)(dl.qrs_off - dl.qrs_on), 30, 220);
    beat->pr_ms = dl.p_on >= 0.0f
        ? ms_or_na((float)dl.qrs_on - dl.p_on, 60, 400) : ECG_NA;
    beat->qt_ms = dl.t_end >= 0.0f
        ? ms_or_na(dl.t_end - (float)dl.qrs_on, 150, 700) : ECG_NA;
    beat->r_uv = rnd(dl.r_amp * 1000.0f);
    beat->st_uv = dl.st_valid ? rnd(dl.st * 1000.0f) : ECG_NA;

    m->have_prev_tend = 1;
    m->prev_tend = dl.t_end >= 0.0f
        ? s0 + (uint32_t)dl.t_end
        : s0 + (uint32_t)dl.qrs_off + ECG_MS(100);
    return 1;
}

int ecg_monitor_process(ecg_monitor_t *m, float x_mv, float *y_mv, ecg_beat_t *beat)
{
    uint32_t n = m->n, det[2];
    float y = ecg_filter_run(&m->filt, x_mv);
    int nd, k, got;

    m->buf[n & BUF_MASK] = y;
    if (y_mv)
        *y_mv = y;

    nd = qrs_process(&m->qrs, y, det);
    for (k = 0; k < nd; ++k)
        on_detection(m, det[k], n);

    got = finalize(m, n, beat);
    m->n = n + 1;
    return got;
}
