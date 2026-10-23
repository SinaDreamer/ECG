#include "ecg_report.h"
#include <stdio.h>

typedef struct {
    char *buf;
    size_t cap, len;
    int ok;
} line_t;

static void ln_put(line_t *l, const char *s)
{
    int w;

    if (!l->ok)
        return;
    w = snprintf(l->buf + l->len, l->cap - l->len, "%s", s);
    if (w < 0 || (size_t)w >= l->cap - l->len)
        l->ok = 0;
    else
        l->len += (size_t)w;
}

/* ",<v>" with a fixed number of decimals (v scaled by 10^dec), or "," if NA */
static void ln_field(line_t *l, int32_t v, int dec)
{
    char tmp[16];
    int64_t a;
    int32_t scale = 1;
    int k;

    if (v == ECG_NA) {
        ln_put(l, ",");
        return;
    }
    for (k = 0; k < dec; ++k)
        scale *= 10;
    a = v < 0 ? -(int64_t)v : (int64_t)v;
    if (dec == 0)
        snprintf(tmp, sizeof(tmp), ",%s%ld", v < 0 ? "-" : "", (long)a);
    else
        snprintf(tmp, sizeof(tmp), ",%s%ld.%0*ld", v < 0 ? "-" : "",
                 (long)(a / scale), dec, (long)(a % scale));
    ln_put(l, tmp);
}

static size_t ln_end(line_t *l)
{
    char tail[8];
    unsigned cs = 0;
    size_t i;

    for (i = 1; l->ok && i < l->len; ++i)
        cs ^= (unsigned char)l->buf[i];
    snprintf(tail, sizeof(tail), "*%02X\r\n", cs & 0xFFu);
    ln_put(l, tail);
    if (!l->ok) {
        if (l->cap)
            l->buf[0] = '\0';
        return 0;
    }
    return l->len;
}

static line_t ln_begin(char *buf, size_t cap, const char *tag)
{
    line_t l;

    l.buf = buf;
    l.cap = cap;
    l.len = 0;
    l.ok = cap > 0;
    ln_put(&l, tag);
    return l;
}

size_t ecg_fmt_info(char *buf, size_t cap, int32_t delay_samples)
{
    line_t l = ln_begin(buf, cap, "$I");
    ln_field(&l, ECG_FS_HZ, 0);
    ln_field(&l, ECG_PROTO_VER, 0);
    ln_field(&l, delay_samples, 0);
    return ln_end(&l);
}

size_t ecg_fmt_sample(char *buf, size_t cap, uint32_t n, int32_t uv)
{
    char tmp[16];
    line_t l = ln_begin(buf, cap, "$S");

    snprintf(tmp, sizeof(tmp), ",%lu", (unsigned long)n);
    ln_put(&l, tmp);
    ln_field(&l, uv, 0);
    return ln_end(&l);
}

size_t ecg_fmt_beat(char *buf, size_t cap, const ecg_beat_t *b)
{
    char tmp[16];
    line_t l = ln_begin(buf, cap, "$B");

    snprintf(tmp, sizeof(tmp), ",%lu", (unsigned long)b->r_index);
    ln_put(&l, tmp);
    ln_field(&l, b->hr_x10, 1);
    ln_field(&l, b->rr_ms, 0);
    ln_field(&l, b->qrs_ms, 0);
    ln_field(&l, b->pr_ms, 0);
    ln_field(&l, b->qt_ms, 0);
    ln_field(&l, b->r_uv, 3);
    ln_field(&l, b->st_uv, 3);
    return ln_end(&l);
}

size_t ecg_fmt_event(char *buf, size_t cap, uint32_t n, const char *code, int32_t value)
{
    char tmp[16];
    line_t l = ln_begin(buf, cap, "$E");

    snprintf(tmp, sizeof(tmp), ",%lu,", (unsigned long)n);
    ln_put(&l, tmp);
    ln_put(&l, code);
    if (value != ECG_NA)
        ln_field(&l, value, 0);
    return ln_end(&l);
}
