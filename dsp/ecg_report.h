#ifndef ECG_REPORT_H
#define ECG_REPORT_H

#include <stddef.h>
#include <stdint.h>
#include "ecg_monitor.h"

/*
 * Line protocol, ASCII, one record per line:
 *
 *   $I,<fs_hz>,<proto_ver>,<delay>*CS       boot banner; delay = filter lag, samples
 *   $S,<n>,<ecg_uv>*CS                      filtered ECG sample
 *   $B,<n>,<hr>,<rr>,<qrs>,<pr>,<qt>,<r>,<st>*CS
 *        n   R peak sample index     hr  bpm, 1 decimal
 *        rr  qrs pr qt  ms            r st  mV, 3 decimals
 *   $E,<n>,<code>[,<value>]*CS              event (LOFF, LON, OVR, DROP, TXD)
 *
 * Unavailable fields are left empty. CS is the XOR of all bytes between
 * '$' and '*', as two hex digits. Lines end with "\r\n".
 */

#define ECG_PROTO_VER   1
#define ECG_LINE_MAX    96

size_t ecg_fmt_info(char *buf, size_t cap, int32_t delay_samples);
size_t ecg_fmt_sample(char *buf, size_t cap, uint32_t n, int32_t uv);
size_t ecg_fmt_beat(char *buf, size_t cap, const ecg_beat_t *b);
size_t ecg_fmt_event(char *buf, size_t cap, uint32_t n, const char *code, int32_t value);

#endif
