#ifndef ECG_DELINEATE_H
#define ECG_DELINEATE_H

/* Fiducial points of one beat, as indices into the analysed segment. */
typedef struct {
    int qrs_on;
    int qrs_off;        /* J point */
    float p_on;         /* < 0 if no P wave found */
    float t_end;        /* < 0 if no T wave found */
    float baseline;     /* isoelectric level (PR segment), mV */
    float r_amp;        /* R peak above baseline, mV */
    float st;           /* ST level at J + ECG_ST_OFFSET_MS, mV */
    int st_valid;
} ecg_delin_t;

/*
 * x:      filtered ECG segment in mV
 * r:      index of the dominant QRS deflection
 * p_min:  earliest index allowed for the P wave (end of previous T)
 * t_max:  latest index allowed for the T wave end (next QRS / RR limit)
 * st_off: ST measurement point, samples after the J point
 * Requires 2 <= r and t_max <= len - 3.
 */
void ecg_delineate(const float *x, int len, int r, int p_min, int t_max,
                   int st_off, ecg_delin_t *out);

#endif
