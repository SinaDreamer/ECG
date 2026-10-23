#ifndef ECG_CONFIG_H
#define ECG_CONFIG_H

/* Sampling */
#define ECG_FS_HZ               250
#define ECG_MS(ms)              ((((ms) * ECG_FS_HZ) + 500) / 1000)

/* Signal conditioning */
#ifndef ECG_HPF_HZ
#define ECG_HPF_HZ              0.5f    /* baseline wander */
#endif
/*
 * 1: linear-phase HPF (delay minus two cascaded moving averages). No ST/T
 *    distortion, adds a constant 0.57 / ECG_HPF_HZ s (1.14 s) delay.
 * 0: 2nd-order Butterworth IIR. Minimal delay, but its phase response
 *    depresses the ST segment by ~0.05 mV and shortens QT.
 */
#ifndef ECG_HPF_LINEAR_PHASE
#define ECG_HPF_LINEAR_PHASE    1
#endif
#define ECG_HPF_MA_MAX          512     /* moving-average buffer, fc >= 0.28 Hz */
#ifndef ECG_LPF_HZ
#define ECG_LPF_HZ              40.0f   /* EMG / muscle noise, 4th-order Butterworth */
#endif
#ifndef ECG_NOTCH_HZ
#define ECG_NOTCH_HZ            60.0f   /* mains; 50.0f outside the Americas */
#endif
#ifndef ECG_NOTCH_Q
#define ECG_NOTCH_Q             30.0f   /* -3 dB bandwidth = f0 / Q = 2 Hz */
#endif

/* Pan-Tompkins QRS detector */
#ifndef QRS_BP_LO_HZ
#define QRS_BP_LO_HZ            5.0f
#endif
#ifndef QRS_BP_HI_HZ
#define QRS_BP_HI_HZ            15.0f
#endif
#ifndef QRS_MWI_MS
#define QRS_MWI_MS              150     /* moving-window integrator width */
#endif
#ifndef QRS_REFRACTORY_MS
#define QRS_REFRACTORY_MS       200
#endif
#ifndef QRS_TWAVE_MS
#define QRS_TWAVE_MS            360     /* T-wave slope check window */
#endif
#ifndef QRS_LEARN_MS
#define QRS_LEARN_MS            2000    /* threshold training period */
#endif
#ifndef QRS_PEAK_CONFIRM_MS
#define QRS_PEAK_CONFIRM_MS     90      /* MWI peak hold before classification */
#endif
#ifndef QRS_SEARCHBACK_RR
#define QRS_SEARCHBACK_RR       1.66f   /* searchback after 166% of mean RR */
#endif
#ifndef QRS_SPKI_DECAY_MS
#define QRS_SPKI_DECAY_MS       2500    /* halve SPKI after this long without a beat */
#endif

/* Delineation */
#ifndef ECG_QRS_SLOPE_FRAC
#define ECG_QRS_SLOPE_FRAC      0.04f   /* QRS on/off: |slope| below this fraction of max */
#endif
#ifndef ECG_ST_OFFSET_MS
#define ECG_ST_OFFSET_MS        60      /* ST level measured at J + 60 ms ... */
#endif
#ifndef ECG_ST_OFFSET_FAST_MS
#define ECG_ST_OFFSET_FAST_MS   40      /* ... or J + 40 ms when HR > 120 bpm */
#endif
#ifndef ECG_PR_SEARCH_MS
#define ECG_PR_SEARCH_MS        300     /* P wave search window before QRS onset */
#endif
#ifndef ECG_T_SEARCH_MAX_MS
#define ECG_T_SEARCH_MAX_MS     600     /* T-end search limit after R */
#endif
#ifndef ECG_T_SEARCH_RR_FRAC
#define ECG_T_SEARCH_RR_FRAC    0.70f   /* ... or this fraction of RR, if shorter */
#endif
#ifndef ECG_P_MIN_MV
#define ECG_P_MIN_MV            0.03f
#endif
#ifndef ECG_T_MIN_MV
#define ECG_T_MIN_MV            0.05f
#endif

/* Filtered-signal history, must be a power of two (4.1 s) */
#define ECG_BUF_LEN             1024

#endif
