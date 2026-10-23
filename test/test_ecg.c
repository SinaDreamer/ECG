/*
 * Host-side verification of the ECG DSP chain.
 *
 *   make -C test && ./test/test_ecg [dump.csv]
 *
 * Synthetic lead-II ECG is built from Gaussian waves with known fiducial
 * points, then corrupted with baseline wander, 60 Hz mains, EMG noise, a
 * front-end DC offset and 12-bit ADC quantisation.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ecg_config.h"
#include "ecg_filter.h"
#include "ecg_monitor.h"
#include "ecg_report.h"

#define PI          3.14159265358979
#define FS          ((double)ECG_FS_HZ)
#define MAX_BEATS   1024
#define GAIN        1100.0      /* front-end gain, electrode mV -> ADC mV */
#define ADC_LSB_MV  (3100.0 / 4095.0)

static int g_fail;

#define CHECK(cond, ...) do { \
    if (!(cond)) { g_fail++; printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

/* ---------------------------------------------------------------- filters */

static double bq_mag(const biquad_t *f, double hz)
{
    double w = 2.0 * PI * hz / FS;
    double nr = f->b0 + f->b1 * cos(w) + f->b2 * cos(2 * w);
    double ni = -f->b1 * sin(w) - f->b2 * sin(2 * w);
    double dr = 1.0 + f->a1 * cos(w) + f->a2 * cos(2 * w);
    double di = -f->a1 * sin(w) - f->a2 * sin(2 * w);
    return sqrt((nr * nr + ni * ni) / (dr * dr + di * di));
}

static double hpf_mag(const ecg_filter_t *f, double hz)
{
#if ECG_HPF_LINEAR_PHASE
    double w = 2.0 * PI * hz / FS, L = f->hp.len;
    double d = sin(w * L / 2.0) / (L * sin(w / 2.0));
    return fabs(1.0 - d * d);
#else
    return bq_mag(&f->hp, hz);
#endif
}

static double chain_db(const ecg_filter_t *f, double hz)
{
    double m = hpf_mag(f, hz) * bq_mag(&f->notch, hz) *
               bq_mag(&f->lp[0], hz) * bq_mag(&f->lp[1], hz);
    return 20.0 * log10(m);
}

static void test_filter(void)
{
    ecg_filter_t f;
    static const double pts[] = { 0.05, 0.1, 0.5, 1, 5, 10, 20, 40, 50, 59, 60, 61, 80, 100 };
    unsigned k;

    printf("filter response (dB)\n");
    ecg_filter_init(&f);
    printf("  HPF %s, delay %d samples (%.0f ms)\n",
           ECG_HPF_LINEAR_PHASE ? "linear-phase" : "IIR", ecg_filter_delay(&f),
           ecg_filter_delay(&f) * 1000.0 / FS);
    for (k = 0; k < sizeof(pts) / sizeof(pts[0]); ++k)
        printf("  %6.2f Hz  %8.2f\n", pts[k], chain_db(&f, pts[k]));

    CHECK(fabs(chain_db(&f, 0.5) + 3.0) < 0.5, "HPF corner not -3 dB");
    CHECK(fabs(chain_db(&f, 40.0) + 3.0) < 0.6, "LPF corner not -3 dB");
    CHECK(chain_db(&f, 0.1) < -25.0, "baseline wander not rejected");
    CHECK(chain_db(&f, 60.0) < -60.0, "60 Hz not rejected");
    CHECK(fabs(chain_db(&f, 10.0)) < 0.2, "passband not flat at 10 Hz");
    for (k = 100; k < 3000; ++k)
        if (chain_db(&f, k / 100.0) < -0.5) {
            CHECK(0, "passband ripple > 0.5 dB at %.2f Hz", k / 100.0);
            break;
        }

    /* time domain: DC offset gives no start-up transient, 60 Hz is removed */
    {
        double peak = 0;
        int i;
        ecg_filter_init(&f);
        for (i = 0; i < 10 * ECG_FS_HZ; ++i) {
            double x = 1.5 + 0.5 * sin(2 * PI * 60.0 * i / FS);
            double y = ecg_filter_run(&f, (float)x);
            if (i > 2 * ECG_FS_HZ && fabs(y) > peak)
                peak = fabs(y);
        }
        printf("  60 Hz 0.5 mV + 1.5 mV DC -> residual %.4f mV\n", peak);
        CHECK(peak < 0.005, "time-domain 60 Hz residual too high");
    }
}

/* -------------------------------------------------------------- synthetic */

typedef struct {
    double hr, rr_jitter_ms;
    double wander_mv, mains_mv, emg_mv;
    double st_mv;
    double seconds;
    unsigned seed;
} scenario_t;

typedef struct {
    double r_t;         /* seconds */
    double rr_ms, qrs_ms, pr_ms, qt_ms, r_mv, st_mv;
} truth_t;

typedef struct {
    double t, rr, pr_off;   /* beat time, RR (s), P centre offset */
    double t_mu, t_sd;
} beat_shape_t;

static double gauss(double t, double mu, double sd)
{
    double z = (t - mu) / sd;
    return exp(-0.5 * z * z);
}

static double randn(void)
{
    double u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

/* Wave set (ms relative to R):
 *   P  0.15 mV  mu -165  sd 18   onset (tangent, mu - 2sd)  -201
 *   Q -0.10     mu  -28  sd  6   QRS onset ~ mu - 2.5sd      -43
 *   R  1.20     mu    0  sd  9
 *   S -0.25     mu  +28  sd  7   QRS offset ~ mu + 2.5sd     +45
 *   T  0.30     mu/sd scale with sqrt(RR)  (Bazett, QTc 400 ms)
 */
#define QRS_ON_MS   (-43.0)
#define QRS_OFF_MS  (45.0)
#define P_MU_MS     (-165.0)
#define P_SD_MS     (18.0)

static double beat_value(const beat_shape_t *b, double t, double st_mv)
{
    double u = (t - b->t) * 1000.0;     /* ms from R */
    double v = 0.0;

    if (u < -400.0 || u > 800.0)
        return 0.0;
    v += 0.15 * gauss(u, P_MU_MS, P_SD_MS);
    v += -0.10 * gauss(u, -28.0, 6.0);
    v += 1.20 * gauss(u, 0.0, 9.0);
    v += -0.25 * gauss(u, 28.0, 7.0);
    v += 0.30 * gauss(u, b->t_mu, b->t_sd);
    if (st_mv != 0.0) {
        /* ST shift from the J point, merging into the T wave */
        double up = 0.5 * (1.0 + tanh((u - QRS_OFF_MS) / 4.0));
        double dn = 0.5 * (1.0 - tanh((u - b->t_mu) / 25.0));
        v += st_mv * up * dn;
    }
    return v;
}

static int synth(const scenario_t *sc, float **out, int *nout, truth_t *truth)
{
    static beat_shape_t beats[MAX_BEATS];
    int nb = 0, n = (int)(sc->seconds * FS), i, k;
    double t = 0.6, rr0 = 60.0 / sc->hr;
    float *x = malloc(sizeof(float) * (size_t)n);

    srand(sc->seed);
    while (t < sc->seconds + 1.0 && nb < MAX_BEATS) {
        double rr = rr0 + sc->rr_jitter_ms * 1e-3 * randn();
        double qt = 400.0 * sqrt(rr);              /* ms, from QRS onset */
        double t_sd = 40.0 * sqrt(rr);
        beat_shape_t *b = &beats[nb];

        b->t = t;
        b->rr = rr;
        b->t_sd = t_sd;
        b->t_mu = QRS_ON_MS + qt - 2.0 * t_sd;   /* tangent T end = mu + 2 sd */

        truth[nb].r_t = t;
        truth[nb].rr_ms = nb ? (t - beats[nb - 1].t) * 1000.0 : 0.0;
        truth[nb].qrs_ms = QRS_OFF_MS - QRS_ON_MS;
        truth[nb].pr_ms = QRS_ON_MS - (P_MU_MS - 2.0 * P_SD_MS);
        truth[nb].qt_ms = qt;
        truth[nb].r_mv = 1.2;
        truth[nb].st_mv = sc->st_mv;
        nb++;
        t += rr;
    }

    for (i = 0; i < n; ++i) {
        double ti = i / FS, v = 0.0, pin;
        for (k = 0; k < nb; ++k)
            if (fabs(beats[k].t - ti) < 1.2)
                v += beat_value(&beats[k], ti, sc->st_mv);
        v += sc->wander_mv * (sin(2 * PI * 0.25 * ti) + 0.5 * sin(2 * PI * 0.07 * ti + 1.0));
        v += sc->mains_mv * sin(2 * PI * 60.0 * ti);
        v += sc->emg_mv * randn();
        /* AD8232-style front end: gain, mid-rail offset, 12-bit ADC */
        pin = v * GAIN + 1650.0;
        pin = floor(pin / ADC_LSB_MV + 0.5) * ADC_LSB_MV;
        x[i] = (float)(pin / GAIN);
    }
    *out = x;
    *nout = n;
    return nb;
}

/* ------------------------------------------------------------- pipeline */

typedef struct {
    int tp, fn, fp;
    double e_rr, e_qrs, e_pr, e_qt, e_r, e_st;    /* mean abs error */
    double m_qrs, m_pr, m_qt, m_r, m_st;          /* mean measured */
    int n_rr, n_qrs, n_pr, n_qt, n_st;
} stats_t;

static void run_scenario(const char *name, const scenario_t *sc, FILE *dump,
                         double tol_pr, double tol_qt, double tol_st)
{
    static truth_t truth[MAX_BEATS];
    static ecg_monitor_t mon;
    static int matched[MAX_BEATS];
    float *x;
    int n, nb, i, k, skip_until, delay;
    stats_t s;
    char line[ECG_LINE_MAX];

    memset(&s, 0, sizeof(s));
    memset(matched, 0, sizeof(matched));
    nb = synth(sc, &x, &n, truth);
    ecg_monitor_init(&mon, 0);
    delay = ecg_filter_delay(&mon.filt);
    /* learning period plus HPF settling */
    skip_until = (int)(3.0 * FS) + delay;

    for (i = 0; i < n; ++i) {
        ecg_beat_t b;
        float y;
        if (ecg_monitor_process(&mon, x[i], &y, &b)) {
            double bt = ((double)b.r_index - delay) / FS;
            int best = -1;
            double bd = 1e9;

            if (dump) {
                ecg_fmt_beat(line, sizeof(line), &b);
                fputs(line, stdout);
            }
            if (b.r_index < (uint32_t)skip_until)
                continue;
            for (k = 0; k < nb; ++k)
                if (fabs(truth[k].r_t - bt) < bd) {
                    bd = fabs(truth[k].r_t - bt);
                    best = k;
                }
            if (best < 0 || bd > 0.05 || matched[best]) {
                s.fp++;
                continue;
            }
            matched[best] = 1;
            s.tp++;
            if (b.rr_ms != ECG_NA && best > 0 && matched[best - 1]) {
                s.e_rr += fabs(b.rr_ms - truth[best].rr_ms);
                s.n_rr++;
            }
            if (b.qrs_ms != ECG_NA) {
                s.e_qrs += fabs(b.qrs_ms - truth[best].qrs_ms);
                s.m_qrs += b.qrs_ms;
                s.n_qrs++;
            }
            if (b.pr_ms != ECG_NA) {
                s.e_pr += fabs(b.pr_ms - truth[best].pr_ms);
                s.m_pr += b.pr_ms;
                s.n_pr++;
            }
            if (b.qt_ms != ECG_NA) {
                s.e_qt += fabs(b.qt_ms - truth[best].qt_ms);
                s.m_qt += b.qt_ms;
                s.n_qt++;
            }
            if (b.st_uv != ECG_NA) {
                s.e_st += fabs(b.st_uv * 1e-3 - truth[best].st_mv);
                s.m_st += b.st_uv * 1e-3;
                s.n_st++;
            }
            s.e_r += fabs(b.r_uv * 1e-3 - truth[best].r_mv);
            s.m_r += b.r_uv * 1e-3;
        }
        if (dump && i < 8 * ECG_FS_HZ)
            fprintf(dump, "%s,%d,%.5f,%.5f\n", name, i, x[i], y);
    }
    /* beats still awaiting their T wave at the end of the record don't count */
    for (k = 0; k < nb; ++k)
        if (!matched[k] && truth[k].r_t * FS >= skip_until - delay &&
            truth[k].r_t + delay / FS < sc->seconds - 1.5)
            s.fn++;

#define AVG(sum, cnt) ((cnt) ? (sum) / (cnt) : NAN)
    printf("\n%s  (HR %.0f, %d truth beats)\n", name, sc->hr, nb);
    printf("  Se %.2f%%  PPV %.2f%%  (TP %d FN %d FP %d)\n",
           100.0 * s.tp / (s.tp + s.fn), 100.0 * s.tp / (s.tp + s.fp), s.tp, s.fn, s.fp);
    printf("  RR  MAE %5.1f ms\n", AVG(s.e_rr, s.n_rr));
    printf("  QRS mean %5.1f ms  truth %5.1f  MAE %5.1f  (%d/%d)\n",
           AVG(s.m_qrs, s.n_qrs), truth[1].qrs_ms, AVG(s.e_qrs, s.n_qrs), s.n_qrs, s.tp);
    printf("  PR  mean %5.1f ms  truth %5.1f  MAE %5.1f  (%d/%d)\n",
           AVG(s.m_pr, s.n_pr), truth[1].pr_ms, AVG(s.e_pr, s.n_pr), s.n_pr, s.tp);
    printf("  QT  mean %5.1f ms  truth %5.1f  MAE %5.1f  (%d/%d)\n",
           AVG(s.m_qt, s.n_qt), truth[1].qt_ms, AVG(s.e_qt, s.n_qt), s.n_qt, s.tp);
    printf("  R   mean %5.3f mV  truth %5.3f  MAE %5.3f\n",
           AVG(s.m_r, s.tp), truth[1].r_mv, AVG(s.e_r, s.tp));
    printf("  ST  mean %6.3f mV truth %6.3f  MAE %5.3f  (%d/%d)\n",
           AVG(s.m_st, s.n_st), truth[1].st_mv, AVG(s.e_st, s.n_st), s.n_st, s.tp);

    CHECK(s.fn == 0 && s.fp == 0, "%s: detection errors", name);
    CHECK(AVG(s.e_rr, s.n_rr) < 4.0, "%s: RR error", name);
    CHECK(s.n_qrs > 0.95 * s.tp && AVG(s.e_qrs, s.n_qrs) < 20.0, "%s: QRS duration", name);
    CHECK(AVG(s.e_r, s.tp) < 0.12, "%s: R amplitude", name);
    if (!ECG_HPF_LINEAR_PHASE)
        tol_pr = tol_qt = tol_st = 0;
    if (tol_pr > 0)
        CHECK(s.n_pr > 0.9 * s.tp && AVG(s.e_pr, s.n_pr) < tol_pr, "%s: PR interval", name);
    if (tol_qt > 0)
        CHECK(s.n_qt > 0.9 * s.tp && AVG(s.e_qt, s.n_qt) < tol_qt, "%s: QT interval", name);
    if (tol_st > 0)
        CHECK(s.n_st > 0.95 * s.tp && AVG(s.e_st, s.n_st) < tol_st, "%s: ST level", name);
    free(x);
}

/* ------------------------------------------------------------- protocol */

static int parse_ok(const char *line)
{
    const char *star = strchr(line, '*');
    unsigned cs = 0, want;
    const char *p;

    if (line[0] != '$' || !star)
        return 0;
    for (p = line + 1; p < star; ++p)
        cs ^= (unsigned char)*p;
    return sscanf(star + 1, "%2X", &want) == 1 && want == cs;
}

static void test_report(void)
{
    char line[ECG_LINE_MAX];
    ecg_beat_t b = { 123456, 723, 830, 92, 156, 384, 1234, -50 };
    ecg_beat_t na = { 7, ECG_NA, ECG_NA, 88, ECG_NA, ECG_NA, 980, ECG_NA };

    printf("\nprotocol\n");
    ecg_fmt_info(line, sizeof(line), 285);
    printf("  %s", line);
    CHECK(parse_ok(line), "info checksum");
    ecg_fmt_sample(line, sizeof(line), 4321, -1234);
    printf("  %s", line);
    CHECK(parse_ok(line) && strncmp(line, "$S,4321,-1234*", 14) == 0, "sample line");
    ecg_fmt_beat(line, sizeof(line), &b);
    printf("  %s", line);
    CHECK(parse_ok(line) &&
          strncmp(line, "$B,123456,72.3,830,92,156,384,1.234,-0.050*", 43) == 0, "beat line");
    ecg_fmt_beat(line, sizeof(line), &na);
    printf("  %s", line);
    CHECK(parse_ok(line) && strncmp(line, "$B,7,,,88,,,0.980,*", 19) == 0, "beat line NA");
    ecg_fmt_event(line, sizeof(line), 99, "LOFF", ECG_NA);
    printf("  %s", line);
    CHECK(parse_ok(line), "event line");
    CHECK(ecg_fmt_beat(line, 10, &b) == 0, "overflow not detected");
}

int main(int argc, char **argv)
{
    FILE *dump = argc > 1 ? fopen(argv[1], "w") : NULL;
    scenario_t clean   = { 72, 20, 0.0, 0.0, 0.0, 0.0, 60, 1 };
    scenario_t noisy   = { 72, 40, 0.5, 0.3, 0.02, 0.0, 120, 2 };
    scenario_t brady   = { 45, 30, 0.5, 0.3, 0.02, 0.0, 120, 3 };
    scenario_t tachy   = { 150, 10, 0.5, 0.3, 0.02, 0.0, 120, 4 };
    scenario_t emg     = { 80, 40, 0.3, 0.3, 0.05, 0.0, 120, 5 };
    scenario_t st_up   = { 72, 20, 0.3, 0.3, 0.02, 0.20, 60, 6 };
    scenario_t st_down = { 72, 20, 0.3, 0.3, 0.02, -0.15, 60, 7 };
    scenario_t irreg   = { 90, 150, 0.4, 0.3, 0.03, 0.0, 120, 8 };

    if (dump)
        fprintf(dump, "scenario,n,raw_mv,filt_mv\n");
    if (!ECG_HPF_LINEAR_PHASE)
        printf("IIR HPF build: PR/QT/ST reported but not checked (phase distortion)\n");

    test_filter();
    run_scenario("clean", &clean, dump, 15, 25, 0.03);
    run_scenario("noisy", &noisy, dump, 20, 30, 0.05);
    run_scenario("brady", &brady, dump, 20, 30, 0.05);
    run_scenario("tachy", &tachy, dump, 0, 40, 0);
    run_scenario("emg", &emg, dump, 25, 35, 0);
    run_scenario("st_elev", &st_up, dump, 0, 30, 0.03);
    run_scenario("st_depr", &st_down, dump, 0, 30, 0.03);
    run_scenario("irregular", &irreg, dump, 0, 40, 0.05);
    test_report();

    if (dump)
        fclose(dump);
    printf("\n%s (%d failure%s)\n", g_fail ? "FAILED" : "PASSED", g_fail, g_fail == 1 ? "" : "s");
    return g_fail ? 1 : 0;
}
