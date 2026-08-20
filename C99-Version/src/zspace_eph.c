/*
 * zspace_eph.c  ·  Harmonics / Alias Ephemeris Resolution (C99)
 * =============================================================
 * Direct port of zspace_engine/ephemeris.py (EphemerisResolver).
 *
 * Resolution policy: a candidate with exactly ONE significant dip in its own
 * fold is classified via the 2P/3P fold multiplicities:
 *     N3 == 1               -> P_TRUE/3  (physical = 3*P)
 *     N3 == 3 and N2 == 1   -> P_TRUE/2  (physical = 2*P)
 *     N3 == 3 and N2 == 2   -> FUNDAMENTAL
 *     anything else         -> FUNDAMENTAL (no resolution)
 * Over-harmonics (P_best == n*P_true, n = 2,3) are resolved downward only
 * when folding at P/n restores exactly one significant dip.
 */
#include "zspace_eph.h"

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DIP_MIN_SNR    5.0
#define REL_MIN_FRAC   0.25
#define ABS_SNR_MAD    3.0
#define INCLUDE_SNR    1.0
#define DECIDE_MIN_SNR 4.0
#define EQUAL_DEPTH_RATIO 1.7
#define DIP_N_BINS     200
#define FINE_N_BINS    600
#define DOWN_N_BINS    400
#define CENTER_N_BINS  400

/* ── Small helpers ─────────────────────────────────────────────────────────── */

static int cmp_dbl(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

static double median_dbl(const double *v, size_t n) {
    double *tmp = (double *)malloc(n * sizeof(double));
    if (!tmp) return 0.0;
    memcpy(tmp, v, n * sizeof(double));
    qsort(tmp, n, sizeof(double), cmp_dbl);
    double med = (n % 2) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
    free(tmp);
    return med;
}

/* phase in [-0.5, 0.5] centred on transit at 0 (mirrors phase_fold) */
static double phase_fold(double time, double period, double t0) {
    double ph = fmod((time - t0) / period, 1.0);
    if (ph < 0.0) ph += 1.0;
    if (ph > 0.5) ph -= 1.0;
    return ph;
}

/* ── fold_and_bin (mirrors BLSDetector.fold_and_bin) ───────────────────────── */

int zs_eph_fold_and_bin(const double *time, size_t n_time,
                        const double *flux,
                        double period, double t0, int n_bins,
                        double *bin_phase, double *bin_flux,
                        unsigned char *valid) {
    if (!time || !flux || !bin_phase || !bin_flux || !valid || n_bins < 2) return -1;
    double *sums = (double *)calloc((size_t)n_bins, sizeof(double));
    int *cnts = (int *)calloc((size_t)n_bins, sizeof(int));
    if (!sums || !cnts) { free(sums); free(cnts); return -1; }

    for (size_t i = 0; i < n_time; i++) {
        double ph = phase_fold(time[i], period, t0);
        int b = (int)floor((ph + 0.5) * (double)n_bins);
        if (b >= n_bins) b = n_bins - 1;
        if (b < 0) b = 0;
        sums[b] += flux[i];
        cnts[b]++;
    }
    for (int b = 0; b < n_bins; b++) {
        bin_phase[b] = (double)b / (double)n_bins - 0.5 + 0.5 / (double)n_bins;
        if (cnts[b] >= 3) {
            bin_flux[b] = sums[b] / (double)cnts[b];
            valid[b] = 1;
        } else {
            bin_flux[b] = 0.0;
            valid[b] = 0;
        }
    }
    free(sums); free(cnts);
    return 0;
}

/* ── dip signature (mirrors _fold_dip_signature) ───────────────────────────── */

/* forward decl (defined below, part 2) */
int zs_eph_dip_signature_groups(const double *bf, const double *ph,
                                const unsigned char *valid, int n_bins,
                                const int *below, int n_below,
                                double baseline, double per_bin, double bin_w,
                                ZSEphSignature *sig, int n_valid);

int zs_eph_dip_signature(const double *time, size_t n_time,
                         const double *flux,
                         double period, double t0, int n_bins,
                         ZSEphSignature *sig) {
    if (!sig) return -1;
    memset(sig, 0, sizeof(*sig));
    sig->period = period;
    sig->n_bins = n_bins;
    sig->dips = NULL;
    sig->n_dips = 0;

    double *bf = (double *)malloc((size_t)n_bins * sizeof(double));
    double *ph = (double *)malloc((size_t)n_bins * sizeof(double));
    unsigned char *valid = (unsigned char *)malloc((size_t)n_bins);
    if (!bf || !ph || !valid) { free(bf); free(ph); free(valid); return -1; }

    int rc = zs_eph_fold_and_bin(time, n_time, flux, period, t0, n_bins, ph, bf, valid);
    if (rc != 0) { free(bf); free(ph); free(valid); return rc; }

    int n_valid = 0;
    for (int b = 0; b < n_bins; b++) if (valid[b]) n_valid++;
    sig->covered_fraction = (double)n_valid / (double)n_bins;
    if (n_valid == 0) { free(bf); free(ph); free(valid); return 0; }

    /* baseline = nanmedian(bf); per_bin = 1.4826 * median(|bf - baseline|) */
    double *vals = (double *)malloc((size_t)n_valid * sizeof(double));
    double *devs = (double *)malloc((size_t)n_valid * sizeof(double));
    double *absdevs = (double *)malloc((size_t)n_valid * sizeof(double));
    if (!vals || !devs || !absdevs) {
        free(vals); free(devs); free(absdevs); free(bf); free(ph); free(valid);
        return -1;
    }
    int v = 0;
    for (int b = 0; b < n_bins; b++) if (valid[b]) vals[v++] = bf[b];
    double baseline = median_dbl(vals, (size_t)n_valid);
    for (int b = 0, k = 0; b < n_bins; b++) {
        if (valid[b]) devs[k++] = bf[b] - baseline;
    }
    for (int k = 0; k < n_valid; k++) absdevs[k] = fabs(devs[k]);
    double per_bin = 1.4826 * median_dbl(absdevs, (size_t)n_valid);
    if (!(isfinite(per_bin) && per_bin > 1e-15)) {
        double s = 0.0;
        for (int k = 0; k < n_valid; k++) s += devs[k] * devs[k];
        per_bin = sqrt(s / (double)n_valid);
        if (!(isfinite(per_bin) && per_bin > 1e-15)) per_bin = 1.0;
    }
    sig->baseline = baseline;
    sig->per_bin_noise = per_bin;

    /* Pass 1: bins >= 3 MAD below the robust median */
    double thr = baseline - ABS_SNR_MAD * per_bin;
    int *below = (int *)malloc((size_t)n_bins * sizeof(int));
    if (!below) { free(below); free(absdevs); free(vals); free(devs); free(bf); free(ph); free(valid); return -1; }
    int n_below = 0;
    for (int b = 0; b < n_bins; b++)
        if (valid[b] && bf[b] < thr) below[n_below++] = b;
    double bin_w = 1.0 / (double)n_bins;

    if (n_below == 0) {
        free(below); free(absdevs); free(vals); free(devs); free(bf); free(ph); free(valid);
        return 0;
    }
    free(absdevs); free(vals); free(devs);
    return zs_eph_dip_signature_groups(bf, ph, valid, n_bins, below, n_below,
                                       baseline, per_bin, bin_w, sig, n_valid);
}

#include "zspace_eph.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int zs_eph_dip_signature_groups(const double *bf, const double *ph,
                                const unsigned char *valid, int n_bins,
                                const int *below, int n_below,
                                double baseline, double per_bin, double bin_w,
                                ZSEphSignature *sig, int n_valid);

int zs_eph_dip_signature_groups(const double *bf, const double *ph,
                                const unsigned char *valid, int n_bins,
                                const int *below, int n_below,
                                double baseline, double per_bin, double bin_w,
                                ZSEphSignature *sig, int n_valid) {
    int n_groups = 0;
    for (int j = 0; j < n_below; j++) {
        if (j == 0 || below[j] - below[j - 1] > 2) n_groups++;
    }
    if (n_groups > n_bins) n_groups = n_bins;
    unsigned char *member = (unsigned char *)calloc((size_t)n_groups * n_bins, 1);
    int *gcount = (int *)calloc((size_t)n_groups, sizeof(int));
    if (!member || !gcount) { free(member); free(gcount); return -1; }

    int g = 0;
    member[g * n_bins + below[0]] = 1;
    gcount[g] = 1;
    for (int j = 1; j < n_below; j++) {
        if (below[j] - below[j - 1] > 2) g++;
        member[g * n_bins + below[j]] = 1;
        gcount[g]++;
    }

    double include_thr = baseline - 1.0 * per_bin;
    int grown = 1;
    while (grown) {
        grown = 0;
        for (int gg = 0; gg < n_groups; gg++) {
            for (int b = 0; b < n_bins; b++) {
                if (!member[gg * n_bins + b]) continue;
                int nb1 = (b + 1) % n_bins;
                int nb2 = (b - 1 + n_bins) % n_bins;
                for (int s = 0; s < 2; s++) {
                    int nb = (s == 0) ? nb1 : nb2;
                    if (valid[nb] && bf[nb] <= include_thr && !member[gg * n_bins + nb]) {
                        member[gg * n_bins + nb] = 1;
                        gcount[gg]++;
                        grown = 1;
                    }
                }
            }
        }
    }
    (void)n_valid;

    int cap = n_groups;
    ZSEphDip *dips = (ZSEphDip *)calloc((size_t)(cap > 0 ? cap : 1), sizeof(ZSEphDip));
    if (!dips) { free(member); free(gcount); return -1; }
    int n_dips = 0;
    for (int gg = 0; gg < n_groups; gg++) {
        if (gcount[gg] == 0) continue;
        double wsum = 0.0, csum = 0.0, dsum = 0.0;
        double gmin = 1e9, gmax = -1e9;
        for (int b = 0; b < n_bins; b++) {
            if (!member[gg * n_bins + b]) continue;
            double w = (baseline - bf[b]) > 0.0 ? (baseline - bf[b]) : 0.0;
            w += 1e-12;
            csum += ph[b] * w;
            wsum += w;
            dsum += baseline - bf[b];
            if (ph[b] < gmin) gmin = ph[b];
            if (ph[b] > gmax) gmax = ph[b];
        }
        if (wsum <= 0.0) continue;
        double center = csum / wsum;
        double depth = dsum / (double)gcount[gg];
        if (depth < 0.0) depth = 0.0;
        double lo = gmin - 0.5 * bin_w;
        double hi = gmax + 0.5 * bin_w;
        double depth_snr = depth * sqrt((double)gcount[gg]) / (per_bin > 1e-15 ? per_bin : 1e-15);
        dips[n_dips].phase_center = center;
        dips[n_dips].phase_lo = lo;
        dips[n_dips].phase_hi = hi;
        dips[n_dips].depth = depth;
        dips[n_dips].depth_snr = depth_snr;
        dips[n_dips].width_phase = hi - lo;
        n_dips++;
    }

    double max_depth = 0.0;
    for (int i = 0; i < n_dips; i++)
        if (dips[i].depth > max_depth) max_depth = dips[i].depth;
    ZSEphDip *confirmed = (ZSEphDip *)malloc((size_t)(n_dips > 0 ? n_dips : 1) * sizeof(ZSEphDip));
    if (!confirmed) { free(member); free(gcount); free(dips); return -1; }
    int n_conf = 0;
    for (int i = 0; i < n_dips; i++) {
        if (dips[i].depth_snr < DIP_MIN_SNR) continue;
        if (max_depth > 1e-15 && dips[i].depth < REL_MIN_FRAC * max_depth) continue;
        confirmed[n_conf++] = dips[i];
    }
    for (int i = 1; i < n_conf; i++) {
        ZSEphDip key = confirmed[i];
        int j = i - 1;
        while (j >= 0 && confirmed[j].phase_center > key.phase_center) {
            confirmed[j + 1] = confirmed[j];
            j--;
        }
        confirmed[j + 1] = key;
    }

    free(member); free(gcount); free(dips);
    sig->dips = confirmed;
    sig->n_dips = n_conf;
    return 0;
}/* part 3: resolver logic (mirrors EphemerisResolver) */
#include "zspace_eph.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void zs_eph_signature_free(ZSEphSignature *sig) {
    if (!sig) return;
    free(sig->dips);
    sig->dips = NULL;
    sig->n_dips = 0;
}

/* mirrors _dips_equal_depth */
static int dips_equal_depth(const ZSEphDip *dips, int n, double max_ratio) {
    if (n == 0) return 0;
    double dmax = 0.0, dmin = 1e18;
    for (int i = 0; i < n; i++) {
        double d = dips[i].depth > 0.0 ? dips[i].depth : 0.0;
        if (d > dmax) dmax = d;
        if (d < dmin) dmin = d;
    }
    return (dmin > 1e-12) && (dmax / dmin) <= max_ratio;
}

static double min_dip_snr(const ZSEphDip *dips, int n) {
    double m = 0.0;
    for (int i = 0; i < n; i++) {
        if (i == 0 || dips[i].depth_snr < m) m = dips[i].depth_snr;
    }
    return m;
}

/* mirrors _spacing_presence */
static int spacing_presence(const ZSEphDip *dips, int n, double expected) {
    if (n < 2) return 1;
    double *c = (double *)malloc((size_t)n * sizeof(double));
    if (!c) return 1;
    for (int i = 0; i < n; i++) c[i] = dips[i].phase_center;
    for (int i = 1; i < n; i++) {
        double key = c[i];
        int j = i - 1;
        while (j >= 0 && c[j] > key) { c[j + 1] = c[j]; j--; }
        c[j + 1] = key;
    }
    double tol = 0.5 * expected;
    int ok = 1;
    for (int i = 0; i < n - 1; i++) {
        double d = fabs(c[i + 1] - c[i]);
        d = fmin(d, 1.0 - d);
        if (fabs(d - expected) > tol) { ok = 0; break; }
    }
    free(c);
    return ok;
}

/* mirrors _merge_near_dips (keep deeper of near-identical events) */
static int merge_near_dips(const ZSEphDip *in, int n, ZSEphDip *out, double tol) {
    int n_out = 0;
    for (int i = 0; i < n; i++) {
        if (n_out == 0) {
            out[n_out++] = in[i];
            continue;
        }
        double sep = fmin(fabs(in[i].phase_center - out[n_out - 1].phase_center),
                          1.0 - fabs(in[i].phase_center - out[n_out - 1].phase_center));
        if (sep < tol) {
            if (in[i].depth > out[n_out - 1].depth) out[n_out - 1] = in[i];
        } else {
            out[n_out++] = in[i];
        }
    }
    /* sort by phase_center */
    for (int i = 1; i < n_out; i++) {
        ZSEphDip key = out[i];
        int j = i - 1;
        while (j >= 0 && out[j].phase_center > key.phase_center) {
            out[j + 1] = out[j];
            j--;
        }
        out[j + 1] = key;
    }
    return n_out;
}

/* mirrors _center_t0: shift t0 so deepest bin sits at phase 0 */
static double center_t0(const double *time, size_t n_time, const double *flux,
                        double period, double t0) {
    double *ph = (double *)malloc((size_t)CENTER_N_BINS * sizeof(double));
    double *bf = (double *)malloc((size_t)CENTER_N_BINS * sizeof(double));
    unsigned char *valid = (unsigned char *)malloc((size_t)CENTER_N_BINS);
    if (!ph || !bf || !valid) { free(ph); free(bf); free(valid); return t0; }
    if (zs_eph_fold_and_bin(time, n_time, flux, period, t0, CENTER_N_BINS, ph, bf, valid) != 0) {
        free(ph); free(bf); free(valid);
        return t0;
    }
    int imin = -1;
    double bmin = 1e300;
    for (int b = 0; b < CENTER_N_BINS; b++) {
        if (valid[b] && bf[b] < bmin) { bmin = bf[b]; imin = b; }
    }
    double res = t0;
    if (imin >= 0) res = t0 + ph[imin] * period;
    free(ph); free(bf); free(valid);
    return res;
}

/* forward decl (part 2) */
int zs_eph_dip_signature_groups(const double *bf, const double *ph,
                                const unsigned char *valid, int n_bins,
                                const int *below, int n_below,
                                double baseline, double per_bin, double bin_w,
                                ZSEphSignature *sig, int n_valid);

/* filled in part 3b */
static void fill_result(ZSEphResult *out, double period_best, double mult,
                        const char *classifier, const char *evidence,
                        double conf, const char *flag1, const char *flag2);

static int resolve_down(const double *time, size_t n_time, const double *flux,
                        double period_best, double t0, int n1,
                        const ZSEphDip *own_dips, int n_own,
                        double period_min, double *p_down_out, double *conf_out) {
    if (n1 != 2 && n1 != 3) return 0;
    if (!dips_equal_depth(own_dips, n_own, EQUAL_DEPTH_RATIO)) return 0;
    if (!spacing_presence(own_dips, n_own, 1.0 / (double)n1)) return 0;
    double p_down = period_best / (double)n1;
    if (p_down < period_min * 0.98) return 0;

    double t0d = center_t0(time, n_time, flux, p_down, t0);
    ZSEphSignature sig;
    if (zs_eph_dip_signature(time, n_time, flux, p_down, t0d, DOWN_N_BINS, &sig) != 0)
        return 0;
    ZSEphDip *merged = (ZSEphDip *)malloc((size_t)(sig.n_dips > 0 ? sig.n_dips : 1) * sizeof(ZSEphDip));
    if (!merged) { zs_eph_signature_free(&sig); return 0; }
    int n_merged = merge_near_dips(sig.dips, sig.n_dips, merged, 0.02);
    int ok = (n_merged == 1);
    double conf = min_dip_snr(sig.dips, sig.n_dips);
    free(merged);
    zs_eph_signature_free(&sig);
    if (!ok || conf < DECIDE_MIN_SNR) return 0;
    *p_down_out = p_down;
    *conf_out = conf;
    return 1;
}

int zs_eph_resolve(const double *time, size_t n_time,
                   const double *flux,
                   double period_best, double t0,
                   double period_min, double period_max,
                   ZSEphResult *out) {
    if (!out) return -1;
    memset(out, 0, sizeof(*out));
    out->candidate_period = period_best;
    out->physical_period = period_best;
    out->multiple = 1;
    snprintf(out->classifier, sizeof(out->classifier), "%s", "FUNDAMENTAL");

    ZSEphSignature sig_own;
    if (zs_eph_dip_signature(time, n_time, flux, period_best, t0, DIP_N_BINS, &sig_own) != 0)
        return -1;
    int n1 = sig_own.n_dips;
    snprintf(out->pattern, sizeof(out->pattern), "{'own':%d}", n1);

    /* over-harmonic probe with centred t0 + fine bins */
    double t0_fine = center_t0(time, n_time, flux, period_best, t0);
    ZSEphSignature sig_fine;
    if (zs_eph_dip_signature(time, n_time, flux, period_best, t0_fine, FINE_N_BINS, &sig_fine) == 0) {
        ZSEphDip *fine_merged = (ZSEphDip *)malloc(
            (size_t)(sig_fine.n_dips > 0 ? sig_fine.n_dips : 1) * sizeof(ZSEphDip));
        if (fine_merged) {
            int n_fine = merge_near_dips(sig_fine.dips, sig_fine.n_dips, fine_merged, 0.02);
            if (1 < n_fine && n_fine <= 3) {
                double p_down = 0.0, conf_down = 0.0;
                if (resolve_down(time, n_time, flux, period_best, t0, n_fine,
                                 fine_merged, n_fine, period_min, &p_down, &conf_down)) {
                    char ev[512];
                    snprintf(ev, sizeof(ev),
                             "own fold %d equally-spaced dips; P/%d fold restores 1 dip "
                             "(SNR=%.1f) -> P_true = P/%d",
                             n_fine, n_fine, conf_down, n_fine);
                    free(fine_merged);
                    zs_eph_signature_free(&sig_fine);
                    zs_eph_signature_free(&sig_own);
                    out->candidate_period = period_best;
                    out->physical_period = p_down;
                    out->multiple = n_fine;
                    snprintf(out->classifier, sizeof(out->classifier), "%s", "OVER_HARMONIC");
                    snprintf(out->evidence, sizeof(out->evidence), "%s", ev);
                    out->confidence = conf_down;
                    snprintf(out->pattern, sizeof(out->pattern),
                             "{'own':%d,'at_down_%dp':1}", n1, n_fine);
                    out->n_flags = 2;
                    snprintf(out->flags[0], sizeof(out->flags[0]), "%s", "EPHEMERIS_RESOLVED");
                    snprintf(out->flags[1], sizeof(out->flags[1]), "RESOLVE_DOWN_%d", n_fine);
                    return 0;
                }
            }
            free(fine_merged);
        }
        zs_eph_signature_free(&sig_fine);
    }

    if (n1 != 1) {
        char ev[512];
        snprintf(ev, sizeof(ev),
                 "own fold has %d dips (need exactly 1 to test sub-harmonics) -> no resolution", n1);
        zs_eph_signature_free(&sig_own);
        fill_result(out, period_best, 1.0, "FUNDAMENTAL", ev, 0.0, "", "");
        return 0;
    }

    /* count dips in the 2P and 3P folds (physical multiples inside window) */
    ZSEphSignature sig2, sig3;
    int n2 = -1, n3 = -1;
    char pat2[32] = "", pat3[32] = "";
    if (2.0 * period_best <= period_max * 1.02) {
        if (zs_eph_dip_signature(time, n_time, flux, 2.0 * period_best, t0, DIP_N_BINS, &sig2) == 0) {
            n2 = sig2.n_dips;
            snprintf(pat2, sizeof(pat2), "%d", n2);
        }
    }
    if (3.0 * period_best <= period_max * 1.02) {
        if (zs_eph_dip_signature(time, n_time, flux, 3.0 * period_best, t0, DIP_N_BINS, &sig3) == 0) {
            n3 = sig3.n_dips;
            snprintf(pat3, sizeof(pat3), "%d", n3);
        }
    }
    snprintf(out->pattern, sizeof(out->pattern), "{'own':%d,'at_2p':%s,'at_3p':%s}",
             n1, pat2[0] ? pat2 : "None", pat3[0] ? pat3 : "None");

    if (pat3[0] == '\0') {
        zs_eph_signature_free(&sig_own);
        fill_result(out, period_best, 1.0, "FUNDAMENTAL",
                    "3P outside search window -> cannot verify sub-harmonic -> fundamental",
                    0.0, "", "");
        return 0;
    }

    if (n3 == 1) {
        double conf = min_dip_snr(sig3.dips, sig3.n_dips);
        if (conf < DECIDE_MIN_SNR) {
            zs_eph_signature_free(&sig_own); zs_eph_signature_free(&sig3);
            char ev[512];
            snprintf(ev, sizeof(ev), "3P fold weak/ambiguous (N3=1, min dip SNR=%.1f < %.1f) -> no resolution (keep P)",
                     conf, DECIDE_MIN_SNR);
            fill_result(out, period_best, 1.0, "FUNDAMENTAL", ev, conf,
                        "LOW_CONFIDENCE_MULTIPLE_FOLD", "");
            return 0;
        }
        char ev[512];
        snprintf(ev, sizeof(ev),
                 "own fold 1 dip; 3P fold restores 1 dip (SNR=%.1f) -> P_true = 3*P", conf);
        zs_eph_signature_free(&sig_own); zs_eph_signature_free(&sig3);
        fill_result(out, period_best, 3.0, "P_TRUE/3_ALIAS", ev, conf,
                    "EPHEMERIS_RESOLVED", "RESOLVE_X3");
        return 0;
    }

    if (n3 == 3) {
        if (n2 == 1) {
            double conf = min_dip_snr(sig2.dips, sig2.n_dips);
            if (conf < DECIDE_MIN_SNR) {
                zs_eph_signature_free(&sig_own); zs_eph_signature_free(&sig3);
                char ev[512];
                snprintf(ev, sizeof(ev), "3P fold weak/ambiguous (N3=3,N2=1, min dip SNR=%.1f < %.1f) -> no resolution (keep P)",
                         conf, DECIDE_MIN_SNR);
                fill_result(out, period_best, 1.0, "FUNDAMENTAL", ev, conf,
                            "LOW_CONFIDENCE_MULTIPLE_FOLD", "");
                return 0;
            }
            char ev[512];
            snprintf(ev, sizeof(ev),
                     "own fold 1 dip; 3P fold shows 3 dips, 2P fold restores 1 dip (SNR=%.1f) -> P_true = 2*P", conf);
            zs_eph_signature_free(&sig_own); zs_eph_signature_free(&sig3);
            fill_result(out, period_best, 2.0, "P_TRUE/2_ALIAS", ev, conf,
                        "EPHEMERIS_RESOLVED", "RESOLVE_X2");
            return 0;
        }
        int equal3 = dips_equal_depth(sig3.dips, sig3.n_dips, EQUAL_DEPTH_RATIO) &&
                     spacing_presence(sig3.dips, sig3.n_dips, 1.0 / 3.0);
        double conf3 = min_dip_snr(sig3.dips, sig3.n_dips);
        char ev[512];
        snprintf(ev, sizeof(ev),
                 "own fold 1 dip; 3P fold shows 3 dip(s) %s; 2P fold shows %d dip(s) -> P_best is fundamental",
                 equal3 ? "equal+spaced" : "irregular", n2);
        zs_eph_signature_free(&sig_own); zs_eph_signature_free(&sig3);
        fill_result(out, period_best, 1.0, "FUNDAMENTAL", ev,
                    conf3, "", "");
        return 0;
    }

    {
        char ev[512];
        snprintf(ev, sizeof(ev), "3P fold weak/ambiguous (N3=%d) -> no resolution (keep P)", n3);
        zs_eph_signature_free(&sig_own); zs_eph_signature_free(&sig3);
        fill_result(out, period_best, 1.0, "FUNDAMENTAL", ev, 0.0,
                    "LOW_CONFIDENCE_MULTIPLE_FOLD", "");
        return 0;
    }
}

static void fill_result(ZSEphResult *out, double period_best, double mult,
                        const char *classifier, const char *evidence,
                        double conf, const char *flag1, const char *flag2) {
    out->candidate_period = period_best;
    out->physical_period = mult * period_best;
    out->multiple = (int)(mult + 0.5);
    snprintf(out->classifier, sizeof(out->classifier), "%s", classifier);
    snprintf(out->evidence, sizeof(out->evidence), "%s", evidence);
    out->confidence = conf;
    out->n_flags = 0;
    if (flag1 && flag1[0]) snprintf(out->flags[out->n_flags++], sizeof(out->flags[0]), "%s", flag1);
    if (flag2 && flag2[0]) snprintf(out->flags[out->n_flags++], sizeof(out->flags[1]), "%s", flag2);
}