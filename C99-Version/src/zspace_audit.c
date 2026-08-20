/*
 * zspace_audit.c  ·  Transit vitality audits (C99)
 * =================================================
 * Port of zspace_engine/auditors.py: extract_individual_transit_depths,
 * even_odd_test, depth_consistency_score, secondary_eclipse_test,
 * ingress_egress_test.
 */
#include "zspace_audit.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define EVEN_ODD_SIGMA_THRESHOLD 3.0
#define EVEN_ODD_P_VALUE_THRESHOLD 0.01
#define CHI2_NORMALISATION       4.0
#define MIN_TRANSITS_FOR_EOTEST  4
#define INGRESS_FRACTION_THRESHOLD 0.45
#define TEMPL_X_B                100

/* ── small helpers ─────────────────────────────────────────────────────────── */

static int cmp_dbl(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da < db) ? -1 : (da > db) ? 1 : 0;
}

static double median_dbl(const double *v, size_t n) {
    if (n == 0) return 0.0;
    double *tmp = (double *)malloc(n * sizeof(double));
    if (!tmp) return 0.0;
    memcpy(tmp, v, n * sizeof(double));
    qsort(tmp, n, sizeof(double), cmp_dbl);
    double med = (n % 2) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
    free(tmp);
    return med;
}

static double mean_dbl(const double *v, size_t n) {
    if (n == 0) return 0.0;
    double s = 0.0;
    for (size_t i = 0; i < n; i++) s += v[i];
    return s / (double)n;
}

static double std_dbl_ddof1(const double *v, size_t n) {
    if (n < 2) return 0.0;
    double m = mean_dbl(v, n), s = 0.0;
    for (size_t i = 0; i < n; i++) { double d = v[i] - m; s += d * d; }
    return sqrt(s / (double)(n - 1));
}

/* phase in [-0.5, 0.5] (mirrors BLSDetector.phase_fold) */
static double phase_fold(double time, double period, double t0) {
    double ph = fmod((time - t0) / period, 1.0);
    if (ph < 0.0) ph += 1.0;
    if (ph > 0.5) ph -= 1.0;
    return ph;
}

/* linear interpolation on sorted x (mirrors np.interp, 0 outside range) */
static double interp_lin(const double *xs, const double *ys, size_t n,
                         double x) {
    if (n == 0) return 0.0;
    if (x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (xs[mid] <= x) lo = mid; else hi = mid;
    }
    double t = (x - xs[lo]) / (xs[hi] - xs[lo]);
    return ys[lo] + t * (ys[hi] - ys[lo]);
}

/* ── incomplete beta / Student-t CDF (for Welch p-value) ───────────────────── */

static double betacf(double a, double b, double x) {
    const int MAXIT = 200;
    const double EPS = 3.0e-12, FPMIN = 1.0e-300;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap;
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= MAXIT; m++) {
        int m2 = 2 * m;
        double aa = (double)m * (b - (double)m) * x /
                    ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + (double)m) * (qab + (double)m) * x /
             ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }
    return h;
}

static double betai(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    double bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b) +
                    a * log(x) + b * log1p(-x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return bt * betacf(a, b, x) / a;
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

/* two-sided Student-t survival: 2 * P(T > |t|)  (scipy ttest_ind) */
static double ttest_twosided_p(double t, double df) {
    double x = df / (df + t * t);
    double p = betai(0.5 * df, 0.5, x);
    return 2.0 * p;
}

/* ── extract_individual_transit_depths ─────────────────────────────────────── */

int zs_extract_depths(const double *time, size_t n_time, const double *flux,
                      double period, double t0, double duration,
                      ZSExtractResult *out) {
    if (!time || !flux || !out) return -1;
    memset(out, 0, sizeof(*out));
    if (n_time == 0 || period <= 0.0 || duration <= 0.0) return 0;

    int n_min = (int)floor((time[0] - t0) / period);
    int n_max = (int)ceil((time[n_time - 1] - t0) / period) + 1;
    int n_epochs = n_max - n_min + 1;
    if (n_epochs < 1) return 0;

    /* cadence: median of positive diffs of sorted time */
    double *ts_sorted = (double *)malloc(n_time * sizeof(double));
    if (!ts_sorted) return -1;
    memcpy(ts_sorted, time, n_time * sizeof(double));
    qsort(ts_sorted, n_time, sizeof(double), cmp_dbl);
    double *dts = (double *)malloc(n_time * sizeof(double));
    if (!dts) { free(ts_sorted); return -1; }
    size_t n_dts = 0;
    for (size_t i = 1; i < n_time; i++) {
        double d = ts_sorted[i] - ts_sorted[i - 1];
        if (d > 0.0) dts[n_dts++] = d;
    }
    double cadence = (n_dts > 0) ? median_dbl(dts, n_dts) : duration;
    free(dts); free(ts_sorted);
    if (!(cadence > 0.0)) cadence = duration;
    if (cadence <= 1e-9) cadence = 1e-9;

    double half_dur  = duration / 2.0;
    double oot_half  = duration * 3.0;
    double temp_half = half_dur * 2.0;
    double temp_phase_h = temp_half / period;

    /* ── folded template on 100 bins across [-temp_phase_h, +temp_phase_h] ─── */
    double bin_w = 2.0 * temp_phase_h / (double)TEMPL_X_B;
    double t_dip[TEMPL_X_B], t_cnt[TEMPL_X_B];
    memset(t_dip, 0, sizeof(t_dip));
    memset(t_cnt, 0, sizeof(t_cnt));
    for (size_t i = 0; i < n_time; i++) {
        double ph = phase_fold(time[i], period, t0);
        if (fabs(ph) <= temp_phase_h) {
            int b = (int)floor((ph + temp_phase_h) / bin_w);
            if (b < 0) b = 0;
            if (b >= TEMPL_X_B) b = TEMPL_X_B - 1;
            t_dip[b] += flux[i];
            t_cnt[b] += 1.0;
        }
    }
    int n_good = 0;
    for (int b = 0; b < TEMPL_X_B; b++) {
        if (t_cnt[b] > 0.0) { t_dip[b] /= t_cnt[b]; n_good++; }
    }
    if ((double)n_good < TEMPL_X_B * 0.4) return 0;   /* fall back: empty */

    double t_xs[TEMPL_X_B];
    for (int b = 0; b < TEMPL_X_B; b++)
        t_xs[b] = -temp_phase_h + (double)b * bin_w + 0.5 * bin_w;
    /* interpolate empty bins from the good ones */
    int *gidx = (int *)malloc((size_t)n_good * sizeof(int));
    if (!gidx) return -1;
    int k = 0;
    for (int b = 0; b < TEMPL_X_B; b++) if (t_cnt[b] > 0.0) gidx[k++] = b;
    for (int b = 0; b < TEMPL_X_B; b++) {
        if (t_cnt[b] > 0.0) continue;
        double x = t_xs[b];
        size_t lo = 0, hi = (size_t)n_good - 1;
        while (hi - lo > 1) {
            size_t mid = (lo + hi) / 2;
            if (t_xs[gidx[mid]] <= x) lo = mid; else hi = mid;
        }
        double x0 = t_xs[gidx[lo]], x1 = t_xs[gidx[hi]];
        double tt = (x - x0) / (x1 - x0);
        t_dip[b] = t_dip[gidx[lo]] + tt * (t_dip[gidx[hi]] - t_dip[gidx[lo]]);
    }
    free(gidx);

    double depth_ref = 0.0;
    for (int b = 0; b < TEMPL_X_B; b++) depth_ref += 1.0 - t_dip[b];
    depth_ref /= (double)TEMPL_X_B;
    if (depth_ref <= 0.0) depth_ref = 1e-8;
    double t_dip_norm[TEMPL_X_B];
    for (int b = 0; b < TEMPL_X_B; b++)
        t_dip_norm[b] = (1.0 - t_dip[b]) / depth_ref;

    /* template_at(phase): clip(interp(phase, t_xs, t_dip_norm, 0, 0), 0, None) */
    #define TEMPLATE_AT(p) \
        (((p) < t_xs[0] || (p) > t_xs[TEMPL_X_B - 1]) ? 0.0 : \
         fmax(interp_lin(t_xs, t_dip_norm, TEMPL_X_B, (p)), 0.0))

    double shift_step = cadence;
    int n_shifts = (int)(2.0 * temp_half / shift_step + 1.0 + 1e-9);
    if (n_shifts < 1) n_shifts = 1;

    size_t cap = 16;
    double *depths = (double *)malloc(cap * sizeof(double));
    double *errs = (double *)malloc(cap * sizeof(double));
    int *ns = (int *)malloc(cap * sizeof(int));
    if (!depths || !errs || !ns) {
        free(depths); free(errs); free(ns);
        return -1;
    }
    int n_out = 0;

    for (int n = n_min; n <= n_max; n++) {
        double t_centre = t0 + (double)n * period;

        /* in/search masks */
        int *inm = (int *)malloc(n_time * sizeof(int));
        int *srch = (int *)malloc(n_time * sizeof(int));
        if (!inm || !srch) { free(inm); free(srch); continue; }
        int n_in = 0, n_srch = 0;
        for (size_t i = 0; i < n_time; i++) {
            double d = fabs(time[i] - t_centre);
            inm[i] = (d <= temp_half) ? 1 : 0;
            srch[i] = (d <= oot_half) ? 1 : 0;
            n_in += inm[i];
            n_srch += srch[i];
        }
        if (n_in < 3 || n_srch < 6) { free(inm); free(srch); continue; }

        /* cross-correlation over candidate shifts */
        double best = -1.0;
        double best_shift = 0.0;
        for (int si = 0; si < n_shifts; si++) {
            double s = -temp_half + (double)si * shift_step;
            double num = 0.0, den = 0.0;
            for (size_t i = 0; i < n_time; i++) {
                if (!inm[i]) continue;
                double rel = (time[i] - (t_centre + s)) / period;
                double tpl = TEMPLATE_AT(rel);
                num += tpl * (1.0 - flux[i]);
                den += tpl * tpl;
            }
            den += 1e-12;
            if (den > 0.0 && num / den > best) {
                best = num / den;
                best_shift = s;
            }
        }
        if (best < 0.0) { free(inm); free(srch); continue; }
        double t_align = t_centre + best_shift;

        /* matched-filter depth over aligned window */
        int *inm2 = (int *)malloc(n_time * sizeof(int));
        int *ootm = (int *)malloc(n_time * sizeof(int));
        if (!inm2 || !ootm) { free(inm2); free(ootm); free(inm); free(srch); continue; }
        int n_in2 = 0, n_oot = 0;
        for (size_t i = 0; i < n_time; i++) {
            double d = fabs(time[i] - t_align);
            inm2[i] = (d <= temp_half) ? 1 : 0;
            ootm[i] = (d > temp_half && d <= oot_half) ? 1 : 0;
            n_in2 += inm2[i];
            n_oot += ootm[i];
        }
        if (n_in2 < 3 || n_oot < 3) {
            free(inm2); free(ootm); free(inm); free(srch); continue;
        }

        double tpl2_sum = 0.0;
        double *tpl2 = (double *)malloc((size_t)n_in2 * sizeof(double));
        if (!tpl2) { free(inm2); free(ootm); free(inm); free(srch); continue; }
        int j = 0;
        for (size_t i = 0; i < n_time; i++) {
            if (!inm2[i]) continue;
            double rel = (time[i] - t_align) / period;
            tpl2[j] = TEMPLATE_AT(rel);
            tpl2_sum += tpl2[j] * tpl2[j];
            j++;
        }
        double den2 = tpl2_sum + 1e-12;

        double *f_oot = (double *)malloc((size_t)n_oot * sizeof(double));
        if (!f_oot) { free(tpl2); free(inm2); free(ootm); free(inm); free(srch); continue; }
        j = 0;
        for (size_t i = 0; i < n_time; i++)
            if (ootm[i]) f_oot[j++] = flux[i];
        double baseline = median_dbl(f_oot, (size_t)n_oot);
        free(f_oot);
        if (baseline <= 0.0) {
            free(tpl2); free(inm2); free(ootm); free(inm); free(srch); continue;
        }

        double yy_sum = 0.0;
        j = 0;
        for (size_t i = 0; i < n_time; i++) {
            if (!inm2[i]) continue;
            yy_sum += tpl2[j] * (1.0 - flux[i] / baseline);
            j++;
        }
        double depth = yy_sum / den2;
        if (depth < -0.01 || depth > 0.5) {
            free(tpl2); free(inm2); free(ootm); free(inm); free(srch); continue;
        }

        /* uncertainty propagation */
        double *f_oot2 = (double *)malloc((size_t)n_oot * sizeof(double));
        if (!f_oot2) { free(tpl2); free(inm2); free(ootm); free(inm); free(srch); continue; }
        j = 0;
        for (size_t i = 0; i < n_time; i++)
            if (ootm[i]) f_oot2[j++] = flux[i];
        double sigma_oot = std_dbl_ddof1(f_oot2, (size_t)n_oot);
        free(f_oot2);
        double depth_err = (sigma_oot / baseline) / den2 * sqrt(tpl2_sum);
        depth_err = sqrt(depth_err * depth_err +
                         (sigma_oot / baseline) * (sigma_oot / baseline));
        if (!isfinite(depth_err) || depth_err <= 0.0) {
            free(tpl2); free(inm2); free(ootm); free(inm); free(srch); continue;
        }

        if (n_out >= (int)cap) {
            cap *= 2;
            depths = (double *)realloc(depths, cap * sizeof(double));
            errs = (double *)realloc(errs, cap * sizeof(double));
            ns = (int *)realloc(ns, cap * sizeof(int));
            if (!depths || !errs || !ns) {
                free(depths); free(errs); free(ns);
                free(tpl2); free(inm2); free(ootm); free(inm); free(srch);
                return -1;
            }
        }
        depths[n_out] = depth;
        errs[n_out] = depth_err;
        ns[n_out] = n;
        n_out++;
        free(tpl2); free(inm2); free(ootm); free(inm); free(srch);
    }

    out->depths = depths;
    out->depth_errs = errs;
    out->ns = ns;
    out->n = n_out;
    return 0;
    #undef TEMPLATE_AT
}

void zs_extract_free(ZSExtractResult *r) {
    if (!r) return;
    free(r->depths);
    free(r->depth_errs);
    free(r->ns);
    memset(r, 0, sizeof(*r));
}/* part 2: even_odd / depth_consistency / secondary_eclipse / ingress_egress */
#include "zspace_audit.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static int cmp_dbl2(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da < db) ? -1 : (da > db) ? 1 : 0;
}

static double median2(const double *v, size_t n) {
    if (n == 0) return 0.0;
    double *tmp = (double *)malloc(n * sizeof(double));
    if (!tmp) return 0.0;
    memcpy(tmp, v, n * sizeof(double));
    qsort(tmp, n, sizeof(double), cmp_dbl2);
    double med = (n % 2) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
    free(tmp);
    return med;
}

/* ── even_odd (Welch t-test) ───────────────────────────────────────────────── */

void zs_even_odd(const double *depths, const int *ns, int n, ZSEvenOdd *out) {
    memset(out, 0, sizeof(*out));
    if (n < MIN_TRANSITS_FOR_EOTEST) return;

    int n_e = 0, n_o = 0;
    for (int i = 0; i < n; i++) (ns[i] % 2 == 0) ? n_e++ : n_o++;
    if (n_e < 2 || n_o < 2) return;

    double *even = (double *)malloc((size_t)n_e * sizeof(double));
    double *odd = (double *)malloc((size_t)n_o * sizeof(double));
    if (!even || !odd) { free(even); free(odd); return; }
    int je = 0, jo = 0;
    for (int i = 0; i < n; i++) {
        if (ns[i] % 2 == 0) even[je++] = depths[i];
        else odd[jo++] = depths[i];
    }
    double mu_e = 0.0, mu_o = 0.0;
    for (int i = 0; i < n_e; i++) mu_e += even[i];
    for (int i = 0; i < n_o; i++) mu_o += odd[i];
    mu_e /= (double)n_e;
    mu_o /= (double)n_o;
    double sig_e = 0.0, sig_o = 0.0;
    for (int i = 0; i < n_e; i++) { double d = even[i] - mu_e; sig_e += d * d; }
    for (int i = 0; i < n_o; i++) { double d = odd[i] - mu_o; sig_o += d * d; }
    sig_e = sqrt(sig_e / (double)(n_e - 1));
    sig_o = sqrt(sig_o / (double)(n_o - 1));
    free(even); free(odd);

    double se_e = sig_e / sqrt((double)n_e);
    double se_o = sig_o / sqrt((double)n_o);
    double combined = sqrt(se_e * se_e + se_o * se_o);
    double delta_sigma = fabs(mu_e - mu_o) / fmax(combined, 1e-12);

    /* Welch t-test, two-sided p-value (scipy ttest_ind, equal_var=False) */
    double t_stat = (mu_e - mu_o) / fmax(combined, 1e-12);
    double a = sig_e * sig_e / (double)n_e;
    double b = sig_o * sig_o / (double)n_o;
    double df = (a + b) * (a + b) /
                (a * a / (double)(n_e - 1) + b * b / (double)(n_o - 1));
    double p_value = ttest_twosided_p(fabs(t_stat), df);
    if (!(p_value >= 0.0 && p_value <= 1.0)) p_value = 1.0;

    out->n_even = n_e;
    out->n_odd = n_o;
    out->depth_even = mu_e;
    out->depth_odd = mu_o;
    out->depth_even_err = se_e;
    out->depth_odd_err = se_o;
    out->delta_sigma = delta_sigma;
    out->t_stat = fabs(t_stat);
    out->p_value = p_value;
    out->is_eb_flag = (delta_sigma >= EVEN_ODD_SIGMA_THRESHOLD &&
                       p_value < EVEN_ODD_P_VALUE_THRESHOLD) ? 1 : 0;
}

/* ── depth consistency (chi2-reduced CV) ───────────────────────────────────── */

void zs_depth_consistency(const double *depths, const double *depth_errs, int n,
                          ZSDepthCons *out) {
    memset(out, 0, sizeof(*out));
    if (n < 2) {
        out->n = n;
        out->s_depth = 0.50;
        return;
    }
    double mu = 0.0;
    for (int i = 0; i < n; i++) mu += depths[i];
    mu /= (double)n;
    double std = 0.0;
    for (int i = 0; i < n; i++) { double d = depths[i] - mu; std += d * d; }
    std = sqrt(std / (double)(n - 1));
    double cv = std / fmax(mu, 1e-12);

    double sigma_med = median2(depth_errs, (size_t)n);
    double chi2_red;
    if (sigma_med > 0.0) {
        chi2_red = 0.0;
        for (int i = 0; i < n; i++) {
            double d = depths[i] - mu;
            chi2_red += d * d / (depth_errs[i] * depth_errs[i]);
        }
        chi2_red /= (double)n;
    } else {
        chi2_red = 0.0;
        for (int i = 0; i < n; i++) {
            double d = depths[i] - mu;
            chi2_red += d * d;
        }
        chi2_red /= (double)n / fmax(std * std, 1e-12);
    }
    double s_delta = (chi2_red >= 1.0)
        ? fmax(0.0, 1.0 - (chi2_red - 1.0) / (CHI2_NORMALISATION - 1.0))
        : 1.0;
    if (s_delta < 0.0) s_delta = 0.0;
    if (s_delta > 1.0) s_delta = 1.0;

    out->n = n;
    out->mean_depth = mu;
    out->std_depth = std;
    out->cv = cv;
    out->sigma_med = sigma_med;
    out->chi2_red = chi2_red;
    out->s_depth = s_delta;
}

/* ── secondary eclipse ─────────────────────────────────────────────────────── */

void zs_secondary_eclipse(const double *time, size_t n_time, const double *flux,
                          double period, double t0, double duration,
                          ZSSecondary *out) {
    memset(out, 0, sizeof(*out));
    double half_h = (duration / period) / 2.0;

    double sum_pri = 0.0, sum_sec = 0.0, sum_oot = 0.0;
    int n_pri = 0, n_sec = 0, n_oot = 0;
    for (size_t i = 0; i < n_time; i++) {
        double ph = phase_fold(time[i], period, t0);
        double ap = fabs(ph);
        if (ap <= half_h) { sum_pri += flux[i]; n_pri++; }
        else if (fabs(ap - 0.5) <= half_h) { sum_sec += flux[i]; n_sec++; }
        else if (ap > half_h * 2.5 && fabs(ap - 0.5) > half_h * 2.5) {
            sum_oot += flux[i]; n_oot++;
        }
    }
    out->n_primary = n_pri;
    out->n_secondary = n_sec;
    out->n_oot = n_oot;
    if (n_pri < 3 || n_oot < 3) { out->ok = 0; return; }

    double baseline = sum_oot / (double)n_oot;
    double prim = baseline - sum_pri / (double)n_pri;
    double sec = baseline - sum_sec / (double)n_sec;
    prim = fmax(prim, 0.0);
    sec = fmax(sec, 0.0);
    out->primary_depth = prim;
    out->secondary_depth = sec;
    if (prim <= 0.0) { out->ok = 0; return; }

    /* secondary SNR against scatter inside the phase-0.5 window */
    double m = sum_sec / (double)n_sec;
    double s2 = 0.0;
    for (size_t i = 0; i < n_time; i++) {
        double ph = phase_fold(time[i], period, t0);
        double ap = fabs(ph);
        if (fabs(ap - 0.5) <= half_h) {
            double d = flux[i] - m;
            s2 += d * d;
        }
    }
    double sec_std = (n_sec > 1) ? sqrt(s2 / (double)(n_sec - 1)) : 0.0;
    out->secondary_snr = sec / fmax(sec_std / sqrt((double)n_sec), 1e-12);
    out->secondary_ratio = fmin(sec / prim, 3.0);
    out->ok = 1;
}

/* ── ingress/egress trapezoid fit (mirrors scipy curve_fit on 3 params) ────── */

static void trapezoid_model(const double *phase, int n, double depth,
                            double ingress_f, double flat_f, double half_dur,
                            double *out_m) {
    double ingress_w = ingress_f * half_dur;
    double flat_w = flat_f * half_dur;
    for (int i = 0; i < n; i++) {
        double ap = fabs(phase[i]);
        if (ap >= half_dur) {
            out_m[i] = 1.0;
        } else if (ap >= flat_w) {
            double frac = (ap - flat_w) / fmax(ingress_w, 1e-9);
            if (frac < 0.0) frac = 0.0;
            if (frac > 1.0) frac = 1.0;
            out_m[i] = 1.0 - depth * (1.0 - frac);
        } else {
            out_m[i] = 1.0 - depth;
        }
    }
}

void zs_ingress_egress(const double *bin_phase, const double *bin_flux,
                       size_t n_bins, double period, double duration,
                       double transit_depth, ZSIngressEgress *out) {
    memset(out, 0, sizeof(*out));
    out->fp_risk = -1;
    out->fit_ok = 0;

    /* collect valid (finite) bins */
    double *ph = (double *)malloc(n_bins * sizeof(double));
    double *fl = (double *)malloc(n_bins * sizeof(double));
    if (!ph || !fl) { free(ph); free(fl); return; }
    int n = 0;
    for (size_t i = 0; i < n_bins; i++) {
        if (isfinite(bin_flux[i])) { ph[n] = bin_phase[i]; fl[n] = bin_flux[i]; n++; }
    }
    if (n < 10) { free(ph); free(fl); return; }

    double half_dur = (duration / period) / 2.0;

    /* box-constrained least squares via damped Gauss-Newton with bounded
       steps (mirrors curve_fit/trf behaviour: local minimiser from p0) */
    double p[3] = { transit_depth, 0.2, 0.6 };
    const double lo[3] = { 0.0, 0.01, 0.01 };
    const double hi[3] = { 0.5, 0.99, 0.99 };
    double *model = (double *)malloc((size_t)n * sizeof(double));
    double *res = (double *)malloc((size_t)n * sizeof(double));
    double *jac = (double *)malloc((size_t)n * 3 * sizeof(double));
    if (!model || !res || !jac) {
        free(model); free(res); free(jac); free(ph); free(fl);
        return;
    }

    double lam = 1e-3;
    for (int it = 0; it < 200; it++) {
        for (int j = 0; j < 3; j++) {
            if (p[j] < lo[j]) p[j] = lo[j];
            if (p[j] > hi[j]) p[j] = hi[j];
        }
        trapezoid_model(ph, n, p[0], p[1], p[2], half_dur, model);
        double cost = 0.0;
        for (int i = 0; i < n; i++) {
            res[i] = fl[i] - model[i];
            cost += res[i] * res[i];
        }

        /* numeric Jacobian (central differences, clipped) */
        double dp = 1e-6;
        for (int j = 0; j < 3; j++) {
            double p_hi = p[j] + dp, p_lo = p[j] - dp;
            if (p_hi > hi[j]) p_hi = hi[j];
            if (p_lo < lo[j]) p_lo = lo[j];
            double *m_hi = (double *)malloc((size_t)n * sizeof(double));
            double *m_lo = (double *)malloc((size_t)n * sizeof(double));
            if (!m_hi || !m_lo) { free(m_hi); free(m_lo); break; }
            double par[3] = { p[0], p[1], p[2] };
            par[j] = p_hi;
            trapezoid_model(ph, n, par[0], par[1], par[2], half_dur, m_hi);
            par[j] = p_lo;
            trapezoid_model(ph, n, par[0], par[1], par[2], half_dur, m_lo);
            double den = (p_hi - p_lo);
            for (int i = 0; i < n; i++)
                jac[i * 3 + j] = (m_hi[i] - m_lo[i]) / den;
            free(m_hi); free(m_lo);
        }

        /* normal equations with damping */
        double A[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
        double g[3] = { 0, 0, 0 };
        for (int i = 0; i < n; i++) {
            for (int a = 0; a < 3; a++) {
                g[a] += jac[i * 3 + a] * res[i];
                for (int b = 0; b < 3; b++)
                    A[a][b] += jac[i * 3 + a] * jac[i * 3 + b];
            }
        }
        for (int a = 0; a < 3; a++) A[a][a] *= (1.0 + lam);

        /* solve 3x3 via Gauss-Jordan */
        double M[3][4] = {
            { A[0][0], A[0][1], A[0][2], g[0] },
            { A[1][0], A[1][1], A[1][2], g[1] },
            { A[2][0], A[2][1], A[2][2], g[2] },
        };
        for (int c = 0; c < 3; c++) {
            double piv = M[c][c];
            if (fabs(piv) < 1e-14) { piv = (piv < 0) ? -1e-14 : 1e-14; M[c][c] = piv; }
            for (int k = 0; k < 4; k++) M[c][k] /= piv;
            for (int r = 0; r < 3; r++) {
                if (r == c) continue;
                double f = M[r][c];
                for (int k = 0; k < 4; k++) M[r][k] -= f * M[c][k];
            }
        }
        double step[3] = { M[0][3], M[1][3], M[2][3] };
        /* bound each step component (trust-region-like behaviour) */
        for (int j = 0; j < 3; j++) {
            double cap = (j == 0) ? 0.05 : 0.05;
            if (step[j] > cap) step[j] = cap;
            if (step[j] < -cap) step[j] = -cap;
        }
        double step_norm = fabs(step[0]) + fabs(step[1]) + fabs(step[2]);
        if (step_norm < 1e-10) break;

        /* Marquardt accept/reject */
        double trial[3];
        int ok_trial = 1;
        for (int j = 0; j < 3; j++) {
            trial[j] = p[j] + step[j];
            if (trial[j] < lo[j]) trial[j] = lo[j];
            if (trial[j] > hi[j]) trial[j] = hi[j];
        }
        for (int j = 0; j < 3; j++) {
            if (p[j] != trial[j] && trial[j] == p[j]) ok_trial = 0;
        }
        if (ok_trial) {
            double *m_t = (double *)malloc((size_t)n * sizeof(double));
            if (!m_t) { break; }
            trapezoid_model(ph, n, trial[0], trial[1], trial[2], half_dur, m_t);
            double cost_t = 0.0;
            for (int i = 0; i < n; i++) {
                double d = fl[i] - m_t[i];
                cost_t += d * d;
            }
            free(m_t);
            if (cost_t < cost) {
                p[0] = trial[0]; p[1] = trial[1]; p[2] = trial[2];
                lam = fmax(lam / 3.0, 1e-10);
            } else {
                lam = fmin(lam * 10.0, 1e6);
            }
        } else {
            lam = fmin(lam * 10.0, 1e6);
        }
        if (step_norm < 1e-8) break;
    }

    for (int j = 0; j < 3; j++) {
        if (p[j] < lo[j]) p[j] = lo[j];
        if (p[j] > hi[j]) p[j] = hi[j];
    }
    double depth_fit = p[0], ingress_fit = p[1], flat_fit = p[2];

    if (ingress_fit + flat_fit > 1.0) {
        double scale = 0.99 / (ingress_fit + flat_fit);
        ingress_fit *= scale;
        flat_fit *= scale;
    }

    double total_dur_hrs = duration * 24.0;
    double ingress_hrs = ingress_fit * total_dur_hrs;
    double flat_hrs = flat_fit * total_dur_hrs;
    double ingress_fraction = ingress_fit;
    double flat_fraction = flat_fit;
    int is_v_shape = (ingress_fraction > INGRESS_FRACTION_THRESHOLD) ? 1 : 0;
    int fp_risk;
    if (ingress_fraction > 0.45) fp_risk = 2;
    else if (ingress_fraction > 0.35) fp_risk = 1;
    else fp_risk = 0;

    out->depth_fit = depth_fit;
    out->ingress_fraction = ingress_fraction;
    out->flat_fraction = flat_fraction;
    out->ingress_hrs = ingress_hrs;
    out->flat_hrs = flat_hrs;
    out->is_v_shape = is_v_shape;
    out->fp_risk = fp_risk;
    out->fit_ok = 1;

    free(model); free(res); free(jac); free(ph); free(fl);
}