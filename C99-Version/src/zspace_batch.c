/*
 * zspace_batch.c  ·  Full-pipeline batch runner (C99)
 * ====================================================
 * Runs the complete hot path per target:
 *   load/generate LC -> normalize -> flatten(hint) -> BLS ->
 *   re-flatten -> depth extraction -> even/odd -> depth consistency ->
 *   secondary eclipse -> ingress/egress -> ephemeris -> sovereign validate
 *
 * Targets come from a manifest (one per line):
 *   /path/to/lc.csv                    -> real light curve
 *   syn:<P>:<depth>:<snr>:<seed>       -> synthetic target
 *
 * OpenMP parallelises across targets (one thread per target).
 */
#define _POSIX_C_SOURCE 200809L

#include "zspace_bls.h"
#include "zspace_eph.h"
#include "zspace_audit.h"
#include "zspace_core.h"
#include "zspace_ingestion.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ── CSV loading (time,flux[,flux_err,model_flux]) ─────────────────────────── */
static int load_csv(const char *path, double **time, double **flux, size_t *n_out) {
    FILE *fh = fopen(path, "r");
    if (!fh) return -1;
    size_t cap = 65536, n = 0;
    double *t = malloc(cap * sizeof(double));
    double *f = malloc(cap * sizeof(double));
    if (!t || !f) { fclose(fh); free(t); free(f); return -1; }
    char line[1024];
    int header = 1;
    while (fgets(line, sizeof(line), fh)) {
        char *save = NULL;
        char *tok = strtok_r(line, ",", &save);
        if (header) { header = 0; continue; }   /* skip header row */
        if (!tok) continue;
        double tv = atof(tok);
        tok = strtok_r(NULL, ",", &save);
        if (!tok) continue;
        double fv = atof(tok);
        if (n == cap) {
            cap *= 2;
            t = realloc(t, cap * sizeof(double));
            f = realloc(f, cap * sizeof(double));
            if (!t || !f) { fclose(fh); free(t); free(f); return -1; }
        }
        t[n] = tv; f[n] = fv; n++;
    }
    fclose(fh);
    *time = t; *flux = f; *n_out = n;
    return 0;
}

/* ── Median (copies) ───────────────────────────────────────────────────────── */
static int cmp_dbl(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
static double median_dbl(const double *v, size_t n) {
    double *c = malloc(n * sizeof(double));
    if (!c) return 0.0;
    memcpy(c, v, n * sizeof(double));
    qsort(c, n, sizeof(double), cmp_dbl);
    double m = (n % 2) ? c[n / 2] : 0.5 * (c[n / 2 - 1] + c[n / 2]);
    free(c);
    return m;
}

/* ── Mergesort by time (stable) ────────────────────────────────────────────── */
static void mergesort_by_time(double *time, double *flux, size_t n,
                              double *tbuf, double *fbuf) {
    if (n < 2) return;
    size_t mid = n / 2;
    mergesort_by_time(time, flux, mid, tbuf, fbuf);
    mergesort_by_time(time + mid, flux + mid, n - mid, tbuf, fbuf);
    size_t i = 0, j = mid, k = 0;
    while (i < mid && j < n) {
        if (time[j] < time[i]) { tbuf[k] = time[j]; fbuf[k] = flux[j]; j++; }
        else                   { tbuf[k] = time[i]; fbuf[k] = flux[i]; i++; }
        k++;
    }
    while (i < mid) { tbuf[k] = time[i]; fbuf[k] = flux[i]; i++; k++; }
    while (j < n)   { tbuf[k] = time[j]; fbuf[k] = flux[j]; j++; k++; }
    memcpy(time, tbuf, n * sizeof(double));
    memcpy(flux, fbuf, n * sizeof(double));
}

/* ── In-memory flatten (mirrors ingestion._savgol_flatten) ─────────────────── */
static int flatten_inmem(const double *time, const double *flux, size_t n,
                         double period_days, double **flat_out,
                         int *window_out, int *poly_out) {
    if (n < 5) return -2;
    double *ts = malloc(n * sizeof(double));
    double *fs = malloc(n * sizeof(double));
    double *fs_orig = malloc(n * sizeof(double));
    double *tbuf = malloc(n * sizeof(double));
    double *fbuf = malloc(n * sizeof(double));
    double *flat = malloc(n * sizeof(double));
    if (!ts || !fs || !fs_orig || !tbuf || !fbuf || !flat) {
        free(ts); free(fs); free(fs_orig); free(tbuf); free(fbuf); free(flat);
        return -1;
    }
    memcpy(ts, time, n * sizeof(double));
    memcpy(fs, flux, n * sizeof(double));
    mergesort_by_time(ts, fs, n, tbuf, fbuf);

    double cadence = 0.0;
    if (n > 1) {
        double *diffs = malloc((n - 1) * sizeof(double));
        for (size_t i = 1; i < n; i++) diffs[i - 1] = ts[i] - ts[i - 1];
        cadence = median_dbl(diffs, n - 1);
        free(diffs);
    }

    double window_days = (period_days > 0) ? 0.75 * period_days : 3.0;
    int window_pts = (int)lround(window_days / (cadence > 0 ? cadence : 1e-9));
    if (window_pts < 51) window_pts = 51;
    if (window_pts % 2 == 0) window_pts += 1;
    if ((size_t)window_pts > n) window_pts = (n % 2 == 1) ? (int)n : (int)n - 1;
    if (window_pts < 5) window_pts = 5;
    int polyorder = (3 < window_pts - 1) ? 3 : window_pts - 1;

    memcpy(fs_orig, fs, n * sizeof(double));
    ZSLightCurve lc;
    lc.time = ts; lc.flux = fs; lc.flux_err = NULL; lc.quality = NULL; lc.n = n;
    if (zs_savgol_flatten(&lc, window_pts, polyorder) != 0) {
        free(ts); free(fs); free(fs_orig); free(tbuf); free(fbuf); free(flat);
        return -1;
    }
    for (size_t i = 0; i < n; i++) flat[i] = fs_orig[i] / lc.flux[i];

    free(ts); free(fs); free(fs_orig); free(tbuf); free(fbuf);
    *flat_out = flat;
    if (window_out) *window_out = window_pts;
    if (poly_out) *poly_out = polyorder;
    return 0;
}

/* ── One target: full pipeline → JSON string (malloc'd, caller frees) ────────
 * Ladder-aware: mirrors benchmarks_controlled/run_controlled.py evaluate_target:
 *   - single 3.0-day flatten for all audits (faithful to LightCurveIngester)
 *   - BLS with production frequency_factor=20
 *   - local-maxima ladder (k=20, min_rel_snr=0.05) with alias/τ/P filters
 *   - validate ladder in power order until sovereign certifies; else first_status
 */
static char *run_target(const double *time, size_t n_lc, const double *flux,
                        double pmin, double pmax, double flatten_hint,
                        const char *label, int verbosity) {
    double med = median_dbl(flux, n_lc);
    if (med <= 0) return NULL;
    double *norm = malloc(n_lc * sizeof(double));
    if (!norm) return NULL;
    for (size_t i = 0; i < n_lc; i++) norm[i] = flux[i] / med;

    double *flat1 = NULL;
    int w1 = 0, p1 = 0;
    if (flatten_inmem(time, norm, n_lc, flatten_hint, &flat1, &w1, &p1) != 0) {
        free(norm); return NULL;
    }

    ZSBLSConfig cfg = ZS_BLS_CONFIG_DEFAULT;
    cfg.period_min_days = pmin;
    cfg.period_max_days = pmax;
    cfg.frequency_factor = 20.0;
    ZSBLSResult res;
    if (zs_bls_search(time, n_lc, flat1, NULL, &cfg, &res) != 0) {
        free(flat1); free(norm); return NULL;
    }

    /* ── Build ladder (local maxima of power_spectrum, like Python top_candidates) */
    ZSBLSCandidate ladder[64];
    int n_ladder = zs_bls_top_candidates(&res, time, n_lc, flat1, 20, 0.05, ladder, 64);
    if (n_ladder == 0) {
        ladder[0].period_days = res.best_period_days;
        ladder[0].power = res.best_power;
        ladder[0].snr = isfinite(res.best_snr) ? res.best_snr : 0.0;
        ladder[0].fap = isfinite(res.best_fap) ? res.best_fap : 1.0;
        ladder[0].duration_hrs = res.best_duration_hrs;
        ladder[0].t0_days = res.best_t0_days;
        ladder[0].depth = res.best_depth;
        n_ladder = 1;
    } else {
        int has_global = 0;
        for (int i = 0; i < n_ladder; i++) {
            double r = ladder[i].period_days / fmax(res.best_period_days, 1e-9);
            if (fabs(log(r)) < 0.05) { has_global = 1; break; }
        }
        if (!has_global && n_ladder < 64) {
            memmove(&ladder[1], &ladder[0], (size_t)n_ladder * sizeof(ZSBLSCandidate));
            ladder[0].period_days = res.best_period_days;
            ladder[0].power = res.best_power;
            ladder[0].snr = isfinite(res.best_snr) ? res.best_snr : 0.0;
            ladder[0].fap = isfinite(res.best_fap) ? res.best_fap : 1.0;
            ladder[0].duration_hrs = res.best_duration_hrs;
            ladder[0].t0_days = res.best_t0_days;
            ladder[0].depth = res.best_depth;
            n_ladder++;
        }
    }

    /* NO_DETECTION gate: any ladder SNR > 5.5 ? */
    int any_above = 0;
    for (int i = 0; i < n_ladder; i++) if (ladder[i].snr > 5.5) { any_above = 1; break; }

    double *bin_phase = malloc(200 * sizeof(double));
    double *bin_flux = malloc(200 * sizeof(double));
    unsigned char *valid = malloc(200);
    if (!bin_phase || !bin_flux || !valid) {
        zs_bls_result_free(&res); free(bin_phase); free(bin_flux); free(valid);
        free(flat1); free(norm); return NULL;
    }

    /* ── Validate ladder in order until certified ─────────────────────────── */
    char first_status[32] = "FALSE_POSITIVE";
    ZSSovereignCard first_card; memset(&first_card, 0, sizeof(first_card));
    ZSEphResult first_eph; int first_eph_ok = 0;
    ZSEvenOdd first_eo; memset(&first_eo, 0, sizeof(first_eo));
    ZSDepthCons first_dc; memset(&first_dc, 0, sizeof(first_dc));
    ZSSecondary first_se; memset(&first_se, 0, sizeof(first_se));
    ZSIngressEgress first_ie; memset(&first_ie, 0, sizeof(first_ie));
    int has_first = 0;

    ZSSovereignCard certified_card; memset(&certified_card, 0, sizeof(certified_card));
    ZSEphResult certified_eph; int certified_eph_ok = 0;
    ZSEvenOdd certified_eo; ZSDepthCons certified_dc; ZSSecondary certified_se; ZSIngressEgress certified_ie;
    ZSBLSCandidate certified_cand; memset(&certified_cand, 0, sizeof(certified_cand));
    int has_certified = 0;
    int n_tested = 0;

    if (any_above) {
        for (int ci = 0; ci < n_ladder; ci++) {
            ZSBLSCandidate *cc = &ladder[ci];
            double dur_days = cc->duration_hrs / 24.0;
            if (!(dur_days > 1e-9)) continue;
            ZSExtractResult ex;
            if (zs_extract_depths(time, n_lc, flat1, cc->period_days, cc->t0_days, dur_days, &ex) != 0) continue;
            ZSEvenOdd eo; zs_even_odd(ex.depths, ex.ns, ex.n, &eo);
            ZSDepthCons dc; zs_depth_consistency(ex.depths, ex.depth_errs, ex.n, &dc);
            ZSSecondary se; zs_secondary_eclipse(time, n_lc, flat1, cc->period_days, cc->t0_days, dur_days, &se);
            zs_eph_fold_and_bin(time, n_lc, flat1, cc->period_days, cc->t0_days, 200, bin_phase, bin_flux, valid);
            for (int b = 0; b < 200; b++) if (!valid[b]) bin_flux[b] = NAN;
            ZSIngressEgress ie; zs_ingress_egress(bin_phase, bin_flux, 200, cc->period_days, dur_days, fabs(cc->depth), &ie);
            ZSEphResult eph; int eph_ok = zs_eph_resolve(time, n_lc, flat1, cc->period_days, cc->t0_days, pmin, pmax, &eph) == 0;
            ZSCandidate cand; memset(&cand, 0, sizeof(cand));
            cand.period_days = cc->period_days;
            cand.transit_depth = fabs(cc->depth);
            cand.transit_duration_hrs = cc->duration_hrs;
            cand.t0_days = cc->t0_days;
            cand.bls_snr = isfinite(cc->snr) ? cc->snr : 0.0;
            cand.bls_fap = isfinite(cc->fap) ? cc->fap : 1.0;
            cand.stellar_mass_solar = 1.0; cand.stellar_radius_solar = 1.0;
            cand.stellar_teff_k = 5772.0; cand.stellar_logg = 4.44;
            cand.limb_dark_u1 = 0.45; cand.limb_dark_u2 = 0.15;
            cand.even_odd_delta_sigma = eo.delta_sigma;
            cand.shape_ratio = ie.is_v_shape ? 0.3 : 0.9;
            cand.secondary_snr = se.secondary_snr;
            cand.secondary_depth_ratio = se.secondary_ratio;
            cand.alias_secondary_ratio = 0.0;
            cand.coherent_evidence = 0;
            cand.centroid_sigma = 0.0;
            cand.s_periodicity = 0.8; cand.s_depth = 0.8; cand.s_limb = 0.8; cand.s_stellar = 0.8;
            ZSSovereignCard card; zs_sovereign_validate(&cand, time, n_lc, flat1, n_lc, NULL, NULL, &card);
            n_tested++;
            const char *vs = (strcmp(card.sovereign_verdict, "SOVEREIGN_PASS") == 0 ||
                              strcmp(card.sovereign_verdict, "CONDITIONAL_PASS") == 0) ? "OFFLINE_NEW_DISCOVERY" : "FALSE_POSITIVE";
            if (!has_first) {
                first_card = card; first_eph = eph; first_eph_ok = eph_ok;
                first_eo = eo; first_dc = dc; first_se = se; first_ie = ie;
                strncpy(first_status, vs, sizeof(first_status) - 1);
                first_status[sizeof(first_status) - 1] = '\0';
                has_first = 1;
            }
            if (!has_certified && (strcmp(card.sovereign_verdict, "SOVEREIGN_PASS") == 0 ||
                                   strcmp(card.sovereign_verdict, "CONDITIONAL_PASS") == 0)) {
                has_certified = 1;
                certified_cand = *cc;
                certified_card = card; certified_eph = eph; certified_eph_ok = eph_ok;
                certified_eo = eo; certified_dc = dc; certified_se = se; certified_ie = ie;
                zs_extract_free(&ex);
                break;
            }
            zs_extract_free(&ex);
        }
    }

    /* ── Choose output candidate (certified ? certified : global-best path) ── */
    ZSBLSCandidate out_cand;
    ZSSovereignCard out_card;
    ZSEphResult out_eph; int out_eph_ok;
    ZSEvenOdd out_eo; ZSDepthCons out_dc; ZSSecondary out_se; ZSIngressEgress out_ie;
    const char *validation_status;
    int ladder_rank = 0;
    if (!any_above) {
        validation_status = "NO_DETECTION";
        out_cand.period_days = res.best_period_days; out_cand.power = res.best_power;
        out_cand.snr = isfinite(res.best_snr) ? res.best_snr : 0.0;
        out_cand.fap = isfinite(res.best_fap) ? res.best_fap : 1.0;
        out_cand.duration_hrs = res.best_duration_hrs; out_cand.t0_days = res.best_t0_days; out_cand.depth = res.best_depth;
        /* need audits for global for JSON; compute quickly if not already */
        if (has_first) { out_card = first_card; out_eph = first_eph; out_eph_ok = first_eph_ok; out_eo = first_eo; out_dc = first_dc; out_se = first_se; out_ie = first_ie; }
        else {
            double dur_days = out_cand.duration_hrs / 24.0;
            ZSExtractResult ex; int ok = zs_extract_depths(time, n_lc, flat1, out_cand.period_days, out_cand.t0_days, dur_days, &ex) == 0;
            if (ok) { zs_even_odd(ex.depths, ex.ns, ex.n, &out_eo); zs_depth_consistency(ex.depths, ex.depth_errs, ex.n, &out_dc); zs_secondary_eclipse(time,n_lc,flat1,out_cand.period_days,out_cand.t0_days,dur_days,&out_se); zs_eph_fold_and_bin(time,n_lc,flat1,out_cand.period_days,out_cand.t0_days,200,bin_phase,bin_flux,valid); for(int b=0;b<200;b++) if(!valid[b]) bin_flux[b]=NAN; zs_ingress_egress(bin_phase,bin_flux,200,out_cand.period_days,dur_days,fabs(out_cand.depth),&out_ie); out_eph_ok = zs_eph_resolve(time,n_lc,flat1,out_cand.period_days,out_cand.t0_days,pmin,pmax,&out_eph)==0; ZSCandidate cand; memset(&cand,0,sizeof(cand)); cand.period_days=out_cand.period_days; cand.transit_depth=fabs(out_cand.depth); cand.transit_duration_hrs=out_cand.duration_hrs; cand.t0_days=out_cand.t0_days; cand.bls_snr=out_cand.snr; cand.bls_fap=out_cand.fap; cand.stellar_mass_solar=1.0; cand.stellar_radius_solar=1.0; cand.stellar_teff_k=5772.0; cand.stellar_logg=4.44; cand.limb_dark_u1=0.45; cand.limb_dark_u2=0.15; cand.even_odd_delta_sigma=out_eo.delta_sigma; cand.shape_ratio=out_ie.is_v_shape?0.3:0.9; cand.secondary_snr=out_se.secondary_snr; cand.secondary_depth_ratio=out_se.secondary_ratio; cand.alias_secondary_ratio=0.0; cand.coherent_evidence=0; cand.centroid_sigma=0.0; cand.s_periodicity=0.8; cand.s_depth=0.8; cand.s_limb=0.8; cand.s_stellar=0.8; zs_sovereign_validate(&cand,time,n_lc,flat1,n_lc,NULL,NULL,&out_card); zs_extract_free(&ex);} else { memset(&out_card,0,sizeof(out_card)); strcpy(out_card.sovereign_verdict,"FALSE_POSITIVE"); strcpy(out_card.cvs_verdict,"NOT_PLANET"); out_eph_ok=0; }
        }
    } else if (has_certified) {
        validation_status = "OFFLINE_NEW_DISCOVERY";
        out_cand = certified_cand; out_card = certified_card; out_eph = certified_eph; out_eph_ok = certified_eph_ok; out_eo = certified_eo; out_dc = certified_dc; out_se = certified_se; out_ie = certified_ie;
        for (int i = 0; i < n_ladder; i++) if (fabs(log(ladder[i].period_days / out_cand.period_days)) < 0.05) { ladder_rank = i + 1; break; }
    } else {
        validation_status = first_status;
        /* period_found = bls global per Python */
        out_cand.period_days = res.best_period_days; out_cand.power = res.best_power;
        out_cand.snr = isfinite(res.best_snr) ? res.best_snr : 0.0;
        out_cand.fap = isfinite(res.best_fap) ? res.best_fap : 1.0;
        out_cand.duration_hrs = res.best_duration_hrs; out_cand.t0_days = res.best_t0_days; out_cand.depth = res.best_depth;
        out_card = first_card; out_eph = first_eph; out_eph_ok = first_eph_ok; out_eo = first_eo; out_dc = first_dc; out_se = first_se; out_ie = first_ie;
        ladder_rank = 0;
    }

    /* Assemble JSON */
    char buf[8192];
    int blen = snprintf(buf, sizeof(buf),
            "{\"target\":\"%s\","
            "\"period_days\":%.8f,\"power\":%.8f,\"snr\":%.4f,\"fap\":%.8e,"
            "\"duration_hrs\":%.4f,\"t0_days\":%.8f,\"depth\":%.8f,"
            "\"validation_status\":\"%s\","
            "\"ladder_rank\":%d,\"n_candidates_tested\":%d,\"ladder_size\":%d,"
            "\"eph\":{\"ok\":%d,\"classifier\":\"%s\",\"physical_period\":%.8f,"
            "\"multiple\":%d,\"confidence\":%.4f},"
            "\"audit\":{\"even_odd_delta_sigma\":%.6f,\"depth_cv\":%.6f,"
            "\"secondary_snr\":%.4f,\"secondary_ratio\":%.6f,"
            "\"ingress_fraction\":%.6f,\"v_shape\":%d,\"fp_risk\":%d},"
            "\"card\":{\"cvs\":%.6f,\"cvs_verdict\":\"%s\","
            "\"sovereign_verdict\":\"%s\",\"n_transits\":%d,"
            "\"fp_verdict\":\"%s\"},"
            "\"n_points\":%zu}\n",
            label,
            out_cand.period_days, out_cand.power,
            isfinite(out_cand.snr) ? out_cand.snr : 0.0,
            isfinite(out_cand.fap) ? out_cand.fap : 1.0,
            out_cand.duration_hrs, out_cand.t0_days, out_cand.depth,
            validation_status, ladder_rank, n_tested, n_ladder,
            out_eph_ok, out_eph_ok ? out_eph.classifier : "ERROR",
            out_eph_ok ? out_eph.physical_period : 0.0,
            out_eph_ok ? out_eph.multiple : 0,
            out_eph_ok ? out_eph.confidence : 0.0,
            out_eo.delta_sigma, out_dc.cv, out_se.secondary_snr, out_se.secondary_ratio,
            out_ie.ingress_fraction, out_ie.is_v_shape, out_ie.fp_risk,
            out_card.cvs, out_card.cvs_verdict, out_card.sovereign_verdict,
            out_card.n_transits, out_card.fp.overall_verdict,
            n_lc);
    char *out = malloc((size_t)blen + 1);
    if (out) memcpy(out, buf, (size_t)blen + 1);

    free(bin_phase); free(bin_flux); free(valid);
    free(flat1); free(norm);
    zs_bls_result_free(&res);
    (void)verbosity; (void)w1; (void)p1;
    return out;
}

/* ── Batch entry: zspace_card batch manifest.txt out.jsonl [pmin] [pmax] ───── */
int zs_batch_run(const char *manifest_path, const char *out_path,
                 double pmin, double pmax) {
    FILE *mh = fopen(manifest_path, "r");
    if (!mh) { fprintf(stderr, "batch: cannot open manifest\n"); return 2; }
    char lines[4096][512];
    int n_targets = 0;
    while (n_targets < 4096 && fgets(lines[n_targets], sizeof(lines[n_targets]), mh)) {
        size_t len = strlen(lines[n_targets]);
        while (len && (lines[n_targets][len-1] == '\n' || lines[n_targets][len-1] == '\r'))
            lines[n_targets][--len] = 0;
        if (len > 0) n_targets++;
    }
    fclose(mh);
    if (n_targets == 0) { fprintf(stderr, "batch: empty manifest\n"); return 2; }

    double t_start = 0.0;
#ifdef _OPENMP
    omp_set_max_active_levels(1);
    omp_set_nested(0);
    t_start = omp_get_wtime();
#else
    {
        struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
        t_start = (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
    }
#endif

    /* Parallel full-pipeline per target */
    char **results = calloc((size_t)n_targets, sizeof(char *));
    if (!results) return 2;

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < n_targets; i++) {
        const char *spec = lines[i];
        double *time = NULL, *flux = NULL;
        size_t n_lc = 0;
        double flatten_hint = 0.0;
        char label[256];
        int ok = 0;

        if (strncmp(spec, "syn:", 4) == 0) {
            /* syn:<P>:<depth>:<snr>:<seed> */
            double P = atof(spec + 4);
            char *rest = strchr(spec + 4, ':');
            double depth = rest ? atof(rest + 1) : 0.008;
            rest = rest ? strchr(rest + 1, ':') : NULL;
            double snr = rest ? atof(rest + 1) : 8.0;
            rest = rest ? strchr(rest + 1, ':') : NULL;
            long seed = rest ? atol(rest + 1) : (long)(1000 + i);
            ZSLightCurve lc;
            if (zs_generate_synthetic(P, snr, depth, 2.5, 1.0, 60.0, 30.0, &lc) == 0) {
                time = lc.time; flux = lc.flux; n_lc = lc.n;
                flatten_hint = P;
                ok = 1;
            }
            snprintf(label, sizeof(label), "syn_%ld", seed);
            (void)snr;
        } else {
            /* CSV path */
            if (load_csv(spec, &time, &flux, &n_lc) == 0) ok = 1;
            snprintf(label, sizeof(label), "%s", spec);
        }

        if (ok && n_lc >= 50) {
            results[i] = run_target(time, n_lc, flux, pmin, pmax, flatten_hint, label, 0);
        }
        free(time); free(flux);
        if (results[i] == NULL) {
            char tmp[512];
            snprintf(tmp, sizeof(tmp), "{\"target\":\"%s\",\"error\":\"pipeline failed\",\"n_points\":%zu}\n",
                     label, n_lc);
            size_t tlen = strlen(tmp);
            results[i] = malloc(tlen + 1);
            if (results[i]) memcpy(results[i], tmp, tlen + 1);
        }
    }

    double t_end = 0.0;
#ifdef _OPENMP
    t_end = omp_get_wtime();
#else
    {
        struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
        t_end = (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
    }
#endif

    FILE *oh = fopen(out_path, "w");
    if (!oh) { fprintf(stderr, "batch: cannot open output\n"); return 2; }
    for (int i = 0; i < n_targets; i++) {
        if (results[i]) { fputs(results[i], oh); free(results[i]); }
    }
    fclose(oh);
    free(results);

    double total = t_end - t_start;
    double per_target_ms = total * 1000.0 / n_targets;
    fprintf(stdout,
            "{\"batch\":{\"n_targets\":%d,\"wall_sec\":%.4f,\"per_target_ms\":%.4f,"
            "\"targets_per_sec\":%.2f,\"output\":\"%s\"}}\n",
            n_targets, total, per_target_ms, n_targets / (total > 0 ? total : 1e-9),
            out_path);
    return 0;
}