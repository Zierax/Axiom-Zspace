/*
 * ═══════════════════════════════════════════════════════════════════════════
 * zspace_card.c  ·  C99-Version Sovereign Card CLI
 * ═══════════════════════════════════════════════════════════════════════════
 * Reads a candidate bundle (key=value lines) plus an optional light-curve
 * CSV (time,flux[,flux_err,model_flux]), runs the sovereign proof engine
 * and prints the Sovereign Logic Card as JSON lines.

 * Usage:
 *   ./zspace_card [candidate.txt] [lc.csv]
 * Candidate file keys (ZSCandidate):
 *   period_days transit_depth transit_duration_hrs t0_days
 *   stellar_mass_solar stellar_radius_solar stellar_teff_k stellar_logg
 *   planet_radius_earth bls_snr bls_fap even_odd_delta_sigma shape_ratio
 *   secondary_snr secondary_depth_ratio alias_secondary_ratio coherent_evidence
 *   centroid_sigma limb_dark_u1 limb_dark_u2 cvs_score
 *   s_periodicity s_depth s_limb s_stellar
 * ═══════════════════════════════════════════════════════════════════════════
 */
#include "zspace_core.h"
#include "zspace_bls.h"
#include "zspace_eph.h"
#include "zspace_audit.h"
#include "zspace_ingestion.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

static void die(const char *msg) {
    fprintf(stderr, "zspace_card: %s\n", msg);
    exit(2);
}

static double parse_dbl(const char *s) {
    char *end = NULL;
    double v = strtod(s, &end);
    if (end == s) die("invalid numeric value");
    return v;
}

static int cmp_dbl_asc(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static int parse_int(const char *s) {
    char *end = NULL;
    long v = strtol(s, &end, 10);
    if (end == s) die("invalid integer value");
    return (int)v;
}

static void parse_candidate(const char *path, ZSCandidate *c) {
    FILE *fh = fopen(path, "r");
    if (!fh) die("cannot open candidate file");
    char line[512];
    while (fgets(line, sizeof(line), fh)) {
        char *eq = strchr(line, '=');
        if (!eq) continue;
        *eq = '\0';
        char *key = line;
        char *val = eq + 1;
        /* trim */
        while (*key == ' ' || *key == '\t') key++;
        char *ke = key + strlen(key);
        while (ke > key && (ke[-1] == ' ' || ke[-1] == '\t' || ke[-1] == '\r' || ke[-1] == '\n')) *--ke = '\0';
        char *ve = val + strlen(val);
        while (ve > val && (ve[-1] == ' ' || ve[-1] == '\t' || ve[-1] == '\r' || ve[-1] == '\n')) *--ve = '\0';

        if      (!strcmp(key, "period_days"))            c->period_days = parse_dbl(val);
        else if (!strcmp(key, "transit_depth"))          c->transit_depth = parse_dbl(val);
        else if (!strcmp(key, "transit_duration_hrs"))   c->transit_duration_hrs = parse_dbl(val);
        else if (!strcmp(key, "t0_days"))                c->t0_days = parse_dbl(val);
        else if (!strcmp(key, "stellar_mass_solar"))     c->stellar_mass_solar = parse_dbl(val);
        else if (!strcmp(key, "stellar_radius_solar"))   c->stellar_radius_solar = parse_dbl(val);
        else if (!strcmp(key, "stellar_teff_k"))         c->stellar_teff_k = parse_dbl(val);
        else if (!strcmp(key, "stellar_logg"))           c->stellar_logg = parse_dbl(val);
        else if (!strcmp(key, "planet_radius_earth"))    c->planet_radius_earth = parse_dbl(val);
        else if (!strcmp(key, "bls_snr"))                c->bls_snr = parse_dbl(val);
        else if (!strcmp(key, "bls_fap"))                c->bls_fap = parse_dbl(val);
        else if (!strcmp(key, "even_odd_delta_sigma"))   c->even_odd_delta_sigma = parse_dbl(val);
        else if (!strcmp(key, "shape_ratio"))            c->shape_ratio = parse_dbl(val);
        else if (!strcmp(key, "secondary_snr"))          c->secondary_snr = parse_dbl(val);
        else if (!strcmp(key, "secondary_depth_ratio"))  c->secondary_depth_ratio = parse_dbl(val);
        else if (!strcmp(key, "alias_secondary_ratio"))  c->alias_secondary_ratio = parse_dbl(val);
        else if (!strcmp(key, "coherent_evidence"))      c->coherent_evidence = parse_int(val);
        else if (!strcmp(key, "centroid_sigma"))         c->centroid_sigma = parse_dbl(val);
        else if (!strcmp(key, "limb_dark_u1"))           c->limb_dark_u1 = parse_dbl(val);
        else if (!strcmp(key, "limb_dark_u2"))           c->limb_dark_u2 = parse_dbl(val);
        else if (!strcmp(key, "cvs_score"))              c->cvs_score = parse_dbl(val);
        else if (!strcmp(key, "s_periodicity"))          c->s_periodicity = parse_dbl(val);
        else if (!strcmp(key, "s_depth"))                c->s_depth = parse_dbl(val);
        else if (!strcmp(key, "s_limb"))                 c->s_limb = parse_dbl(val);
        else if (!strcmp(key, "s_stellar"))              c->s_stellar = parse_dbl(val);
        /* unknown keys ignored */
    }
    fclose(fh);
}

/* CSV: time,flux[,flux_err,model_flux] — one header line skipped.
   Returns n rows and sets *has_full = 1 when 4 columns are present. */
static size_t load_lc(const char *path,
                      double **time, double **flux,
                      double **flux_err, double **model_flux,
                      int *has_full) {
    FILE *fh = fopen(path, "r");
    if (!fh) die("cannot open light curve file");
    size_t cap = 4096, n = 0;
    *has_full = 0;
    *time = (double *)malloc(cap * sizeof(double));
    *flux = (double *)malloc(cap * sizeof(double));
    *flux_err = (double *)malloc(cap * sizeof(double));
    *model_flux = (double *)malloc(cap * sizeof(double));
    char line[512];
    int first = 1;
    while (fgets(line, sizeof(line), fh)) {
        if (first) { first = 0; continue; }   /* skip header */
        double t, f, e = 0.0, m = 0.0;
        int ncol = sscanf(line, "%lf,%lf,%lf,%lf", &t, &f, &e, &m);
        if (ncol < 2) continue;
        const char *p = line;
        int commas = 0;
        while (*p) { if (*p == ',') commas++; p++; }
        if (commas >= 3 && ncol == 4) *has_full = 1;
        if (n >= cap) {
            cap *= 2;
            *time = (double *)realloc(*time, cap * sizeof(double));
            *flux = (double *)realloc(*flux, cap * sizeof(double));
            *flux_err = (double *)realloc(*flux_err, cap * sizeof(double));
            *model_flux = (double *)realloc(*model_flux, cap * sizeof(double));
        }
        (*time)[n] = t;
        (*flux)[n] = f;
        (*flux_err)[n] = e;
        (*model_flux)[n] = m;
        n++;
    }
    fclose(fh);
    return n;
}

static void print_card(const ZSCandidate *c, const ZSSovereignCard *card,
                       size_t n_lc, const char *lc_path) {
    printf("{\n");
    printf("  \"schema\": \"Axiom-ZSpace Sovereign Logic Card v1.0 (C99-Version)\",\n");
    printf("  \"sovereign_verdict\": \"%s\",\n", card->sovereign_verdict);
    printf("  \"cvs\": %.8f,\n", card->cvs);
    printf("  \"cvs_verdict\": \"%s\",\n", card->cvs_verdict);
    printf("  \"all_sections_pass\": %s,\n", card->all_sections_pass ? "true" : "false");
    printf("  \"n_transits\": %d,\n", card->n_transits);
    printf("  \"section_1_kepler\": {\n");
    printf("    \"a_au\": %.7f,\n", card->kepler.a_au);
    printf("    \"residual_si_pct\": %.9f,\n", card->kepler.residual_si_pct);
    printf("    \"residual_solar_pct\": %.6f,\n", card->kepler.residual_solar_pct);
    printf("    \"verdict\": \"%s\"\n", card->kepler.verdict_pass ? "PASS" : "WARN");
    printf("  },\n");
    printf("  \"section_2_geometry\": {\n");
    printf("    \"k\": %.8f,\n", card->geometry.k);
    printf("    \"delta_ld_corrected\": %.8f,\n", card->geometry.delta_ld_corrected);
    printf("    \"consistency_residual_pct\": %.4f,\n", card->geometry.consistency_residual_pct);
    printf("    \"rp_earth\": %.4f,\n", card->geometry.rp_earth);
    printf("    \"verdict\": \"%s\"\n", card->geometry.verdict_pass ? "PASS" : "WARN");
    printf("  },\n");
    printf("  \"section_3_density\": {\n");
    printf("    \"rho_transit_gcc\": %.5f,\n", card->density.rho_transit_gcc);
    printf("    \"rho_tic_gcc\": %.5f,\n", card->density.rho_tic_gcc);
    printf("    \"density_ratio\": %.5f,\n", card->density.density_ratio);
    printf("    \"logg_calc\": %.4f,\n", card->density.logg_calc);
    printf("    \"logg_residual\": %.4f,\n", card->density.logg_residual);
    printf("    \"is_tdur_placeholder\": %d,\n", card->density.is_tdur_placeholder);
    printf("    \"is_eb_density_flag\": %d,\n", card->density.is_eb_density_flag);
    printf("    \"verdict\": \"%s\"\n",
           card->density.verdict == 0 ? "PASS" : card->density.verdict == 1 ? "WARN" : "FAIL");
    printf("  },\n");
    printf("  \"section_4_probability\": {\n");
    printf("    \"P_tr\": %.6f,\n", card->probability.P_tr);
    printf("    \"impact_parameter_b\": %.5f,\n", card->probability.impact_parameter_b);
    printf("    \"ingress_hrs\": %.4f,\n", card->probability.ingress_hrs);
    printf("    \"i_min_deg\": %.4f,\n", card->probability.i_min_deg);
    printf("    \"is_grazing\": %d,\n", card->probability.is_grazing);
    printf("    \"verdict\": \"%s\"\n", card->probability.verdict_pass ? "PASS" : "WARN");
    printf("  },\n");
    printf("  \"section_5_fp_ruling\": {\n");
    printf("    \"n_tests\": %d,\n", card->fp.n_tests);
    printf("    \"n_pass\": %d,\n", card->fp.n_pass);
    printf("    \"n_fail\": %d,\n", card->fp.n_fail);
    printf("    \"n_critical\": %d,\n", card->fp.n_critical);
    printf("    \"n_critical_pass\": %d,\n", card->fp.n_critical_pass);
    printf("    \"conflict_snr_density\": %d,\n", card->fp.conflict_snr_density);
    printf("    \"conflict_snr_shape\": %d,\n", card->fp.conflict_snr_shape);
    printf("    \"fp_verdicts\": [");
    for (int i = 0; i < 12; i++) {
        printf("%s%d", i ? "," : "", card->fp.fp_verdicts[i]);
    }
    printf("],\n");
    printf("    \"overall_verdict\": \"%s\"\n", card->fp.overall_verdict);
    printf("  },\n");
    if (n_lc > 0 && card->chi_squared.dof > 0) {
        printf("  \"section_6_chi_squared\": {\n");
        printf("    \"chi2\": %.6f,\n", card->chi_squared.chi2);
        printf("    \"reduced_chi2\": %.6f,\n", card->chi_squared.reduced_chi2);
        printf("    \"dof\": %d\n", card->chi_squared.dof);
        printf("  },\n");
    }
    printf("  \"inputs\": {\n");
    printf("    \"period_days\": %.7f,\n", c->period_days);
    printf("    \"transit_depth\": %.8f,\n", c->transit_depth);
    printf("    \"transit_duration_hrs\": %.4f,\n", c->transit_duration_hrs);
    printf("    \"bls_snr\": %.4f,\n", c->bls_snr);
    printf("    \"bls_fap\": %.6f,\n", c->bls_fap);
    printf("    \"lc_points\": %zu,\n", n_lc);
    printf("    \"lc_file\": \"%s\"\n", lc_path ? lc_path : "-");
    printf("  }\n");
    printf("}\n");
}

static int run_bls(int argc, char **argv);
static int run_eph(int argc, char **argv);
static int run_audit(int argc, char **argv);
static int run_flatten(int argc, char **argv);
static int run_pipeline_all(int argc, char **argv);

int zs_batch_run(const char *manifest_path, const char *out_path,
                 double pmin, double pmax);

static int run_batch(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: zspace_card batch manifest.txt out.jsonl [pmin] [pmax]\n");
        return 2;
    }
    double pmin = (argc > 4) ? atof(argv[4]) : 0.5;
    double pmax = (argc > 5) ? atof(argv[5]) : 13.5;
    return zs_batch_run(argv[2], argv[3], pmin, pmax);
}

int main(int argc, char **argv) {
    if (argc > 1 && strcmp(argv[1], "bls") == 0)
        return run_bls(argc, argv);
    if (argc > 1 && strcmp(argv[1], "eph") == 0)
        return run_eph(argc, argv);
    if (argc > 1 && strcmp(argv[1], "audit") == 0)
        return run_audit(argc, argv);
    if (argc > 1 && strcmp(argv[1], "flatten") == 0)
        return run_flatten(argc, argv);
    if (argc > 1 && strcmp(argv[1], "pipeline") == 0)
        return run_pipeline_all(argc, argv);
    if (argc > 1 && strcmp(argv[1], "batch") == 0)
        return run_batch(argc, argv);

    ZSCandidate c;
    memset(&c, 0, sizeof(c));
    /* sane defaults (typical F/G dwarf, central transit) */
    c.stellar_mass_solar = 1.0;
    c.stellar_radius_solar = 1.0;
    c.stellar_teff_k = 5772.0;
    c.stellar_logg = 4.44;
    c.limb_dark_u1 = 0.45;
    c.limb_dark_u2 = 0.15;
    c.s_periodicity = 0.8;
    c.s_depth = 0.8;
    c.s_limb = 0.8;
    c.s_stellar = 0.8;

    const char *cand_path = (argc > 1) ? argv[1] : NULL;
    const char *lc_path = (argc > 2 && argv[2][0] != '\0') ? argv[2] : NULL;

    if (cand_path) parse_candidate(cand_path, &c);

    double *time = NULL, *flux = NULL, *flux_err = NULL, *model_flux = NULL;
    size_t n_lc = 0;
    int has_full = 0;
    if (lc_path) {
        n_lc = load_lc(lc_path, &time, &flux, &flux_err, &model_flux, &has_full);
        if (n_lc == 0) die("empty light curve");
    }
    if (!has_full) { flux_err = NULL; model_flux = NULL; }

    ZSSovereignCard card;
    zs_sovereign_validate(&c, time, n_lc, flux, n_lc, flux_err, model_flux, &card);

    print_card(&c, &card, n_lc, lc_path);

    free(time); free(flux); free(flux_err); free(model_flux);
    return 0;
}

/* ── BLS subcommand: zspace_card bls lc.csv [period_min] [period_max] ───────── */
static int run_bls(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: zspace_card bls lc.csv period_min period_max\n");
        return 2;
    }
    const char *lc_path = argv[2];
    const double pmin = atof(argv[3]);
    const double pmax = atof(argv[4]);

    double *time = NULL, *flux = NULL, *flux_err = NULL, *model_flux = NULL;
    int has_full = 0;
    size_t n_lc = load_lc(lc_path, &time, &flux, &flux_err, &model_flux, &has_full);
    if (n_lc == 0) die("empty light curve");

    ZSBLSConfig cfg = ZS_BLS_CONFIG_DEFAULT;
    cfg.period_min_days = pmin;
    cfg.period_max_days = pmax;

    ZSBLSResult res;
    if (zs_bls_search(time, n_lc, flux, NULL, &cfg, &res) != 0)
        die("BLS search failed");

    printf("{\n");
    printf("  \"schema\": \"Axiom-ZSpace BLS Periodogram v1.0 (C99-Version)\",\n");
    printf("  \"period_days\": %.8f,\n", res.best_period_days);
    printf("  \"power\": %.8f,\n", res.best_power);
    printf("  \"snr\": %.4f,\n", isfinite(res.best_snr) ? res.best_snr : 0.0);
    printf("  \"fap\": %.8e,\n", isfinite(res.best_fap) ? res.best_fap : 1.0);
    printf("  \"duration_hrs\": %.4f,\n", res.best_duration_hrs);
    printf("  \"t0_days\": %.8f,\n", res.best_t0_days);
    printf("  \"depth\": %.8f,\n", res.best_depth);
    printf("  \"lc_points\": %zu\n", n_lc);
    printf("}\n");

    free(time); free(flux); free(flux_err); free(model_flux);
    zs_bls_result_free(&res);
    return 0;
}

/* ── EPH subcommand: zspace_card eph lc.csv period t0 period_min period_max ─── */
static int run_eph(int argc, char **argv) {
    if (argc < 7) {
        fprintf(stderr, "usage: zspace_card eph lc.csv period t0 period_min period_max\n");
        return 2;
    }
    const char *lc_path = argv[2];
    const double period = atof(argv[3]);
    const double t0 = atof(argv[4]);
    const double pmin = atof(argv[5]);
    const double pmax = atof(argv[6]);

    double *time = NULL, *flux = NULL, *flux_err = NULL, *model_flux = NULL;
    int has_full = 0;
    size_t n_lc = load_lc(lc_path, &time, &flux, &flux_err, &model_flux, &has_full);
    if (n_lc == 0) die("empty light curve");

    ZSEphResult res;
    if (zs_eph_resolve(time, n_lc, flux, period, t0, pmin, pmax, &res) != 0)
        die("ephemeris resolution failed");

    printf("{\n");
    printf("  \"schema\": \"Axiom-ZSpace Ephemeris Resolution v1.0 (C99-Version)\",\n");
    printf("  \"candidate_period\": %.8f,\n", res.candidate_period);
    printf("  \"physical_period\": %.8f,\n", res.physical_period);
    printf("  \"multiple\": %d,\n", res.multiple);
    printf("  \"classifier\": \"%s\",\n", res.classifier);
    printf("  \"confidence\": %.4f,\n", res.confidence);
    printf("  \"evidence\": \"%s\",\n", res.evidence);
    printf("  \"pattern\": \"%s\",\n", res.pattern);
    printf("  \"flags\": [");
    for (int i = 0; i < res.n_flags; i++)
        printf("%s\"%s\"", i ? ", " : "", res.flags[i]);
    printf("],\n");
    printf("  \"n_transits\": %d,\n", (int)((pmax - pmin) / period + 0.5));
    printf("  \"lc_points\": %zu\n", n_lc);
    printf("}\n");

    free(time); free(flux); free(flux_err); free(model_flux);
    return 0;
}

/* ── AUDIT subcommand: zspace_card audit lc.csv period t0 dur_hrs depth ────── */
static int run_audit(int argc, char **argv) {
    if (argc < 7) {
        fprintf(stderr,
                "usage: zspace_card audit lc.csv period t0 duration_hrs transit_depth\n");
        return 2;
    }
    const char *lc_path = argv[2];
    const double period = atof(argv[3]);
    const double t0 = atof(argv[4]);
    const double duration = atof(argv[5]) / 24.0;
    const double transit_depth = atof(argv[6]);

    double *time = NULL, *flux = NULL, *flux_err = NULL, *model_flux = NULL;
    int has_full = 0;
    size_t n_lc = load_lc(lc_path, &time, &flux, &flux_err, &model_flux, &has_full);
    if (n_lc == 0) die("empty light curve");

    ZSExtractResult ex;
    if (zs_extract_depths(time, n_lc, flux, period, t0, duration, &ex) != 0)
        die("depth extraction failed");

    ZSEvenOdd eo;
    zs_even_odd(ex.depths, ex.ns, ex.n, &eo);

    ZSDepthCons dc;
    zs_depth_consistency(ex.depths, ex.depth_errs, ex.n, &dc);

    ZSSecondary se;
    zs_secondary_eclipse(time, n_lc, flux, period, t0, duration, &se);

    double *bin_phase = (double *)malloc(200 * sizeof(double));
    double *bin_flux = (double *)malloc(200 * sizeof(double));
    unsigned char *valid = (unsigned char *)malloc(200);
    if (!bin_phase || !bin_flux || !valid) die("alloc failed");
    zs_eph_fold_and_bin(time, n_lc, flux, period, t0, 200, bin_phase, bin_flux, valid);
    for (int b = 0; b < 200; b++)
        if (!valid[b]) bin_flux[b] = NAN;

    ZSIngressEgress ie;
    zs_ingress_egress(bin_phase, bin_flux, 200, period, duration, transit_depth, &ie);

    static const char *fp_risk_names[] = { "LOW", "MEDIUM", "HIGH", "UNKNOWN" };

    printf("{\n");
    printf("  \"schema\": \"Axiom-ZSpace Transit Audit v1.0 (C99-Version)\",\n");
    printf("  \"even_odd\": {\n");
    printf("    \"n_even\": %d,\n", eo.n_even);
    printf("    \"n_odd\": %d,\n", eo.n_odd);
    printf("    \"depth_even\": %.8f,\n", eo.depth_even);
    printf("    \"depth_odd\": %.8f,\n", eo.depth_odd);
    printf("    \"depth_even_err\": %.8f,\n", eo.depth_even_err);
    printf("    \"depth_odd_err\": %.8f,\n", eo.depth_odd_err);
    printf("    \"delta_sigma\": %.6f,\n", eo.delta_sigma);
    printf("    \"t_stat\": %.6f,\n", eo.t_stat);
    printf("    \"p_value\": %.6e,\n", eo.p_value);
    printf("    \"is_eb_flag\": %s\n", eo.is_eb_flag ? "true" : "false");
    printf("  },\n");
    printf("  \"depth_consistency\": {\n");
    printf("    \"n_transits\": %d,\n", dc.n);
    printf("    \"mean_depth\": %.8f,\n", dc.mean_depth);
    printf("    \"std_depth\": %.8f,\n", dc.std_depth);
    printf("    \"cv\": %.6f,\n", dc.cv);
    printf("    \"sigma_med\": %.8f,\n", dc.sigma_med);
    printf("    \"chi2_red\": %.6f,\n", dc.chi2_red);
    printf("    \"s_depth\": %.6f\n", dc.s_depth);
    printf("  },\n");
    printf("  \"secondary_eclipse\": {\n");
    printf("    \"primary_depth\": %.8f,\n", se.primary_depth);
    printf("    \"secondary_depth\": %.8f,\n", se.secondary_depth);
    printf("    \"secondary_ratio\": %.6f,\n", se.secondary_ratio);
    printf("    \"secondary_snr\": %.4f,\n", se.secondary_snr);
    printf("    \"n_primary\": %d,\n", se.n_primary);
    printf("    \"n_secondary\": %d,\n", se.n_secondary);
    printf("    \"ok\": %s\n", se.ok ? "true" : "false");
    printf("  },\n");
    printf("  \"ingress_egress\": {\n");
    printf("    \"depth_fit\": %.8f,\n", ie.depth_fit);
    printf("    \"ingress_fraction\": %.6f,\n", ie.ingress_fraction);
    printf("    \"flat_fraction\": %.6f,\n", ie.flat_fraction);
    printf("    \"ingress_hrs\": %.4f,\n", ie.ingress_hrs);
    printf("    \"flat_hrs\": %.4f,\n", ie.flat_hrs);
    printf("    \"is_v_shape\": %s,\n", ie.is_v_shape ? "true" : "false");
    printf("    \"fp_risk\": \"%s\",\n",
           ie.fp_risk < 0 ? fp_risk_names[3] : fp_risk_names[ie.fp_risk]);
    printf("    \"fit_ok\": %s\n", ie.fit_ok ? "true" : "false");
    printf("  },\n");
    printf("  \"lc_points\": %zu\n", n_lc);
    printf("}\n");

    free(bin_phase); free(bin_flux); free(valid);
    free(time); free(flux); free(flux_err); free(model_flux);
    zs_extract_free(&ex);
    return 0;
}

/* ── flatten: Savitzky-Golay flatten matching ingestion._savgol_flatten ───── */

static void mergesort_by_time(double *time, double *flux, size_t n, double *tbuf, double *fbuf) {
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

static double median_dbl(const double *v, size_t n) {
    double *c = malloc(n * sizeof(double));
    if (!c) return 0.0;
    memcpy(c, v, n * sizeof(double));
    qsort(c, n, sizeof(double), cmp_dbl_asc);
    double m = (n % 2) ? c[n / 2] : 0.5 * (c[n / 2 - 1] + c[n / 2]);
    free(c);
    return m;
}

/* In-memory flatten: sorted-SG (window from period hint) → flux/trend, restored
 * to original time order.  Matches ingestion._savgol_flatten. */
static int flatten_inmem(const double *time, const double *flux, size_t n,
                         double period_days, double **flat_out,
                         int *window_out, int *poly_out, double *cadence_out) {
    if (n < 5) return -2;
    size_t *order = malloc(n * sizeof(size_t));
    double *ts = malloc(n * sizeof(double));
    double *fs = malloc(n * sizeof(double));
    double *fs_orig = malloc(n * sizeof(double));
    double *tbuf = malloc(n * sizeof(double));
    double *fbuf = malloc(n * sizeof(double));
    double *flat = malloc(n * sizeof(double));
    if (!order || !ts || !fs || !fs_orig || !tbuf || !fbuf || !flat) {
        free(order); free(ts); free(fs); free(fs_orig); free(tbuf); free(fbuf); free(flat);
        return -1;
    }
    for (size_t i = 0; i < n; i++) { order[i] = i; ts[i] = time[i]; fs[i] = flux[i]; }
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
    lc.time = ts;
    lc.flux = fs;
    lc.flux_err = NULL;
    lc.quality = NULL;
    lc.n = n;
    if (zs_savgol_flatten(&lc, window_pts, polyorder) != 0) {
        free(order); free(ts); free(fs); free(fs_orig); free(tbuf); free(fbuf); free(flat);
        return -1;
    }
    for (size_t i = 0; i < n; i++) flat[order[i]] = fs_orig[i] / lc.flux[i];

    free(order); free(ts); free(fs); free(fs_orig); free(tbuf); free(fbuf);
    *flat_out = flat;
    *window_out = window_pts;
    *poly_out = polyorder;
    if (cadence_out) *cadence_out = cadence;
    return 0;
}

static int run_flatten(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: zspace_card flatten in.csv out.csv period_days\n");
        return 2;
    }
    const char *in_path = argv[2];
    const char *out_path = argv[3];
    const double period_days = atof(argv[4]);

    double *time = NULL, *flux = NULL, *flux_err = NULL, *model_flux = NULL;
    int has_full = 0;
    size_t n = load_lc(in_path, &time, &flux, &flux_err, &model_flux, &has_full);
    if (n == 0) die("empty light curve");
    if (n < 5) {
        FILE *of = fopen(out_path, "w");
        if (!of) die("cannot open output");
        fprintf(of, "time,flux\n");
        for (size_t i = 0; i < n; i++) fprintf(of, "%.8f,%.10f\n", time[i], flux[i]);
        fclose(of);
        free(time); free(flux); free(flux_err); free(model_flux);
        printf("{\"points\": %zu, \"window_pts\": 0, \"polyorder\": 0}\n", n);
        return 0;
    }

    size_t *order = malloc(n * sizeof(size_t));
    double *ts = malloc(n * sizeof(double));
    double *fs = malloc(n * sizeof(double));
    double *tbuf = malloc(n * sizeof(double));
    double *fbuf = malloc(n * sizeof(double));
    if (!order || !ts || !fs || !tbuf || !fbuf) die("alloc failed");

    for (size_t i = 0; i < n; i++) { order[i] = i; ts[i] = time[i]; fs[i] = flux[i]; }
    mergesort_by_time(ts, fs, n, tbuf, fbuf);

    double cadence = 0.0;
    if (n > 1) {
        double *diffs = malloc((n - 1) * sizeof(double));
        for (size_t i = 1; i < n; i++) diffs[i - 1] = ts[i] - ts[i - 1];
        cadence = median_dbl(diffs, n - 1);
        free(diffs);
    }

    double *flat_orig = NULL;
    int window_pts = 0, polyorder = 0;
    int rc = flatten_inmem(time, flux, n, period_days, &flat_orig, &window_pts,
                           &polyorder, &cadence);
    if (rc != 0) die("SG flatten failed");

    FILE *of = fopen(out_path, "w");
    if (!of) die("cannot open output");
    fprintf(of, "time,flux\n");
    for (size_t i = 0; i < n; i++) fprintf(of, "%.8f,%.10f\n", time[i], flat_orig[i]);
    fclose(of);

    printf("{\"points\": %zu, \"window_pts\": %d, \"polyorder\": %d, \"cadence_days\": %.10f}\n",
           n, window_pts, polyorder, cadence);

    free(flat_orig);
    free(time); free(flux); free(flux_err); free(model_flux);
    return 0;
}

/* ── pipeline: flatten + BLS + audits in a single process ──────────────────── */

static void print_audit_json(const ZSExtractResult *ex, const ZSEvenOdd *eo,
                             const ZSDepthCons *dc, const ZSSecondary *se,
                             const ZSIngressEgress *ie) {
    static const char *fp_risk_names[] = { "LOW", "MEDIUM", "HIGH", "UNKNOWN" };
    printf("  \"even_odd\": {\n");
    printf("    \"n_even\": %d,\n", eo->n_even);
    printf("    \"n_odd\": %d,\n", eo->n_odd);
    printf("    \"depth_even\": %.8f,\n", eo->depth_even);
    printf("    \"depth_odd\": %.8f,\n", eo->depth_odd);
    printf("    \"depth_even_err\": %.8f,\n", eo->depth_even_err);
    printf("    \"depth_odd_err\": %.8f,\n", eo->depth_odd_err);
    printf("    \"delta_sigma\": %.6f,\n", eo->delta_sigma);
    printf("    \"t_stat\": %.6f,\n", eo->t_stat);
    printf("    \"p_value\": %.6e,\n", eo->p_value);
    printf("    \"is_eb_flag\": %s\n", eo->is_eb_flag ? "true" : "false");
    printf("  },\n");
    printf("  \"depth_consistency\": {\n");
    printf("    \"n_transits\": %d,\n", dc->n);
    printf("    \"mean_depth\": %.8f,\n", dc->mean_depth);
    printf("    \"std_depth\": %.8f,\n", dc->std_depth);
    printf("    \"cv\": %.6f,\n", dc->cv);
    printf("    \"sigma_med\": %.8f,\n", dc->sigma_med);
    printf("    \"chi2_red\": %.6f,\n", dc->chi2_red);
    printf("    \"s_depth\": %.6f\n", dc->s_depth);
    printf("  },\n");
    printf("  \"secondary_eclipse\": {\n");
    printf("    \"primary_depth\": %.8f,\n", se->primary_depth);
    printf("    \"secondary_depth\": %.8f,\n", se->secondary_depth);
    printf("    \"secondary_ratio\": %.6f,\n", se->secondary_ratio);
    printf("    \"secondary_snr\": %.4f,\n", se->secondary_snr);
    printf("    \"n_primary\": %d,\n", se->n_primary);
    printf("    \"n_secondary\": %d,\n", se->n_secondary);
    printf("    \"ok\": %s\n", se->ok ? "true" : "false");
    printf("  },\n");
    printf("  \"ingress_egress\": {\n");
    printf("    \"depth_fit\": %.8f,\n", ie->depth_fit);
    printf("    \"ingress_fraction\": %.6f,\n", ie->ingress_fraction);
    printf("    \"flat_fraction\": %.6f,\n", ie->flat_fraction);
    printf("    \"ingress_hrs\": %.4f,\n", ie->ingress_hrs);
    printf("    \"flat_hrs\": %.4f,\n", ie->flat_hrs);
    printf("    \"is_v_shape\": %s,\n", ie->is_v_shape ? "true" : "false");
    printf("    \"fp_risk\": \"%s\",\n",
           ie->fp_risk < 0 ? fp_risk_names[3] : fp_risk_names[ie->fp_risk]);
    printf("    \"fit_ok\": %s\n", ie->fit_ok ? "true" : "false");
    printf("  },\n");
    (void)ex;
}

static int run_pipeline_all(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr,
                "usage: zspace_card pipeline lc.csv period_min period_max flatten_hint_period\n");
        return 2;
    }
    const char *lc_path = argv[2];
    const double pmin = atof(argv[3]);
    const double pmax = atof(argv[4]);
    const double flatten_hint = atof(argv[5]);

    double *time = NULL, *flux = NULL, *flux_err = NULL, *model_flux = NULL;
    int has_full = 0;
    size_t n_lc = load_lc(lc_path, &time, &flux, &flux_err, &model_flux, &has_full);
    if (n_lc == 0) die("empty light curve");

    double med = median_dbl(flux, n_lc);
    if (med <= 0) die("non-positive median flux");
    double *norm = malloc(n_lc * sizeof(double));
    if (!norm) die("alloc failed");
    for (size_t i = 0; i < n_lc; i++) norm[i] = flux[i] / med;

    double *flat1 = NULL, *flat2 = NULL;
    int w1 = 0, w2 = 0, p1 = 0, p2 = 0;
    double cad1 = 0.0, cad2 = 0.0;
    if (flatten_inmem(time, norm, n_lc, flatten_hint, &flat1, &w1, &p1, &cad1) != 0)
        die("flatten failed");

    ZSBLSConfig cfg = ZS_BLS_CONFIG_DEFAULT;
    cfg.period_min_days = pmin;
    cfg.period_max_days = pmax;
    ZSBLSResult res;
    if (zs_bls_search(time, n_lc, flat1, NULL, &cfg, &res) != 0)
        die("BLS search failed");

    if (flatten_inmem(time, norm, n_lc, res.best_period_days, &flat2, &w2, &p2, &cad2) != 0)
        die("re-flatten failed");

    double dur_days = res.best_duration_hrs / 24.0;
    ZSExtractResult ex;
    if (zs_extract_depths(time, n_lc, flat2, res.best_period_days, res.best_t0_days,
                          dur_days, &ex) != 0)
        die("depth extraction failed");
    ZSEvenOdd eo;
    zs_even_odd(ex.depths, ex.ns, ex.n, &eo);
    ZSDepthCons dc;
    zs_depth_consistency(ex.depths, ex.depth_errs, ex.n, &dc);
    ZSSecondary se;
    zs_secondary_eclipse(time, n_lc, flat2, res.best_period_days, res.best_t0_days,
                         dur_days, &se);

    double *bin_phase = (double *)malloc(200 * sizeof(double));
    double *bin_flux = (double *)malloc(200 * sizeof(double));
    unsigned char *valid = (unsigned char *)malloc(200);
    if (!bin_phase || !bin_flux || !valid) die("alloc failed");
    zs_eph_fold_and_bin(time, n_lc, flat2, res.best_period_days, res.best_t0_days,
                        200, bin_phase, bin_flux, valid);
    for (int b = 0; b < 200; b++)
        if (!valid[b]) bin_flux[b] = NAN;
    ZSIngressEgress ie;
    zs_ingress_egress(bin_phase, bin_flux, 200, res.best_period_days, dur_days,
                      fabs(res.best_depth), &ie);

    printf("{\n");
    printf("  \"schema\": \"Axiom-ZSpace Pipeline v1.0 (C99-Version)\",\n");
    printf("  \"period_days\": %.8f,\n", res.best_period_days);
    printf("  \"power\": %.8f,\n", res.best_power);
    printf("  \"snr\": %.4f,\n", isfinite(res.best_snr) ? res.best_snr : 0.0);
    printf("  \"fap\": %.8e,\n", isfinite(res.best_fap) ? res.best_fap : 1.0);
    printf("  \"duration_hrs\": %.4f,\n", res.best_duration_hrs);
    printf("  \"t0_days\": %.8f,\n", res.best_t0_days);
    printf("  \"depth\": %.8f,\n", res.best_depth);
    printf("  \"window_pts\": %d,\n", w2);
    printf("  \"polyorder\": %d,\n", p2);
    print_audit_json(&ex, &eo, &dc, &se, &ie);
    printf("  \"flattened\": [");
    for (size_t i = 0; i < n_lc; i++)
        printf("%s%.10f", i ? "," : "", flat2[i]);
    printf("],\n");
    printf("  \"lc_points\": %zu\n", n_lc);
    printf("}\n");

    free(bin_phase); free(bin_flux); free(valid);
    free(flat1); free(flat2); free(norm);
    free(time); free(flux); free(flux_err); free(model_flux);
    zs_extract_free(&ex);
    zs_bls_result_free(&res);
    return 0;
}