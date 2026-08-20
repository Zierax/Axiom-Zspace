/*
 * ═══════════════════════════════════════════════════════════════════════════
 * zspace_core.h  ·  Axiom-ZSpace C99-Version Sovereign Proof Core
 * ═══════════════════════════════════════════════════════════════════════════
 * Orchestration layer around Purce-generated math kernels
 * (C99-Version/generated/, prov.json for every kernel).
 *
 * Hand-written parts (documented honestly):
 *   - FP-1..FP-10 gate control flow + verdict aggregation
 *   - epoch loop of FP-10 (count_observed_transits)
 *   - input parsing / card emission
 * Reason: Purce v0.1.0 emits per-expression element-wise kernels only
 * (measured: loops and if/else branches are dropped or rejected).
 * All pure math (Kepler, geometry, density, probability, chi2, CVS)
 * is delegated to the Purce-generated kernels.
 *
 * Thresholds mirror zspace_engine/thresholds.py (production profile).
 * ═══════════════════════════════════════════════════════════════════════════
 */
#ifndef ZSPACE_CORE_H
#define ZSPACE_CORE_H

#include <stdint.h>
#include <stddef.h>

/* ── IAU 2015 constants (identical to zspace_engine/constants.py) ─────────── */
#define ZS_G_SI    6.67430e-11
#define ZS_M_SUN   1.9884e30
#define ZS_R_SUN   6.957e8
#define ZS_AU      1.495978707e11
#define ZS_R_EARTH 6.3781e6
#define ZS_PI      3.14159265358979323846

/* ── Thresholds (zspace_engine/thresholds.py, production profile) ─────────── */
#define ZS_FP1_SNR_MIN            5.5
#define ZS_FP2_FAP_MAX            0.05
#define ZS_FP3_EO_SIGMA_MAX       3.0
#define ZS_FP4_SHAPE_MIN          0.4
#define ZS_FP5_SEC_SNR_MAX        3.0
#define ZS_FP5B_SEC_RATIO_MAX     0.30
#define ZS_FP5C_ALIAS_LO          0.20
#define ZS_FP5C_ALIAS_HI          0.90
#define ZS_FP6_CENTROID_MAX       3.0
#define ZS_FP7_DENSITY_LO         0.2
#define ZS_FP7_DENSITY_HI         5.0
#define ZS_FP8_IMPACT_MAX         0.9
#define ZS_FP10_MIN_TRANSITS      2
#define ZS_VERDICT_MAX_FAIL_PASS   2
#define ZS_VERDICT_MAX_FAIL_COND   3
#define ZS_DENSITY_CONFLICT_LO    0.5
#define ZS_DENSITY_CONFLICT_HI    2.0
#define ZS_DENSITY_CONFLICT_SNR   10.0

#define ZS_CVS_W_PERIODICITY 0.97
#define ZS_CVS_W_DEPTH       0.83
#define ZS_CVS_W_LIMB        0.61
#define ZS_CVS_W_SECONDARY   0.31

#define ZS_PARITY_MIN_N      9
#define ZS_PARITY_MIN_SIGMA  3.0
#define ZS_EO_MIN_TRANSITS   4

/* ── Candidate input bundle (mirrors AxiomValidator.validate()) ───────────── */
typedef struct {
    double period_days;
    double transit_depth;
    double transit_duration_hrs;
    double t0_days;
    double stellar_mass_solar;
    double stellar_radius_solar;
    double stellar_teff_k;
    double stellar_logg;
    double planet_radius_earth;
    double bls_snr;
    double bls_fap;
    double even_odd_delta_sigma;
    double shape_ratio;
    double secondary_snr;
    double secondary_depth_ratio;
    double alias_secondary_ratio;
    int    coherent_evidence;
    double centroid_sigma;
    double limb_dark_u1;
    double limb_dark_u2;
    double cvs_score;
    double s_periodicity;
    double s_depth;
    double s_limb;
    double s_stellar;
} ZSCandidate;

/* ── Section results (numeric core, mirrors ProofEngine dicts) ────────────── */
typedef struct {
    double a_m;
    double a_au;
    double residual_si_pct;
    double residual_solar_pct;
    int    verdict_pass;   /* PASS=1 / WARN=0 */
} ZSKeplerResult;

typedef struct {
    double k;
    double k_sq;
    double delta_ld_corrected;
    double consistency_residual_pct;
    double rp_earth;
    int    verdict_pass;
} ZSGeometryResult;

typedef struct {
    double a_over_rs_transit;
    double a_over_rs_direct;
    double rho_transit_gcc;
    double rho_tic_gcc;
    double density_ratio;
    double logg_calc;
    double logg_residual;
    int    is_tdur_placeholder;
    int    is_eb_density_flag;
    int    verdict;          /* 0 PASS, 1 WARN, 2 FAIL */
} ZSDensityResult;

typedef struct {
    double P_tr;
    double impact_parameter_b;
    double ingress_hrs;
    double i_min_deg;
    int    is_grazing;
    int    verdict_pass;
} ZSProbabilityResult;

typedef struct {
    int    n_tests;
    int    n_pass;
    int    n_fail;
    int    n_critical;
    int    n_critical_pass;
    int    fp_verdicts[12];   /* 0 FAIL, 1 PASS, -1 skipped */
    char   overall_verdict[32];
    int    conflict_snr_density;
    int    conflict_snr_shape;
} ZSFalsePositiveResult;

typedef struct {
    double chi2;
    double reduced_chi2;
    int    dof;
} ZSChiSquaredResult;

/* ── Sovereign card output ────────────────────────────────────────────────── */
typedef struct {
    ZSKeplerResult      kepler;
    ZSGeometryResult    geometry;
    ZSDensityResult     density;
    ZSProbabilityResult probability;
    ZSFalsePositiveResult fp;
    ZSChiSquaredResult  chi_squared;
    double cvs;
    char   cvs_verdict[32];
    char   sovereign_verdict[32];
    int    n_transits;
    int    all_sections_pass;
} ZSSovereignCard;

/* ── Core entry points ────────────────────────────────────────────────────── */
void zs_compute_kepler(const ZSCandidate *c, ZSKeplerResult *out);
void zs_compute_geometry(const ZSCandidate *c, ZSGeometryResult *out);
void zs_compute_density(const ZSCandidate *c, const ZSKeplerResult *kep, ZSDensityResult *out);
void zs_compute_probability(const ZSCandidate *c, const ZSKeplerResult *kep, ZSProbabilityResult *out);
void zs_compute_cvs(const ZSCandidate *c, double *cvs, char *verdict, size_t verdict_sz);

/* FP-10: count independent parity-supporting transit series (0..2).
 * time/flux are the (finite) light-curve arrays; returns -1 when the
 * light curve is absent/degenerate (test skipped). */
int zs_count_observed_transits(const double *time, size_t n_time,
                               const double *flux, size_t n_flux,
                               double period_days, double t0_days,
                               double duration_days);

void zs_false_positive_ruling(const ZSCandidate *c, const ZSDensityResult *den,
                              const ZSProbabilityResult *prob,
                              int n_transits, ZSFalsePositiveResult *out);

void zs_sovereign_validate(const ZSCandidate *c,
                           const double *time, size_t n_time,
                           const double *flux, size_t n_flux,
                           const double *flux_err, const double *model_flux,
                           ZSSovereignCard *out);

#endif /* ZSPACE_CORE_H */