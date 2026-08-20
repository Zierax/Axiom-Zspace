/*
 * ═══════════════════════════════════════════════════════════════════════════
 * zspace_core.c  ·  Sovereign Proof Engine — C99 orchestration layer
 * ═══════════════════════════════════════════════════════════════════════════
 * Calls the Purce-generated math kernels (generated/purce_src_*.c, every
 * one differentially verified against numpy, see tests/) and implements
 * the §1–§6 gate logic of zspace_engine/validator.py.
 *
 * Hand-written here (documented limitation of Purce v0.1.0):
 *   - control flow: gate evaluation, verdict aggregation, epoch loop
 *   - arccos (i_min_deg, informational only): Purce emits #error
 *   - median/std of the light curve for FP-10 (array algorithms)
 *
 * All thresholds mirror zspace_engine/thresholds.py production profile.
 * ═══════════════════════════════════════════════════════════════════════════
 */
#include "zspace_core.h"

#include "purce_src.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── constants (identical to zspace_engine/constants.py) ─────────────────── */
#define ZS_M_SUN_S  1.9884e30
#define ZS_R_SUN_S  6.957e8
#define ZS_R_EARTH_S 6.3781e6

/* single-element helper buffers (kernels are element-wise, n = 1) */
#define DECL1(v)   double v[1]
#define SET1(v, x) ((v)[0] = (x))
#define OUT1(o)    (o)[0]

/* kernel runners (Purce signatures: input(s) + output, constants inlined) */
static double k1(void (*fn)(int, const double *restrict, double *restrict),
                 double a) {
    DECL1(A); DECL1(R);
    SET1(A, a); R[0] = 0.0;
    fn(1, A, R);
    return OUT1(R);
}

static double k2(void (*fn)(int, const double *restrict, const double *restrict, double *restrict),
                 double a, double b) {
    DECL1(A); DECL1(B); DECL1(R);
    SET1(A, a); SET1(B, b); R[0] = 0.0;
    fn(1, A, B, R);
    return OUT1(R);
}

/* ───────────────────────────────────────────────────────────────────────────
 * §1 Keplerian Dynamics
 * ─────────────────────────────────────────────────────────────────────────── */
void zs_compute_kepler(const ZSCandidate *c, ZSKeplerResult *out) {
    double M_clamped = c->stellar_mass_solar;
    if (M_clamped <= 0.0 || M_clamped > 100.0)
        M_clamped = fmin(100.0, fmax(0.08, M_clamped));

    const double P_sec = c->period_days * 86400.0;
    const double inv3  = 1.0 / 3.0;
    const double four  = 4.0;
    const double four_pi_sq = 4.0 * ZS_PI * ZS_PI;
    const double eps_a = 1e-30, eps_e = 1e-10;
    const double hundred = 100.0;
    const double inv365 = 1.0 / 365.25;

    double gm_partial = k2(purce_src_kepler_gm_cce91793, ZS_G_SI, M_clamped);
    double gm          = k2(purce_src_kepler_gm_full_aa78c33d, gm_partial, ZS_M_SUN_S);
    double p_sq        = k1(purce_src_kepler_p_sec_sq_d786d462, P_sec);
    double a3_num      = k2(purce_src_kepler_a3_num_d9ffca74, gm, p_sq);
    double pi_sq       = k1(purce_src_kepler_pi_sq_0e1b6818, ZS_PI);
    double denom       = k2(purce_src_kepler_denominator_ff600e16, four, pi_sq);
    double a3          = k2(purce_src_kepler_a3_08bcd3f7, a3_num, denom);
    double a_m         = k2(purce_src_kepler_a_m_f674c7f2, a3, inv3);
    double a_au        = k2(purce_src_kepler_a_au_e9f15bbc, a_m, ZS_AU);

    /* SI residual */
    double p_sq2   = k1(purce_src_kepler_p_sec_sq_2_ed24ce94, P_sec);
    double a_cubed = k1(purce_src_kepler_a_m_cubed_9f77db4b, a_m);
    double ratio_si = k2(purce_src_kepler_ratio_si_30d9bfbf, p_sq2, a_cubed);
    double exp_si   = k2(purce_src_kepler_expected_si_983db113, four_pi_sq, gm);
    double diff_si  = k2(purce_src_kepler_diff_si_ebc5792c, ratio_si, exp_si);
    double abs_si   = k1(purce_src_kepler_diff_abs_48aeecda, diff_si);
    double norm_si  = k2(purce_src_kepler_norm_diff_si_4ed6e8b4, abs_si, exp_si);
    double res_si   = k2(purce_src_kepler_residual_si_3df9ed0f, norm_si, hundred);

    /* solar residual */
    double p_yr     = k2(purce_src_kepler_p_yr_8a8b968e, c->period_days, inv365);
    double p_yr_sq  = k1(purce_src_kepler_p_yr_sq_dad6e103, p_yr);
    double au_cubed = k1(purce_src_kepler_a_au_cubed_f9db1b7f, a_au);
    double au_max   = k2(purce_src_kepler_a_au_max_2517d228, au_cubed, eps_a);
    double ratio    = k2(purce_src_kepler_ratio_solar_19fb5df0, p_yr_sq, au_max);
    double expect   = k2(purce_src_kepler_expect_solar_bfe8eecc, 1.0, M_clamped);
    double diff     = k2(purce_src_kepler_diff_solar_01a0f6f9, ratio, expect);
    double abs_d    = k1(purce_src_kepler_diff_solar_abs_7f8924f3, diff);
    double exp_max  = k2(purce_src_kepler_expect_max_34d0f1bf, expect, eps_e);
    double norm_d   = k2(purce_src_kepler_norm_solar_167c610a, abs_d, exp_max);
    double res_sol  = k2(purce_src_kepler_residual_solar_cfad1534, norm_d, hundred);

    out->a_m = a_m;
    out->a_au = a_au;
    out->residual_si_pct = res_si;
    out->residual_solar_pct = res_sol;
    out->verdict_pass = (res_si < 0.001) ? 1 : 0;
}

/* ───────────────────────────────────────────────────────────────────────────
 * §2 Geometric Consistency
 * ─────────────────────────────────────────────────────────────────────────── */
void zs_compute_geometry(const ZSCandidate *c, ZSGeometryResult *out) {
    const double eps   = 1e-12;
    const double one   = 1.0;
    const double inv3  = 1.0 / 3.0;
    const double inv6  = 1.0 / 6.0;
    const double floor_ld = 0.1;
    const double hundred = 100.0;
    const double exp1  = 1.0, exp2 = 2.0;

    double d_clamped = k2(purce_src_geometry_delta_clamped_f1c04251, c->transit_depth, eps);
    double k         = k1(purce_src_geometry_k_a714c37c, d_clamped);
    double k_sq      = k1(purce_src_geometry_k_sq_a9e9f55b, k);

    double u1t = k2(purce_src_geometry_u1_term_0e734a16, c->limb_dark_u1, inv3);
    double u2t = k2(purce_src_geometry_u2_term_a239ecee, c->limb_dark_u2, inv6);
    double s1  = k2(purce_src_geometry_ld_sub_5dd1c288, one, u1t);
    double s2  = k2(purce_src_geometry_ld_sub2_f9f8a586, s1, u2t);
    double i_mean = k2(purce_src_geometry_ld_i_mean_73906914, s2, floor_ld);

        /* (1 - mu_c)^e with mu_c = 1 (central transit, b=0) */
    double zpd1 = k2(purce_src_geometry_ld_zeropow_diff_e2b6d51d, one, 1.0);
    double zpd2 = k2(purce_src_geometry_ld_zeropow_diff_e2b6d51d, one, 1.0);
    double zp  = k2(purce_src_geometry_ld_zeropow_3e6ae1b6, zpd1, exp1);
    double zps = k2(purce_src_geometry_ld_zeropow_3e6ae1b6, zpd2, exp2);
    double u1c = k2(purce_src_geometry_ld_u1cen_term_97c1eb73, c->limb_dark_u1, zp);
    double u2c = k2(purce_src_geometry_ld_u2cen_term_86c94295, c->limb_dark_u2, zps);
    double c1  = k2(purce_src_geometry_ld_cen_sub_3e670b09, one, u1c);
    double c2  = k2(purce_src_geometry_ld_cen_sub2_30a34eea, c1, u2c);
    double i_cen = k2(purce_src_geometry_ld_i_cen_1d5b45bc, c2, floor_ld);

    double num = k2(purce_src_geometry_ld_numerator_052155c8, k_sq, i_cen);
    double delta_ld = k2(purce_src_geometry_ld_corrected_d478fad9, num, i_mean);

    /* consistency residual |delta - k²| / max(k², eps) * 100 */
    double cs = k2(purce_src_geometry_cons_resid_sub_2bc0676c, c->transit_depth, k_sq);
    double ca = k1(purce_src_geometry_cons_resid_abs_58db32c0, cs);
    double cm = k2(purce_src_geometry_cons_k_sq_max_72c33749, k_sq, eps);
    double cd = k2(purce_src_geometry_cons_resid_div_7e29062d, ca, cm);
    double cp = k2(purce_src_geometry_cons_resid_pct_cdf32450, cd, hundred);

    double rp_m = k2(purce_src_geometry_rp_m_0bc355fa, k, c->stellar_radius_solar * ZS_R_SUN_S);
    double rp_e = k2(purce_src_geometry_rp_earth_d541fdd4, rp_m, ZS_R_EARTH_S);

    out->k = k;
    out->k_sq = k_sq;
    out->delta_ld_corrected = delta_ld;
    out->consistency_residual_pct = cp;
    out->rp_earth = rp_e;
    out->verdict_pass = (cp < 5.0) ? 1 : 0;
}

/* ───────────────────────────────────────────────────────────────────────────
 * §3 Stellar Density Constraint
 * ─────────────────────────────────────────────────────────────────────────── */
void zs_compute_density(const ZSCandidate *c, const ZSKeplerResult *kep, ZSDensityResult *out) {
    const double P_sec = c->period_days * 86400.0;
    const double T_sec = c->transit_duration_hrs * 3600.0;
    const double one   = 1.0;
    const double three = 3.0;
    const double three_pi = 3.0 * ZS_PI;
    const double four_thirds = 4.0 / 3.0;
    const double inv1000 = 1.0 / 1000.0;
    const double hundred = 100.0;
    const double eps_r = 1e-10;

    double opk = k2(purce_src_density_one_plus_k_aef82505,
                    sqrt(fmax(c->transit_depth, 1e-10)), one);

    double tmax = k2(purce_src_density_tdur_max_730f8bd8, T_sec, one);
    double anum = k2(purce_src_density_a_over_rs_num_dbe138f3, P_sec, opk);
    double aden = k2(purce_src_density_a_over_rs_den_d67f7aa5, ZS_PI, tmax);
    double a_or_rs_transit = k2(purce_src_density_a_over_rs_transit_6ece4a52, anum, aden);

    double a_or_rs_direct = kep->a_m / fmax(c->stellar_radius_solar * ZS_R_SUN_S, 1e-9);

    double p2    = k1(purce_src_density_rho_pow_p_4b908d4e, P_sec);
    double rden  = k2(purce_src_density_rho_den_e36e01d8, ZS_G_SI, p2);
    double coeff = k2(purce_src_density_rho_coeff_4b89eabc, three_pi, rden);
    double rpow  = k2(purce_src_density_rho_pow_a_170d3d0c, a_or_rs_transit, three);
    double rho_kgm3 = k2(purce_src_density_rho_transit_kgm3_24e14fe1, coeff, rpow);
    double rho_transit_gcc = k2(purce_src_density_rho_transit_gcc_654aae3e, rho_kgm3, inv1000);

    double M_star_kg = c->stellar_mass_solar * ZS_M_SUN_S;
    double R_star_m  = c->stellar_radius_solar * ZS_R_SUN_S;
    double vc = k2(purce_src_density_vol_coeff_4509092e, four_thirds, ZS_PI);
    double vp = k1(purce_src_density_vol_pow_3c8c6d79, R_star_m);
    double vol = k2(purce_src_density_vol_5c746dc2, vc, vp);
    double rho_tic_kgm3 = k2(purce_src_density_rho_tic_kgm3_ab7e90f2, M_star_kg, vol);
    double rho_tic_gcc = k2(purce_src_density_rho_tic_gcc_190d4c2c, rho_tic_kgm3, inv1000);

    double rtmax = k2(purce_src_density_rho_tic_max_a44dbdda, rho_tic_gcc, eps_r);
    double dratio = k2(purce_src_density_ratio_67fcbc7c, rho_transit_gcc, rtmax);

    double gnum = k2(purce_src_density_g_num_93b6e99f, ZS_G_SI, M_star_kg);
    double gden = k1(purce_src_density_g_den_a8d7e46a, R_star_m);
    double gdiv = k2(purce_src_density_g_cgs_div_2a470693, gnum, gden);
    double g_cgs = k2(purce_src_density_g_cgs_659a9574, gdiv, hundred);
    double gmax = k2(purce_src_density_g_cgs_max_12048b9e, g_cgs, one);
    double logg_calc = k1(purce_src_density_logg_3eed90ec, gmax);

    int is_tdur_placeholder = (fabs(c->transit_duration_hrs - 1.0) < 0.001) ? 1 : 0;
    int is_eb_density_flag = 0;
    if (!is_tdur_placeholder)
        is_eb_density_flag = (dratio > 2.0 || dratio < 0.5) ? 1 : 0;

    double logg_residual = fabs(logg_calc - c->stellar_logg);

    out->a_over_rs_transit = a_or_rs_transit;
    out->a_over_rs_direct = a_or_rs_direct;
    out->rho_transit_gcc = rho_transit_gcc;
    out->rho_tic_gcc = rho_tic_gcc;
    out->density_ratio = dratio;
    out->logg_calc = logg_calc;
    out->logg_residual = logg_residual;
    out->is_tdur_placeholder = is_tdur_placeholder;
    out->is_eb_density_flag = is_eb_density_flag;

    /* PASS / WARN / FAIL (validator.py:524-528) */
    out->verdict = 0;
    if (is_eb_density_flag && !is_tdur_placeholder)
        out->verdict = 2;
    else if (is_tdur_placeholder)
        out->verdict = 1;
    else
        out->verdict = (!is_eb_density_flag && logg_residual < 0.3) ? 0 : 1;
}

/* ───────────────────────────────────────────────────────────────────────────
 * §4 Probability of Transit
 * ─────────────────────────────────────────────────────────────────────────── */
void zs_compute_probability(const ZSCandidate *c, const ZSKeplerResult *kep, ZSProbabilityResult *out) {
    const double P_sec = c->period_days * 86400.0;
    const double T_sec = c->transit_duration_hrs * 3600.0;
    const double one = 1.0;
    const double eps = 1e-9;
    const double R_star_m = c->stellar_radius_solar * ZS_R_SUN_S;
    const double R_p_m = c->planet_radius_earth * ZS_R_EARTH_S;
    const double k = sqrt(fmax(c->transit_depth, 1e-10));

    double num = k2(purce_src_probability_num_ab1d7dec, R_star_m, R_p_m);
    double den = k2(purce_src_probability_den_dc355a2e, kep->a_m, num);
    double ptr = k2(purce_src_probability_ptr_ca056473, num, den);
    double P_tr = k2(purce_src_probability_ptr_clamped_360acecd, ptr, one);

    double rmax = k2(purce_src_probability_rstar_max_772c4022, R_star_m, eps);
    double a_or_rs = k2(purce_src_probability_a_over_rs_b15534e1, kep->a_m, rmax);
    double pmax = k2(purce_src_probability_p_max_bdde5587, P_sec, one);
    double sratio = k2(purce_src_probability_sin_ratio_bc35f349, T_sec, pmax);
    double sprod = k2(purce_src_probability_sin_prod_08cbc925, ZS_PI, sratio);
    double sarg = k2(purce_src_probability_sin_arg_e4b6b0e4, sprod, ZS_PI);
    double sval = k1(purce_src_probability_sin_9b57fd41, sarg);
    double sterm = k2(purce_src_probability_sin_term_0841a329, a_or_rs, sval);

    double opk = k2(purce_src_probability_one_plus_k_92c81f69, k, one);
    double opk2 = k1(purce_src_probability_opk_sq_e716b82f, opk);
    double st2 = k1(purce_src_probability_st_sq_2ff19abf, sterm);
    double bsq = k2(purce_src_probability_b_sq_feea2f20, opk2, st2);
    double bcl = k1(purce_src_probability_b_clamped_5af44d93, bsq);
    double b = k1(purce_src_probability_b_a73daf1d, bcl);

    double iden = k2(purce_src_probability_ingress_den_ca99298e, opk, eps);
    double irat = k2(purce_src_probability_ingress_ratio_2703541e, k, iden);
    double ingress = k2(purce_src_probability_ingress_hrs_01a860d8, c->transit_duration_hrs, irat);

    double cmax = k2(purce_src_probability_cos_max_9fb76630, num, kep->a_m);
    /* arccos not supported by Purce v0.1.0 — informational only (see header) */
    double i_min_deg = acos(fmin(cmax, 1.0)) * 180.0 / ZS_PI;

    out->P_tr = P_tr;
    out->impact_parameter_b = b;
    out->ingress_hrs = ingress;
    out->i_min_deg = i_min_deg;
    out->is_grazing = (b > 0.9) ? 1 : 0;
    out->verdict_pass = out->is_grazing ? 0 : 1;
}

/* ───────────────────────────────────────────────────────────────────────────
 * CVS Composite Vitality Score
 * ─────────────────────────────────────────────────────────────────────────── */
void zs_compute_cvs(const ZSCandidate *c, double *cvs, char *verdict, size_t verdict_sz) {
    const double w_p = ZS_CVS_W_PERIODICITY;
    const double w_d = ZS_CVS_W_DEPTH;
    const double w_l = ZS_CVS_W_LIMB;
    const double w_s = ZS_CVS_W_SECONDARY;
    const double eps = 1e-12;

    double wp_sp = k2(purce_src_cvs_wp_sp_cc7ba620, w_p, c->s_periodicity);
    double wd_sd = k2(purce_src_cvs_wd_sd_ffed4e1e, w_d, c->s_depth);
    double wl_sl = k2(purce_src_cvs_wl_sl_65d565ca, w_l, c->s_limb);
    double ws_ss = k2(purce_src_cvs_ws_ss_0fcc6119, w_s, c->s_stellar);

    double na = k2(purce_src_cvs_num_a_70589ee8, wp_sp, wd_sd);
    double nb = k2(purce_src_cvs_num_b_b310441a, wl_sl, ws_ss);
    double num = k2(purce_src_cvs_numerator_9533d661, na, nb);

    double da = k2(purce_src_cvs_den_a_b65a4875, w_p, w_d);
    double db = k2(purce_src_cvs_den_b_c209fe64, w_l, w_s);
    double tot = k2(purce_src_cvs_total_119628ce, da, db);
    double tmax = k2(purce_src_cvs_total_max_49a252a0, tot, eps);
    double score = k2(purce_src_cvs_score_5bac9461, num, tmax);

    *cvs = score;
    snprintf(verdict, verdict_sz, "%s", score > 0.80 ? "PLANET" : "NOT_PLANET");
}

/* ───────────────────────────────────────────────────────────────────────────
 * FP-10: count independent parity-supporting transit series (0..2)
 * Mirrors validator.py:1468 count_observed_transits exactly.
 * ─────────────────────────────────────────────────────────────────────────── */
static int cmp_dbl(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static double median_of(const double *x, size_t n) {
    if (n == 0) return 0.0;
    double *copy = (double *)malloc(n * sizeof(double));
    if (!copy) return 0.0;
    memcpy(copy, x, n * sizeof(double));
    qsort(copy, n, sizeof(double), cmp_dbl);
    double med = (n % 2) ? copy[n / 2] : 0.5 * (copy[n / 2 - 1] + copy[n / 2]);
    free(copy);
    return med;
}

int zs_count_observed_transits(const double *time, size_t n_time,
                               const double *flux, size_t n_flux,
                               double period_days, double t0_days,
                               double duration_days) {
    if (!time || !flux || n_time < 20 || n_flux < 20) return -1;
    if (n_time != n_flux) return -1;

    /* keep finite pairs */
    double *t = (double *)malloc(n_time * sizeof(double));
    double *f = (double *)malloc(n_flux * sizeof(double));
    if (!t || !f) { free(t); free(f); return -1; }
    size_t m = 0;
    for (size_t i = 0; i < n_time; i++) {
        if (isfinite(time[i]) && isfinite(flux[i])) {
            t[m] = time[i]; f[m] = flux[i]; m++;
        }
    }
    if (m < 20) { free(t); free(f); return -1; }

    double P = period_days;
    double t0 = t0_days;
    double D = fmax(duration_days, 1e-3);
    double half = 0.5 * D;

    double baseline = median_of(f, m);
    double sq = 0.0;
    for (size_t i = 0; i < m; i++) {
        double d = f[i] - baseline;
        sq += d * d;
    }
    double sigma = sqrt(sq / (double)m);   /* np.std(f - median), population */

    if (!(sigma > 0.0) || !isfinite(sigma)) { free(t); free(f); return 0; }

    double tmin = t[0], tmax = t[0];
    for (size_t i = 1; i < m; i++) {
        if (t[i] < tmin) tmin = t[i];
        if (t[i] > tmax) tmax = t[i];
    }

    long long k_start = (long long)floor((tmin - t0) / P) - 1;
    long long k_end   = (long long)ceil((tmax - t0) / P) + 1;

    double pooled_n[2] = {0.0, 0.0};
    double pooled_sum[2] = {0.0, 0.0};

    for (long long k = k_start; k <= k_end; k++) {
        double e = t0 + (double)k * P;
        double lo = e - half, hi = e + half;
        size_t n_cad = 0;
        double sum_f = 0.0;
        for (size_t i = 0; i < m; i++) {
            if (t[i] >= lo && t[i] <= hi) { n_cad++; sum_f += f[i]; }
        }
        if (n_cad < 3) continue;
        int par = (int)(k % 2);
        if (par < 0) par += 2;
        pooled_n[par] += (double)n_cad;
        pooled_sum[par] += sum_f;
    }

    int n_supporting = 0;
    for (int par = 0; par < 2; par++) {
        if (pooled_n[par] < 9.0) continue;
        double mean_p = pooled_sum[par] / pooled_n[par];
        double deficit = baseline - mean_p;
        if (deficit <= 0.0) continue;
        double sig = deficit / (sigma / sqrt(pooled_n[par]));
        if (sig >= 3.0) n_supporting++;
    }

    free(t); free(f);
    return n_supporting;
}

/* ───────────────────────────────────────────────────────────────────────────
 * §5 False-Positive Ruling (validator.py:661)
 * ─────────────────────────────────────────────────────────────────────────── */
void zs_false_positive_ruling(const ZSCandidate *c, const ZSDensityResult *den,
                              const ZSProbabilityResult *prob,
                              int n_transits, ZSFalsePositiveResult *out) {
    int v[12];
    int n_tests = 11;   /* FP-1..FP-9 (FP-5b, FP-5c included) */
    int n_extra = 0;    /* +1 when FP-10 is active */

    /* FP-1 BLS SNR > 5.5 */
    v[0] = (c->bls_snr > ZS_FP1_SNR_MIN) ? 1 : 0;
    /* FP-2 FAP < 0.05 | overridden by coherent evidence */
    v[1] = ((c->bls_fap < ZS_FP2_FAP_MAX) ||
            (c->coherent_evidence && c->bls_snr > ZS_FP1_SNR_MIN)) ? 1 : 0;
    /* FP-3 Even/Odd Δσ < 3.0 */
    v[2] = (c->even_odd_delta_sigma < ZS_FP3_EO_SIGMA_MAX) ? 1 : 0;
    /* FP-4 Shape ratio > 0.4 */
    v[3] = (c->shape_ratio > ZS_FP4_SHAPE_MIN) ? 1 : 0;
    /* FP-5 Secondary eclipse SNR < 3.0 */
    v[4] = (c->secondary_snr < ZS_FP5_SEC_SNR_MAX) ? 1 : 0;
    /* FP-5b Secondary/primary depth ratio < 0.30 */
    v[5] = (c->secondary_depth_ratio < ZS_FP5B_SEC_RATIO_MAX) ? 1 : 0;
    /* FP-5c alias band 0.20 < ratio < 0.90 → FAIL */
    v[6] = (c->alias_secondary_ratio > ZS_FP5C_ALIAS_LO &&
            c->alias_secondary_ratio < ZS_FP5C_ALIAS_HI) ? 0 : 1;
    /* FP-6 Centroid shift σ < 3.0 */
    v[7] = (c->centroid_sigma < ZS_FP6_CENTROID_MAX) ? 1 : 0;
    /* FP-7 Density ratio in [0.2, 5.0] */
    v[8] = (den->density_ratio > ZS_FP7_DENSITY_LO &&
            den->density_ratio < ZS_FP7_DENSITY_HI) ? 1 : 0;
    /* FP-8 Impact parameter < 0.9 (grazing → 0.95 value) */
    v[9] = ((prob->is_grazing ? 0.95 : 0.0) < ZS_FP8_IMPACT_MAX) ? 1 : 0;
    /* FP-9 catalog: C99 has no network — mirrors Python CATALOG_OFFLINE PASS */
    v[10] = 1;
    /* FP-10 min transits (only when light curve provided) */
    if (n_transits >= 0) {
        v[11] = (n_transits >= ZS_FP10_MIN_TRANSITS) ? 1 : 0;
        n_extra = 1;
    } else {
        v[11] = -1; /* skipped */
    }

    int n_pass = 0, n_fail = 0;
    int n_critical = 0, n_crit_pass = 0;
    /* weight: critical = FP1,2,3,5,5b,5c,7,9,10 ; major = FP4,6 ; moderate = FP8 */
    const int crit[] = {0, 1, 2, 4, 5, 6, 8, 10, 11};
    for (int i = 0; i < 12; i++) {
        if (v[i] < 0) continue;
        if (v[i]) n_pass++; else n_fail++;
    }
    for (size_t j = 0; j < sizeof(crit) / sizeof(crit[0]); j++) {
        int i = crit[j];
        if (v[i] < 0) continue;
        n_critical++;
        if (v[i]) n_crit_pass++;
    }

    int critical_passed = (n_crit_pass == n_critical);
    const char *overall;
    if (!critical_passed)
        overall = "FALSE_POSITIVE";
    else if (n_fail <= ZS_VERDICT_MAX_FAIL_PASS)
        overall = "SOVEREIGN_PASS";
    else if (n_fail <= ZS_VERDICT_MAX_FAIL_COND)
        overall = "CONDITIONAL_PASS";
    else
        overall = "FALSE_POSITIVE";

    out->n_tests = n_tests + n_extra;
    out->n_pass = n_pass;
    out->n_fail = n_fail;
    out->n_critical = n_critical;
    out->n_critical_pass = n_crit_pass;
    for (int i = 0; i < 12; i++) out->fp_verdicts[i] = v[i];
    snprintf(out->overall_verdict, sizeof(out->overall_verdict), "%s", overall);

    out->conflict_snr_density = 0;
    if (c->bls_snr > ZS_DENSITY_CONFLICT_SNR &&
        (den->density_ratio < ZS_DENSITY_CONFLICT_LO ||
         den->density_ratio > ZS_DENSITY_CONFLICT_HI))
        out->conflict_snr_density = 1;
    out->conflict_snr_shape = (c->bls_snr > 10.0 && c->shape_ratio <= 1.0) ? 1 : 0;
}

/* ───────────────────────────────────────────────────────────────────────────
 * Sovereign validation (build_full_proof equivalent, §1–§6)
 * ─────────────────────────────────────────────────────────────────────────── */
void zs_sovereign_validate(const ZSCandidate *c,
                           const double *time, size_t n_time,
                           const double *flux, size_t n_flux,
                           const double *flux_err, const double *model_flux,
                           ZSSovereignCard *out) {
    memset(out, 0, sizeof(*out));

    zs_compute_kepler(c, &out->kepler);
    zs_compute_geometry(c, &out->geometry);
    zs_compute_density(c, &out->kepler, &out->density);
    zs_compute_probability(c, &out->kepler, &out->probability);
    zs_compute_cvs(c, &out->cvs, out->cvs_verdict, sizeof(out->cvs_verdict));

    /* FP-10 when light curve + period present */
    int n_transits = -1;
    if (time && flux && n_time > 0 && n_flux > 0 && c->period_days > 0.0) {
        n_transits = zs_count_observed_transits(
            time, n_time, flux, n_flux,
            c->period_days, c->t0_days, c->transit_duration_hrs / 24.0);
    }
    out->n_transits = n_transits;

    zs_false_positive_ruling(c, &out->density, &out->probability, n_transits, &out->fp);

    /* §6 chi-squared (optional) */
    if (time && flux && flux_err && model_flux && n_time > 0) {
        double *diff = (double *)malloc(n_time * sizeof(double));
        double *resid = (double *)malloc(n_time * sizeof(double));
        double *rsq = (double *)malloc(n_time * sizeof(double));
        if (diff && resid && rsq) {
            for (size_t i = 0; i < n_time; i++)
                diff[i] = flux[i] - model_flux[i];
            for (size_t i = 0; i < n_time; i++)
                resid[i] = diff[i] / flux_err[i];
            for (size_t i = 0; i < n_time; i++)
                rsq[i] = resid[i] * resid[i];
            double chi2 = 0.0;
            for (size_t i = 0; i < n_time; i++) chi2 += rsq[i];
            double dof = fmax((double)n_time - 1.0, 1.0);
            out->chi_squared.chi2 = chi2;
            out->chi_squared.reduced_chi2 = chi2 / dof;
            out->chi_squared.dof = (int)dof;
        }
        free(diff); free(resid); free(rsq);
    }

    /* Aggregate verdict */
    int kepler_ok = out->kepler.verdict_pass;
    int geo_ok = out->geometry.verdict_pass;
    int den_ok = (out->density.verdict != 2);
    int prob_ok = out->probability.verdict_pass;
    int fp_ok = (strcmp(out->fp.overall_verdict, "FALSE_POSITIVE") != 0);

    int all_pass = kepler_ok && geo_ok && den_ok && prob_ok && fp_ok;
    out->all_sections_pass = all_pass;
    snprintf(out->sovereign_verdict, sizeof(out->sovereign_verdict), "%s",
             out->fp.overall_verdict);
}