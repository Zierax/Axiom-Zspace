"""
zspace_kernels.py  ·  C99-Version math kernels (Purce source) — v2 strict
==========================================================================
Module-level, numpy-pure mirror of the Axiom-ZSpace Sovereign Proof
Engine math.  INPUT to the Purce semantic compiler.

v2 REWRITE RULE (measured necessity, Purce v0.1.0):
  * ONE numpy call per line.
  * NO nested calls, NO arithmetic operators inside arguments.
  * Every compound constant (4π², 4/3, 1/3, 1+k, ...) is passed in as
    a FUNCTION PARAMETER, never inlined.
  * Rationale: Purce v0.1.0 corrupts compound expressions (observed:
    kepler_a_au_a3 emitted x/x; kepler_a_au_expr_3 emitted 1.0/3.0 in
    place of 4π²; density_rho_tic_gcc_expr_1 dropped π; generated .c
    signatures diverged from purce_src.h).  Single-op kernels with
    passed-in operands verified correct.
  * The C99 caller supplies the actual constant VALUES (same values as
    zspace_engine/constants.py).

Generated C99 carries provenance (.prov.json); every kernel MUST pass
the differential verifier (tests/verify_*.c vs numpy) before use.
"""

import numpy as np

# ── IAU 2015 constants (identical to zspace_engine/constants.py) ────────────
G_SI    = 6.67430e-11     # m^3 kg^-1 s^-2
M_SUN   = 1.9884e30       # kg
R_SUN   = 6.957e8         # m
AU      = 1.495978707e11  # m
R_EARTH = 6.3781e6        # m
PI      = 3.14159265358979323846


# ─────────────────────────────────────────────────────────────────────────────
# §1 Keplerian Dynamics (validator.py:239)
# ─────────────────────────────────────────────────────────────────────────────

def kepler_gm(G, M_star_solar):
    """gm_partial = G * M_star_solar."""
    return np.multiply(G, M_star_solar)


def kepler_gm_full(gm_partial, m_sun):
    """gm = gm_partial * m_sun."""
    return np.multiply(gm_partial, m_sun)


def kepler_p_sec_sq(P_sec):
    """P_sec^2."""
    return np.power(P_sec, 2.0)


def kepler_a3_num(gm, p_sec_sq):
    """a3 numerator: gm * P_sec^2."""
    return np.multiply(gm, p_sec_sq)


def kepler_pi_sq(pi):
    """π²."""
    return np.power(pi, 2.0)


def kepler_denominator(four, pi_sq):
    """4·π²."""
    return np.multiply(four, pi_sq)


def kepler_a3(a3_num, denom):
    """a3 = numerator / denominator."""
    return np.divide(a3_num, denom)


def kepler_a_m(a3, inv3):
    """a_m = a3^(1/3) (inv3 = 1.0/3.0 passed in)."""
    return np.power(a3, inv3)


def kepler_a_au(a_m, au):
    """a_au = a_m / AU."""
    return np.divide(a_m, au)


def kepler_p_sec_sq_2(P_sec):
    """P_sec^2 (residual path)."""
    return np.power(P_sec, 2.0)


def kepler_a_m_cubed(a_m):
    """a_m^3."""
    return np.power(a_m, 3.0)


def kepler_ratio_si(p_sec_sq, a_m_cubed):
    """P_sec^2 / a_m^3."""
    return np.divide(p_sec_sq, a_m_cubed)


def kepler_expected_si(four_pi_sq, gm):
    """4π² / GM."""
    return np.divide(four_pi_sq, gm)


def kepler_diff_si(ratio_si, expected_si):
    """ratio - expected."""
    return np.subtract(ratio_si, expected_si)


def kepler_diff_abs(diff):
    """|diff|."""
    return np.abs(diff)


def kepler_norm_diff_si(diff_abs, expected_si):
    """diff / expected."""
    return np.divide(diff_abs, expected_si)


def kepler_residual_si(norm_diff, hundred):
    """residual % = norm_diff * 100."""
    return np.multiply(norm_diff, hundred)


def kepler_p_yr(P_days, inv365):
    """P_yr = P_days / 365.25 (inv365 = 1/365.25 passed in)."""
    return np.multiply(P_days, inv365)


def kepler_p_yr_sq(P_yr):
    """P_yr^2."""
    return np.power(P_yr, 2.0)


def kepler_a_au_cubed(a_au):
    """a_au^3."""
    return np.power(a_au, 3.0)


def kepler_a_au_max(a_au_cubed, eps):
    """max(a³, eps)."""
    return np.maximum(a_au_cubed, eps)


def kepler_ratio_solar(p_yr_sq, a_au_max):
    """P_yr² / max(a³, eps)."""
    return np.divide(p_yr_sq, a_au_max)


def kepler_expect_solar(one, M_star_solar):
    """1 / M_star_solar (one = 1.0 passed in — Purce v0.1.0 flips the
    dividend when it is a literal constant: 1.0/x was emitted as x/1.0)."""
    return np.divide(one, M_star_solar)


def kepler_diff_solar(ratio, expect):
    """ratio - expect."""
    return np.subtract(ratio, expect)


def kepler_diff_solar_abs(diff):
    """|diff|."""
    return np.abs(diff)


def kepler_expect_max(expect, eps):
    """max(expect, eps)."""
    return np.maximum(expect, eps)


def kepler_norm_solar(diff_abs, expect_max):
    """diff / max(expect, eps)."""
    return np.divide(diff_abs, expect_max)


def kepler_residual_solar(norm_diff, hundred):
    """residual %."""
    return np.multiply(norm_diff, hundred)


# ─────────────────────────────────────────────────────────────────────────────
# §2 Geometric Consistency (validator.py:324)
# ─────────────────────────────────────────────────────────────────────────────

def geometry_delta_clamped(delta, eps):
    """max(delta, eps)."""
    return np.maximum(delta, eps)


def geometry_k(delta_clamped):
    """k = sqrt(delta)."""
    return np.sqrt(delta_clamped)


def geometry_k_sq(k):
    """k^2."""
    return np.power(k, 2.0)


def geometry_u1_term(u1, inv3):
    """u1/3."""
    return np.multiply(u1, inv3)


def geometry_u2_term(u2, inv6):
    """u2/6."""
    return np.multiply(u2, inv6)


def geometry_ld_sub(one, u1_term):
    """1 - u1/3."""
    return np.subtract(one, u1_term)


def geometry_ld_sub2(one_minus_u1, u2_term):
    """1 - u1/3 - u2/6."""
    return np.subtract(one_minus_u1, u2_term)


def geometry_ld_i_mean(one_minus_both, min_floor):
    """max(..., 0.1)."""
    return np.maximum(one_minus_both, min_floor)


def geometry_ld_zeropow_diff(one, mu_c):
    """1 - mu_c (mu_c = 1 for central transit)."""
    return np.subtract(one, mu_c)


def geometry_ld_zeropow(zeropow_diff, exponent):
    """(1 - mu_c)^exponent."""
    return np.power(zeropow_diff, exponent)


def geometry_ld_u1cen_term(u1, zeropow):
    """u1·(1-1)."""
    return np.multiply(u1, zeropow)


def geometry_ld_u2cen_term(u2, zeropow_sq):
    """u2·(1-1)²."""
    return np.multiply(u2, zeropow_sq)


def geometry_ld_cen_sub(one, u1cen_term):
    """1 - u1·(1-1)."""
    return np.subtract(one, u1cen_term)


def geometry_ld_cen_sub2(one_minus_u1, u2cen_term):
    """1 - u1·(1-1) - u2·(1-1)²."""
    return np.subtract(one_minus_u1, u2cen_term)


def geometry_ld_i_cen(one_minus_both, min_floor):
    """max(..., 0.1)."""
    return np.maximum(one_minus_both, min_floor)


def geometry_ld_numerator(k_sq, i_cen):
    """k² · i_cen."""
    return np.multiply(k_sq, i_cen)


def geometry_ld_corrected(numerator, i_mean):
    """k²·i_cen / i_mean."""
    return np.divide(numerator, i_mean)


def geometry_cons_resid_sub(delta, k_sq):
    """δ - k²."""
    return np.subtract(delta, k_sq)


def geometry_cons_resid_abs(sub):
    """|δ - k²|."""
    return np.abs(sub)


def geometry_cons_k_sq_max(k_sq, eps):
    """max(k², eps)."""
    return np.maximum(k_sq, eps)


def geometry_cons_resid_div(diff_abs, k_sq_max):
    """|δ - k²| / max(k², eps)."""
    return np.divide(diff_abs, k_sq_max)


def geometry_cons_resid_pct(residual, hundred):
    """residual %."""
    return np.multiply(residual, hundred)


def geometry_rp_m(k, R_star):
    """R_p [m] = k · R_★."""
    return np.multiply(k, R_star)


def geometry_rp_earth(rp_m, r_earth):
    """R_p / R_earth."""
    return np.divide(rp_m, r_earth)


# ─────────────────────────────────────────────────────────────────────────────
# §3 Stellar Density Constraint (validator.py:417)
# ─────────────────────────────────────────────────────────────────────────────

def density_one_plus_k(k, one):
    """1 + k."""
    return np.add(k, one)


def density_tdur_max(T_sec, one):
    """max(T_sec, 1)."""
    return np.maximum(T_sec, one)


def density_a_over_rs_num(P_sec, one_plus_k):
    """P_sec · (1+k)."""
    return np.multiply(P_sec, one_plus_k)


def density_a_over_rs_den(pi, tdur_max):
    """π · T_dur."""
    return np.multiply(pi, tdur_max)


def density_a_over_rs_transit(num, den):
    """a/R_★ = num / den."""
    return np.divide(num, den)


def density_rho_pow_p(P_sec):
    """P_sec^2."""
    return np.power(P_sec, 2.0)


def density_rho_den(G, p_sec_sq):
    """G · P_sec^2."""
    return np.multiply(G, p_sec_sq)


def density_rho_coeff(three_pi, rho_den):
    """3π / (G·P²)."""
    return np.divide(three_pi, rho_den)


def density_rho_pow_a(a_over_rs, three):
    """(a/R_★)³."""
    return np.power(a_over_rs, three)


def density_rho_transit_kgm3(coeff, rho_pow):
    """ρ [kg/m³] = coeff · pow."""
    return np.multiply(coeff, rho_pow)


def density_rho_transit_gcc(rho_kgm3, inv1000):
    """ρ [g/cm³]."""
    return np.multiply(rho_kgm3, inv1000)


def density_vol_coeff(four_thirds, pi):
    """(4/3)·π."""
    return np.multiply(four_thirds, pi)


def density_vol_pow(R_star):
    """R_★³."""
    return np.power(R_star, 3.0)


def density_vol(vol_coeff, R_star_cubed):
    """V = (4π/3)·R³."""
    return np.multiply(vol_coeff, R_star_cubed)


def density_rho_tic_kgm3(M_star, vol):
    """ρ_★ = M / V."""
    return np.divide(M_star, vol)


def density_rho_tic_gcc(rho_kgm3, inv1000):
    """ρ [g/cm³]."""
    return np.multiply(rho_kgm3, inv1000)


def density_rho_tic_max(rho_tic, eps):
    """max(ρ_tic, eps)."""
    return np.maximum(rho_tic, eps)


def density_ratio(rho_transit, rho_tic_max):
    """ρ_t / max(ρ_tic, eps)."""
    return np.divide(rho_transit, rho_tic_max)


def density_g_num(G, M_star):
    """G·M."""
    return np.multiply(G, M_star)


def density_g_den(R_star):
    """R_★²."""
    return np.power(R_star, 2.0)


def density_g_cgs_div(g_num, g_den):
    """G·M / R²."""
    return np.divide(g_num, g_den)


def density_g_cgs(g_cgs_div, hundred):
    """G·M/R² · 100."""
    return np.multiply(g_cgs_div, hundred)


def density_g_cgs_max(g_cgs, one):
    """max(g_cgs, 1)."""
    return np.maximum(g_cgs, one)


def density_logg(g_cgs_max):
    """log10(max(g_cgs, 1))."""
    return np.log10(g_cgs_max)


# ─────────────────────────────────────────────────────────────────────────────
# §4 Probability of Transit (validator.py:557)
# ─────────────────────────────────────────────────────────────────────────────

def probability_num(R_star, R_planet):
    """R_★ + R_p."""
    return np.add(R_star, R_planet)


def probability_den(a_m, num):
    """max(a, R_★+R_p)."""
    return np.maximum(a_m, num)


def probability_ptr(num, den):
    """P_tr = num/den."""
    return np.divide(num, den)


def probability_ptr_clamped(ptr, one):
    """min(P_tr, 1)."""
    return np.minimum(ptr, one)


def probability_rstar_max(R_star, eps):
    """max(R_★, eps)."""
    return np.maximum(R_star, eps)


def probability_a_over_rs(a_m, rstar_max):
    """a / max(R_★, eps)."""
    return np.divide(a_m, rstar_max)


def probability_p_max(P_sec, one):
    """max(P, 1)."""
    return np.maximum(P_sec, one)


def probability_sin_ratio(T_sec, p_max):
    """T / max(P, 1)."""
    return np.divide(T_sec, p_max)


def probability_sin_prod(pi, ratio):
    """π · T/P."""
    return np.multiply(pi, ratio)


def probability_sin_arg(sin_prod, pi):
    """min(π·T/P, π)."""
    return np.minimum(sin_prod, pi)


def probability_sin(sin_arg):
    """sin(π·T/P)."""
    return np.sin(sin_arg)


def probability_sin_term(a_over_rs, sin_val):
    """(a/R)·sin(...)."""
    return np.multiply(a_over_rs, sin_val)


def probability_one_plus_k(k, one):
    """1 + k."""
    return np.add(k, one)


def probability_opk_sq(one_plus_k):
    """(1+k)²."""
    return np.power(one_plus_k, 2.0)


def probability_st_sq(sin_term):
    """sin_term²."""
    return np.power(sin_term, 2.0)


def probability_b_sq(opk_sq, st_sq):
    """(1+k)² - sin_term²."""
    return np.subtract(opk_sq, st_sq)


def probability_b_clamped(b_sq):
    """max(b_sq, 0)."""
    return np.maximum(b_sq, 0.0)


def probability_b(b_clamped):
    """b = sqrt(...)."""
    return np.sqrt(b_clamped)


def probability_ingress_den(one_plus_k, eps):
    """max(1+k, eps)."""
    return np.maximum(one_plus_k, eps)


def probability_ingress_ratio(k, den):
    """k/(1+k)."""
    return np.divide(k, den)


def probability_ingress_hrs(T_hrs, ratio):
    """T_dur · k/(1+k)."""
    return np.multiply(T_hrs, ratio)


def probability_cos_max(num, a_m):
    """(R_★+R_p)/a."""
    return np.divide(num, a_m)


def probability_cos_clamped(cos_max, one):
    """min(cos_max, 1)."""
    return np.minimum(cos_max, one)


# NOTE: i_min_deg = degrees(arccos(cos_max)) is informational only
# (validator.py:602, not part of any gate).  np.arccos is NOT supported
# by Purce v0.1.0 (generator emits #error "composite"), so it is
# computed in the hand-written orchestration layer (zspace_core.c)
# with acos() from <math.h> — documented limitation.


# ─────────────────────────────────────────────────────────────────────────────
# CVS Composite Vitality Score (core.py:410)
# ─────────────────────────────────────────────────────────────────────────────

def cvs_wp_sp(w_p, s_p):
    """w_p·S_p."""
    return np.multiply(w_p, s_p)


def cvs_wd_sd(w_d, s_d):
    """w_d·S_d."""
    return np.multiply(w_d, s_d)


def cvs_wl_sl(w_l, s_l):
    """w_l·S_l."""
    return np.multiply(w_l, s_l)


def cvs_ws_ss(w_s, s_s):
    """w_s·S_s."""
    return np.multiply(w_s, s_s)


def cvs_num_a(wp_sp, wd_sd):
    """w_p·S_p + w_d·S_d."""
    return np.add(wp_sp, wd_sd)


def cvs_num_b(wl_sl, ws_ss):
    """w_l·S_l + w_s·S_s."""
    return np.add(wl_sl, ws_ss)


def cvs_numerator(num_a, num_b):
    """Σ w_i·S_i."""
    return np.add(num_a, num_b)


def cvs_den_a(w_p, w_d):
    """w_p + w_d."""
    return np.add(w_p, w_d)


def cvs_den_b(w_l, w_s):
    """w_l + w_s."""
    return np.add(w_l, w_s)


def cvs_total(den_a, den_b):
    """Σ w_i."""
    return np.add(den_a, den_b)


def cvs_total_max(total, eps):
    """max(Σw, eps)."""
    return np.maximum(total, eps)


def cvs_score(numerator, total_max):
    """CVS = num / max(Σw, eps)."""
    return np.divide(numerator, total_max)


# ─────────────────────────────────────────────────────────────────────────────
# §6 Chi-squared (chi_squared.py:261 quick_chi_squared)
# ─────────────────────────────────────────────────────────────────────────────

def chi_diff(flux, model_flux):
    """flux - model."""
    return np.subtract(flux, model_flux)


def chi_residuals(diff, flux_err):
    """(flux - model) / err."""
    return np.divide(diff, flux_err)


def chi_resid_sq(residuals):
    """residuals²."""
    return np.power(residuals, 2.0)


def chi_squared(resid_sq):
    """Σ residuals²."""
    return np.sum(resid_sq)


def chi_n_minus(n, one):
    """n - 1."""
    return np.subtract(n, one)


def chi_dof(n_minus, one):
    """max(n - 1, 1)."""
    return np.maximum(n_minus, one)


def chi_reduced(chi2, dof):
    """χ²_red = χ² / dof."""
    return np.divide(chi2, dof)


def chi_alt_dof(n, n_params):
    """n - n_params."""
    return np.subtract(n, n_params)


def chi_alt_dof_max(dof, one):
    """max(dof, 1)."""
    return np.maximum(dof, one)


def chi_reduced_alt(chi2, dof):
    """χ² / max(n - n_params, 1)."""
    return np.divide(chi2, dof)


# ─────────────────────────────────────────────────────────────────────────────
# Depth-consistency statistics (auditors.py:834)
# ─────────────────────────────────────────────────────────────────────────────

def depths_mean(depths):
    """μ = mean(depths)."""
    return np.mean(depths)


def depths_var_pop(depths):
    """Population variance (Purce reduce_var)."""
    return np.var(depths)


def depths_std_scale(var, n):
    """var * n."""
    return np.multiply(var, n)


def depths_std_scale2(scaled, inv_nm1):
    """var·n/(n-1) = sample variance (inv_nm1 = 1/(n-1) passed in)."""
    return np.multiply(scaled, inv_nm1)


def depths_std_sqrt(sample_var):
    """σ = sqrt(sample variance) — matches np.std(ddof=1)."""
    return np.sqrt(sample_var)


def depths_mu_max(mu, eps):
    """max(μ, eps)."""
    return np.maximum(mu, eps)


def depths_cv(std, mu_max):
    """CV = σ / max(μ, eps)."""
    return np.divide(std, mu_max)


def depths_resid(depths, mu):
    """d - μ."""
    return np.subtract(depths, mu)


def depths_resid_sq(resid):
    """(d - μ)²."""
    return np.power(resid, 2.0)


def depths_err_sq(errs):
    """err²."""
    return np.power(errs, 2.0)


def depths_ratio(resid_sq, err_sq):
    """(d - μ)² / err²."""
    return np.divide(resid_sq, err_sq)


def depths_chi2_reduced(ratio):
    """mean(...)."""
    return np.mean(ratio)


def depth_cons_chi_sub(chi2_red, one):
    """χ²_red - 1."""
    return np.subtract(chi2_red, one)


def depth_cons_norm_sub(norm, one):
    """norm - 1."""
    return np.subtract(norm, one)


def depth_cons_ratio(chi_sub, norm_sub):
    """(χ²_red - 1)/(norm - 1)."""
    return np.divide(chi_sub, norm_sub)


def depth_cons_sub(one, ratio):
    """1 - ratio."""
    return np.subtract(one, ratio)


def depth_consistency_s(sub, zero, one):
    """clip(sub, 0, 1)."""
    return np.clip(sub, zero, one)


# ─────────────────────────────────────────────────────────────────────────────
# Even/odd Welch statistics (auditors.py:798)
# ─────────────────────────────────────────────────────────────────────────────

def welch_sqrt_n(n):
    """√n."""
    return np.sqrt(n)


def welch_se(sigma, sqrt_n):
    """σ / √n."""
    return np.divide(sigma, sqrt_n)


def welch_delta(mu_e, mu_o):
    """μ_e - μ_o."""
    return np.subtract(mu_e, mu_o)


def welch_delta_abs(delta):
    """|μ_e - μ_o|."""
    return np.abs(delta)


def welch_se_e_sq(se_e):
    """se_e²."""
    return np.power(se_e, 2.0)


def welch_se_o_sq(se_o):
    """se_o²."""
    return np.power(se_o, 2.0)


def welch_var_sum(se_e_sq, se_o_sq):
    """se_e² + se_o²."""
    return np.add(se_e_sq, se_o_sq)


def welch_den_sqrt(var_sum):
    """sqrt(se_e² + se_o²)."""
    return np.sqrt(var_sum)


def welch_den_max(den_sqrt, eps):
    """max(sqrt(...), eps)."""
    return np.maximum(den_sqrt, eps)


def welch_delta_sigma(delta_abs, den_max):
    """Δσ = |Δ| / max(...)."""
    return np.divide(delta_abs, den_max)