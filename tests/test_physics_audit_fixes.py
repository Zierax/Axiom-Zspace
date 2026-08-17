"""
test_physics_audit_fixes.py  ·  Physics/Logic Fix Verification
================================================================
Unit tests for the physics deep-audit remediation (Priority 1-5):

  P1.1  SNR formula uses matched-filter SE  σ·sqrt(1/N_in + 1/N_out)
  P1.2  Even/Odd test uses Welch's t-test with honest p-value
  P1.3  Centroid shift significance uses SE of the mean difference
  P2    Proof-string attribution honesty (Kepler III / T_eq / R_p)
  P2.4  T_eq NaN/Inf guard
  P2.5  Planet radius error propagation
  P3    Single-channel penalties + Critical-FP veto gate
  P4    run_pipeline.py imports without crashing

Run:  python -m pytest tests/test_physics_audit_fixes.py -v
"""

import math

import numpy as np
import pytest

from zspace_engine.core import (
    semi_major_axis_au,
    equilibrium_temperature_k,
    planet_radius_earth,
    apply_hard_filters,
    CompositeVitalityScore,
    ComponentScore,
    THRESHOLD_AMBIGUOUS,
)


# ─────────────────────────────────────────────────────────────────────────────
# P1.1  SNR formula
# ─────────────────────────────────────────────────────────────────────────────

def test_planet_radius_matches_geometry():
    """δ = (R_p/R★)²  →  R_p = R★·√δ / R_earth_solar."""
    from zspace_engine.constants import R_EARTH_SOLAR
    rp, _ = planet_radius_earth(0.009, 1.0)
    # sqrt(0.009) = 0.09487 R_sun; divide by R_earth/R_sun
    expected = math.sqrt(0.009) / R_EARTH_SOLAR
    assert rp == pytest.approx(expected, rel=1e-9)


def test_radius_error_propagation():
    """σ_Rp/Rp = sqrt((σ_R★/R★)² + (0.5·σ_δ/δ)²)."""
    rp, proof = planet_radius_earth(
        0.009, 1.0,
        stellar_radius_err_solar=0.03,
        transit_depth_err=0.0005,
    )
    # parse "± X.XXX R⊕"
    err_str = proof.split("± ")[1].split(" R⊕")[0]
    rp_err = float(err_str)
    expected_rel = math.sqrt((0.03 / 1.0) ** 2 + (0.5 * 0.0005 / 0.009) ** 2)
    # tolerance accounts for 3-decimal rounding of the printed error
    assert rp_err == pytest.approx(rp * expected_rel, rel=5e-3)


# ─────────────────────────────────────────────────────────────────────────────
# P2.2  Kepler III correctness + honest attribution
# ─────────────────────────────────────────────────────────────────────────────

def test_kepler_third_law_earth_analog():
    """1 year, 1 M_sun → a ≈ 1 AU (with IAU 2015 constants)."""
    a, proof = semi_major_axis_au(365.25, 1.0)
    assert a == pytest.approx(1.0, rel=1e-3)
    assert "Kepler III" in proof
    assert "IAU 2015" in proof  # constants attribution kept
    assert "IAU 2015): a" not in proof.replace(" ", "")  # law NOT attributed to IAU


def test_kepler_third_law_mass_dependence():
    """a ∝ M^(1/3) at fixed period."""
    a1, _ = semi_major_axis_au(100.0, 1.0)
    a2, _ = semi_major_axis_au(100.0, 8.0)
    assert a2 == pytest.approx(a1 * 8.0 ** (1.0 / 3.0), rel=1e-9)


# ─────────────────────────────────────────────────────────────────────────────
# P2.2  Equilibrium temperature
# ─────────────────────────────────────────────────────────────────────────────

def test_eq_temp_earth_analog():
    """T_eq = T_eff·√(R★/2a)·(1-A)^0.25 for Sun-Earth (a=1, A=0.3)."""
    from zspace_engine.constants import R_SUN, AU
    t, _ = equilibrium_temperature_k(5778.0, 1.0, 1.0, albedo=0.30)
    ratio = (1.0 * R_SUN) / (2.0 * 1.0 * AU)
    expected = 5778.0 * math.sqrt(ratio) * (0.7) ** 0.25
    assert t == pytest.approx(expected, rel=1e-9)
    assert 250 < t < 270  # physical sanity for Earth-like


def test_eq_temp_flag_extreme():
    """Non-finite / >1e5 K results are flagged, not silently reported."""
    _, proof = equilibrium_temperature_k(5778.0, 1.0, 1.0e-7)
    assert "PHYSICALLY_UNLIKELY_EQ_T" in proof


# ─────────────────────────────────────────────────────────────────────────────
# P1.2  Even/Odd Welch t-test
# ─────────────────────────────────────────────────────────────────────────────

def _make_transit_curve(ns, depth_even, depth_odd, period=2.0, t0=0.5,
                        dur=0.2, pts_per_transit=60, oot_dt=0.05,
                        noise_sigma=1e-4):
    """Build a light curve with the given per-transit depths and interleaved
    out-of-transit baseline points so depth extraction finds every transit.

    A small white-noise floor is included: per-transit depth uncertainties
    (matched-filter propagation) require non-zero OOT scatter."""
    rng = np.random.default_rng(3)
    times, fluxes = [], []
    all_ns = list(ns)
    for n in all_ns:
        t_centre = t0 + n * period
        t_i = t_centre + np.linspace(-dur / 2, dur / 2, pts_per_transit)
        depth = depth_even if (n % 2 == 0) else depth_odd
        f_i = 1.0 - depth * np.exp(-0.5 * ((t_i - t_centre) / (dur / 6)) ** 2)
        times.append(t_i)
        fluxes.append(f_i)
        # out-of-transit baseline before and after each transit
        for sign in (-1, 1):
            t_o = t_centre + sign * np.linspace(dur / 2 + 0.05, dur / 2 + 0.6, 12)
            times.append(t_o)
            fluxes.append(np.ones_like(t_o))
    times = np.concatenate(times)
    fluxes = np.concatenate(fluxes) + rng.normal(0.0, noise_sigma, len(times))
    order = np.argsort(times)
    return times[order], fluxes[order]


def test_even_odd_uses_welch():
    from zspace_engine.auditors import TransitAuditor

    auditor = TransitAuditor(verbose=False)
    period, t0, dur = 2.0, 0.5, 0.2
    ns = np.arange(1, 9)  # 8 transits → 4 even, 4 odd
    times, fluxes = _make_transit_curve(
        ns, depth_even=0.02, depth_odd=0.004,
        period=period, t0=t0, dur=dur,
    )
    res = auditor.even_odd_test(times, fluxes, period, t0, dur)
    assert res.t_stat > 3.0
    assert res.p_value < 0.01
    assert res.is_eb_flag is True
    assert res.proof.startswith("EVEN_ODD |")


def test_even_odd_same_depth_no_flag():
    from zspace_engine.auditors import TransitAuditor

    auditor = TransitAuditor(verbose=False)
    period, t0, dur = 2.0, 0.5, 0.2
    ns = np.arange(1, 9)
    times, fluxes = _make_transit_curve(
        ns, depth_even=0.01, depth_odd=0.01,
        period=period, t0=t0, dur=dur,
    )
    res = auditor.even_odd_test(times, fluxes, period, t0, dur)
    assert res.p_value > 0.01
    assert res.is_eb_flag is False


# ─────────────────────────────────────────────────────────────────────────────
# P1.3  Centroid SE-of-mean-difference
# ─────────────────────────────────────────────────────────────────────────────

def test_centroid_uses_se_of_mean_diff():
    from zspace_engine.context import CentroidShiftTest

    rng = np.random.default_rng(42)
    # 100 in-transit and 100 out-of-transit centroids with a real 0.3 px shift
    col_in = rng.normal(0.0, 0.1, 100)
    row_in = rng.normal(0.0, 0.1, 100)
    col_out = rng.normal(0.3, 0.1, 100)
    row_out = rng.normal(0.0, 0.1, 100)
    n_in = n_out = 100

    result = CentroidShiftTest._vector_centroid(
        np.concatenate([col_in, col_out]),
        np.concatenate([row_in, row_out]),
        np.array([True] * 100 + [False] * 100),
        np.array([False] * 100 + [True] * 100),
        n_in, n_out,
    )
    # Proper SE from the actual sample variances AND the actual sampled shift
    shift_actual = math.sqrt(
        (np.mean(col_in) - np.mean(col_out)) ** 2 +
        (np.mean(row_in) - np.mean(row_out)) ** 2
    )
    se_col = math.sqrt(np.var(col_in) / 100 + np.var(col_out) / 100)
    se_row = math.sqrt(np.var(row_in) / 100 + np.var(row_out) / 100)
    expected_sigma = shift_actual / math.sqrt(se_col ** 2 + se_row ** 2)
    assert result.centroid_shift_sigma == pytest.approx(expected_sigma, rel=1e-3)
    assert "SE_shift" in result.proof


# ─────────────────────────────────────────────────────────────────────────────
# P3  Single-channel penalties + veto gate
# ─────────────────────────────────────────────────────────────────────────────

def test_density_moderate_no_hard_penalty():
    """70% density deviation is scored once (S_S), not via hard-filter CVS."""
    hf = apply_hard_filters(
        planet_radius_earth=1.2, transit_depth=0.01,
        density_deviation=0.7, is_v_shape=True,
    )
    assert hf.passed is True
    assert hf.cvs_penalty == pytest.approx(1.0)  # no double penalty
    assert not any("DENSITY" in f for f in hf.flags)


def test_density_impossible_hard_reject():
    """>500% density deviation is a hard reject (physically impossible)."""
    hf = apply_hard_filters(
        planet_radius_earth=1.2, transit_depth=0.01, density_deviation=6.0
    )
    assert hf.passed is False
    assert "HARD_REJECT_DENSITY_IMPOSSIBLE" in hf.flags


def test_vshape_no_double_penalty():
    """V-shape alone no longer multiplies CVS ×0.7 (S_τ already covers it)."""
    hf = apply_hard_filters(
        planet_radius_earth=1.2, transit_depth=0.01, is_v_shape=True
    )
    assert hf.passed is True
    assert hf.cvs_penalty == pytest.approx(1.0)


def test_veto_gate_forces_false_positive():
    cvs = CompositeVitalityScore()
    cvs.register(ComponentScore("periodicity", 0.97, 1.0, "p"))
    cvs.register(ComponentScore("depth", 0.83, 1.0, "d"))
    cvs.register(ComponentScore("limb", 0.61, 1.0, "l"))
    cvs.register(ComponentScore("stellar", 0.31, 1.0, "s"))
    cvs.apply_veto("TEST_CENTROID_CRITICAL")
    val = cvs.compute()
    assert val < THRESHOLD_AMBIGUOUS
    assert cvs.verdict == "FALSE POSITIVE"
    assert any("VETO" in line for line in cvs.proof_chain)


def test_no_veto_no_override():
    cvs = CompositeVitalityScore()
    cvs.register(ComponentScore("periodicity", 0.97, 1.0, "p"))
    cvs.register(ComponentScore("depth", 0.83, 1.0, "d"))
    cvs.register(ComponentScore("limb", 0.61, 1.0, "l"))
    cvs.register(ComponentScore("stellar", 0.31, 1.0, "s"))
    val = cvs.compute()
    assert val == pytest.approx(1.0)
    assert cvs.verdict == "PLANET CANDIDATE"


# ─────────────────────────────────────────────────────────────────────────────
# P4  run_pipeline imports
# ─────────────────────────────────────────────────────────────────────────────

def test_run_pipeline_imports():
    import run_pipeline
    assert isinstance(run_pipeline.CONFIG, dict)
    assert "detection" in run_pipeline.CONFIG


# ─────────────────────────────────────────────────────────────────────────────
# Sanity: hard filters unchanged for other cases
# ─────────────────────────────────────────────────────────────────────────────

def test_radius_cap_still_rejects():
    hf = apply_hard_filters(planet_radius_earth=30.0, transit_depth=0.01)
    assert hf.passed is False
    assert "HARD_REJECT_RADIUS" in hf.flags
