#!/usr/bin/env python3
"""
run_pipeline.py  ·  Axiom-ZSpace Full Pipeline Runner (v1.0)
============================================================================
Orchestrates ingestion → detection → auditing → context → report.

V2.0 NEW Features:
  - Stellar Density Constraint (a/R★ mismatch < 20%)
  - Transit Shape V/U Discriminator (Ingress/Egress ratio)
  - Enhanced Centroid Offset Test (Pixel-Level TPF)
  - Multi-Sector Consistency Check
  - MCMC Posterior Validation
  - Local FITS Caching for faster re-analysis

Usage
-----
  # Run on a real TIC ID (requires internet + lightkurve):
  python run_pipeline.py --tic 260128333

  # Run synthetic self-test (no internet required):
  python run_pipeline.py --synthetic

  # Enable MCMC posterior validation (slower but more rigorous):
  python run_pipeline.py --tic 260128333 --mcmc

  # Enable TPF centroid analysis (downloads pixel data):
  python run_pipeline.py --tic 260128333 --tpf-centroids

  # Both at once:
  python run_pipeline.py --tic 260128333 --synthetic

Outputs
-------
  discovery_card_ZS-T-<TICID>-01.json   (Truthimatics Discovery Card)
"""

# CRITICAL: Import astroquery BEFORE any logging configuration
# This ensures astropy's custom logger class is initialized properly
# and prevents the "'Logger' object has no attribute '_set_defaults'" error
try:
    import astroquery
    from astropy import log as astropy_log
    astropy_log.setLevel('WARNING')
except ImportError:
    pass  # Libraries not installed yet

import argparse
import json
import sys
import traceback
import logging
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import yaml

# ── Engine imports ────────────────────────────────────────────────────────────
from zspace_engine.core      import VitalityMatrix, apply_hard_filters
from types import SimpleNamespace as _NS

from zspace_engine.auditors  import TransitAuditor
from zspace_engine.context   import StellarContextAuditor, StellarMetadata, ContextAuditResult
from zspace_engine.report    import TruthimaticsReport
from zspace_engine.logging_config import setup_logging, get_logger, suppress_astroquery_logger

# ── Logging setup ─────────────────────────────────────────────────────────────
# Initialize logging from production config
setup_logging()
logger = get_logger(__name__)

# Suppress verbose logging from astroquery and lightkurve
try:
    suppress_astroquery_logger()
    logging.getLogger('lightkurve').setLevel(logging.WARNING)
except Exception:
    pass

# ── Load production configuration ─────────────────────────────────────────────
def load_config(config_path: str = "config/production.yaml") -> dict:
    """Load production configuration from YAML file.
    
    Searches multiple locations: given path, config.yaml in repo root.
    """
    candidates = [config_path, "config/production.yaml", "config.yaml"]
    loaded = None
    for candidate in candidates:
        try:
            with open(candidate, 'r', encoding='utf-8') as f:
                loaded = yaml.safe_load(f)
                logger.info(f"Configuration loaded from {candidate}")
                break
        except FileNotFoundError:
            continue
        except Exception as e:
            logger.warning(f"Failed to load {candidate}: {e}")
            continue
    if loaded is not None:
        return loaded
    logger.warning(f"No configuration file found (searched: {candidates}), using defaults")
    return {
        "detection": {
            "bls_snr_threshold": 5.5,
            "fap_threshold": 1.0e-4,
            "cvs_planet_threshold": 0.80
        },
        "output": {
            "base_directory": "axiom_output"
        },
        "fp_filters": {
            "density_mismatch_threshold": 0.20,
            "ingress_fraction_threshold": 0.45,
            "run_mcmc": False,
            "use_tpf_centroids": False,
            "check_multi_sector": False,
        }
    }

# Load configuration at module level
CONFIG = load_config()


# ─────────────────────────────────────────────────────────────────────────────
# Synthetic transit generator (for self-testing without MAST)
# ─────────────────────────────────────────────────────────────────────────────

def generate_synthetic_transit(
    n_days:       float = 27.0,
    cadence_min:  float = 2.0,
    period_days:  float = 3.7,
    depth:        float = 0.009,    # ~0.9% — Neptune-size around G-dwarf
    duration_hrs: float = 2.0,
    t0_offset:    float = 1.2,
    noise_ppm:    float = 300.0,
    seed:         int   = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Inject a box-shaped transit into Gaussian noise.
    Returns (time_days, flux_normalised).
    """
    rng = np.random.default_rng(seed)
    cadence_days = cadence_min / 1440.0
    time = np.arange(0, n_days, cadence_days)
    flux = np.ones_like(time)

    # Inject transits
    half_dur = (duration_hrs / 24.0) / 2.0
    phase    = ((time - t0_offset) / period_days) % 1.0
    phase[phase > 0.5] -= 1.0
    in_transit = np.abs(phase) <= half_dur / period_days
    flux[in_transit] -= depth

    # Add Gaussian noise
    noise_frac = noise_ppm * 1e-6
    flux += rng.normal(0, noise_frac, size=len(flux))

    return time, flux


# ─────────────────────────────────────────────────────────────────────────────
# Core pipeline function (V2.0 — all FP filters integrated)
# ─────────────────────────────────────────────────────────────────────────────

def run_pipeline(
    tic_id:         str,
    lc_product:     _NS,
    planet_order:   int = 1,
    fetch_tic_meta: bool = True,
    run_mcmc:       bool = False,
    use_tpf_centroids: bool = False,
    check_multi_sector: bool = False,
    engine:         str = "python",
) -> dict:
    """
    Execute the full Axiom-ZSpace pipeline on a pre-processed LightCurveProduct.

    Parameters
    ----------
    tic_id          : TIC identifier string
    lc_product      : output of LightCurveIngester.process() or from_arrays()
    planet_order    : planet index for ZSpace ID (default 1)
    fetch_tic_meta  : whether to query MAST for stellar parameters
    run_mcmc        : whether to run MCMC posterior validation
    use_tpf_centroids : whether to download TPF for pixel-level centroids
    check_multi_sector : whether to check transit across multiple sectors

    Returns
    -------
    Discovery Card dictionary (also prints JSON to stdout)
    """
    logger.info("="*70)
    logger.info(f"AXIOM-ZSPACE PIPELINE V2.0  |  TIC {tic_id}  |  Planet {planet_order:02d}")
    logger.info("="*70)

    # Read FP filter config
    fp_config = CONFIG.get("fp_filters", {})
    run_mcmc = run_mcmc or fp_config.get("run_mcmc", False)
    use_tpf_centroids = use_tpf_centroids or fp_config.get("use_tpf_centroids", False)
    check_multi_sector = check_multi_sector or fp_config.get("check_multi_sector", False)

    # ── Phase 1: BLS Detection ────────────────────────────────────────────────
    logger.info("[PHASE 1] BLS Signal Detection ...")
    
    # Get detection thresholds from config
    snr_threshold = CONFIG.get("detection", {}).get("bls_snr_threshold", 5.5)
    fap_threshold = CONFIG.get("detection", {}).get("fap_threshold", 1.0e-4)
    
    if engine == "c99":
        # C99 BLS engine (zspace_card bls, OpenMP-parallel) — ~100x faster
        from c99_bridge import run_c99_bls
        c99bls = run_c99_bls(
            lc_product.time, lc_product.flux_flat,
            period_min=0.5, period_max=13.5,
        )
        bls = _NS(
            period_best=float(c99bls["period_days"]),
            transit_depth=abs(float(c99bls["depth"])),
            transit_duration=float(c99bls["duration_hrs"]) / 24.0,
            t0=float(c99bls["t0_days"]),
            snr=float(c99bls["snr"]),
            fap=float(c99bls["fap"]),
            n_trial_periods=2000,
            s_periodicity=0.0,
            proof="C99 BLS engine (zspace_card bls; OpenMP; "
                  "log-likelihood objective matching astropy)",
            flags=["C99_ENGINE"],
            bls_power_max=float(c99bls["power"]),
            snr_threshold=snr_threshold,
            fap_threshold=fap_threshold,
        )
        bls.passed_detection_gate = lambda: (
            bls.snr > bls.snr_threshold and bls.fap < bls.fap_threshold)
    else:
        from zspace_engine.detectors import BLSDetector
        detector = BLSDetector(
            period_min=0.5, 
            period_max=13.5,
            snr_threshold=snr_threshold,
            fap_threshold=fap_threshold
        )
        bls      = detector.run(lc_product.time, lc_product.flux_flat)

    logger.info(f"  Period:  {bls.period_best:.5f} d")
    logger.info(f"  SNR:     {bls.snr:.2f}  [threshold: {snr_threshold}]")
    logger.info(f"  FAP:     {bls.fap:.2e}  [threshold: {fap_threshold:.0e}]")
    logger.info(f"  Gate:    {'PASS' if bls.passed_detection_gate() else 'FAIL'}")

    if not bls.passed_detection_gate():
        logger.warning("[ABORT] BLS detection gate failed.  Not a credible periodic signal.")
        logger.warning("Emitting partial Discovery Card with FALSE_POSITIVE verdict.")

    # Re-flatten with period hint
    if bls.passed_detection_gate():
        if engine == "c99":
            from c99_bridge import run_c99_flatten
            flux_flat2 = run_c99_flatten(
                lc_product.time, lc_product.flux_norm, bls.period_best
            )
        else:
            from zspace_engine.ingestion import LightCurveIngester
            ingester   = LightCurveIngester(tic_id=tic_id)
            flux_flat2, _, _ = ingester._savgol_flatten(
                lc_product.time, lc_product.flux_norm, bls.period_best
            )
        time_f = lc_product.time
        flux_f = flux_flat2
    else:
        time_f = lc_product.time
        flux_f = lc_product.flux_flat

    # Phase-fold for auditors (Python path; C99 audits fold internally)
    bin_phase = bin_flux = None
    if engine == "python" or run_mcmc:
        from zspace_engine.detectors import BLSDetector
        bin_phase, bin_flux, _ = BLSDetector.fold_and_bin(
            time_f, flux_f, bls.period_best, bls.t0, n_bins=200
        )

    # ── Phase 2: Transit Auditing (Enhanced) ──────────────────────────────────
    logger.info("[PHASE 2] Transit Vitality Auditing ...")
    auditor = TransitAuditor(run_mcmc=run_mcmc)

    if engine == "c99":
        from c99_bridge import run_c99_audit
        c99aud = run_c99_audit(
            time_f, flux_f, bls.period_best, bls.t0,
            bls.transit_duration * 24.0, bls.transit_depth)
        eo_c = c99aud["even_odd"]
        dc_c = c99aud["depth_consistency"]
        ie_c = c99aud["ingress_egress"]
        eo_result = _NS(
            delta_sigma=float(eo_c["delta_sigma"]),
            is_eb_flag=bool(eo_c["is_eb_flag"]),
            t_stat=float(eo_c["t_stat"]),
            p_value=float(eo_c["p_value"]),
            n_even=int(eo_c["n_even"]),
            n_odd=int(eo_c["n_odd"]),
            depth_even=float(eo_c["depth_even"]),
            depth_odd=float(eo_c["depth_odd"]),
            proof="C99 audit engine (zspace_card audit)",
            flags=["C99_ENGINE"],
        )
        depth_result = _NS(
            n_transits=int(dc_c["n_transits"]),
            mean_depth=float(dc_c["mean_depth"]),
            std_depth=float(dc_c["std_depth"]),
            cv=float(dc_c["cv"]),
            sigma_med=float(dc_c["sigma_med"]),
            chi2_red=float(dc_c["chi2_red"]),
            s_depth=float(dc_c["s_depth"]),
            depths=[0.0] * int(dc_c["n_transits"]),
            proof="C99 audit engine (zspace_card audit)",
            flags=["C99_ENGINE"],
        )
        ie_result = _NS(
            depth_fit=float(ie_c["depth_fit"]),
            ingress_fraction=float(ie_c["ingress_fraction"]),
            flat_fraction=float(ie_c["flat_fraction"]),
            ingress_hrs=float(ie_c["ingress_hrs"]),
            flat_hrs=float(ie_c["flat_hrs"]),
            is_v_shape=bool(ie_c["is_v_shape"]),
            fp_risk=str(ie_c["fp_risk"]),
            fit_ok=bool(ie_c["fit_ok"]),
            proof="C99 audit engine (zspace_card audit)",
            flags=["C99_ENGINE"],
        )
    else:
        eo_result = auditor.even_odd_test(
            time_f, flux_f, bls.period_best, bls.t0, bls.transit_duration
        )
        depth_result = auditor.depth_consistency_score(
            time_f, flux_f, bls.period_best, bls.t0, bls.transit_duration, eo_result
        )
        ie_result = auditor.ingress_egress_test(
            bin_phase, bin_flux, bls.period_best, bls.transit_duration, bls.transit_depth
        )

    logger.info(f"  Even/Odd Delta-sigma:    {eo_result.delta_sigma:.3f}  "
          f"[EB flag: {'YES' if eo_result.is_eb_flag else 'NO'}]")

    logger.info(f"  Depth CV:        {depth_result.cv:.4f}  ->  S_delta = {depth_result.s_depth:.4f}")

    limb_result = auditor.limb_shape_score(
        period=bls.period_best, duration=bls.transit_duration,
        transit_depth=bls.transit_depth,
        time=time_f, flux=flux_f, t0=bls.t0,
    )
    logger.info(f"  Shape ratio:     {limb_result.shape_ratio:.3f}  ->  S_tau = {limb_result.s_limb:.4f}")

    # ── Phase 2.5: NEW — Ingress/Egress V/U Shape Test ────────────────────────
    logger.info("[PHASE 2.5] Ingress/Egress Shape Analysis ...")
    logger.info(f"  Ingress fraction:  {ie_result.ingress_fraction:.3f}  "
          f"[V-shape: {'YES' if ie_result.is_v_shape else 'NO'}]  "
          f"[FP risk: {ie_result.fp_risk}]")
    
    # Penalize limb score if V-shape detected
    if ie_result.is_v_shape:
        limb_result.s_limb *= 0.5
        limb_result.flags.append(f"V_SHAPE_PENALTY | S_τ halved | ingress_frac={ie_result.ingress_fraction:.3f}")
        logger.info(f"  -> S_tau penalized to {limb_result.s_limb:.4f} (V-shape)")

    # ── Phase 2.7: NEW — MCMC Posterior Validation ────────────────────────────
    if run_mcmc:
        # Adaptive MCMC step count based on signal quality
        if ie_result.is_v_shape or ie_result.fp_risk == "HIGH":
            mcmc_steps = 50   # Fast-reject: obvious FP candidates
            logger.info("[PHASE 2.7] MCMC Posterior Validation (FAST: 50 steps, V-shape/high FP risk) ...")
        elif bls.snr > 30 and not ie_result.is_v_shape:
            mcmc_steps = 500  # Full analysis for clean high-SNR signals
            logger.info("[PHASE 2.7] MCMC Posterior Validation (FULL: 500 steps, high SNR clean signal) ...")
        else:
            mcmc_steps = 200  # Standard analysis
            logger.info("[PHASE 2.7] MCMC Posterior Validation (STANDARD: 200 steps) ...")
        mcmc_result = auditor.mcmc_validate(
            bin_phase, bin_flux, bls.period_best, bls.transit_duration, bls.transit_depth,
            n_steps=mcmc_steps,
        )
        logger.info(f"  Gaussian:    {'YES' if mcmc_result.is_gaussian else 'NO'}  "
              f"[skew_max: {mcmc_result.skewness_max:.2f}, kurt_max: {mcmc_result.kurtosis_max:.2f}]")
        if not mcmc_result.is_gaussian:
            # Penalize depth score if posteriors are non-Gaussian
            depth_result.s_depth *= 0.7
            depth_result.flags.append("MCMC_NOISE | S_δ × 0.7 due to non-Gaussian posteriors")
            logger.info(f"  -> S_delta penalized to {depth_result.s_depth:.4f} (non-Gaussian posteriors)")
    else:
        mcmc_result = None

    # ── Phase 3: Stellar Context (Enhanced) ──────────────────────────────────
    logger.info("[PHASE 3] Stellar Context Auditing ...")
    ctx_auditor = StellarContextAuditor(
        fetch_tic=fetch_tic_meta,
        use_tpf_centroids=use_tpf_centroids,
        check_multi_sector=check_multi_sector,
    )
    
    # Pass a/R★ from limb shape fit for density constraint check
    a_rs_transit = limb_result.a_rs if limb_result.a_rs > 0 else 0.0
    
    ctx_result = ctx_auditor.audit(
        tic_id, time_f, flux_f, bls.period_best, bls.t0, bls.transit_duration,
        a_rs_transit=a_rs_transit,
        sector=lc_product.sector,
    )
    logger.info(f"  TIC source:      {ctx_result.metadata.source}  [version: {ctx_result.metadata.tic_version}]")
    logger.info(f"  Centroid shift:  {ctx_result.centroid.centroid_shift_sigma:.3f}sigma  "
          f"[{'FLAGGED' if ctx_result.centroid.is_flagged else 'OK'}]  "
          f"[method: {ctx_result.centroid.method}]")
    logger.info(f"  Secondary SNR:   {ctx_result.secondary.snr_at_half_phase:.3f}  "
          f"[{'FLAGGED' if ctx_result.secondary.is_flagged else 'OK'}]")
    
    # NEW: Density constraint
    logger.info(f"  Density check:   deviation={ctx_result.density_check.fractional_deviation*100:.1f}%  "
          f"[{'FLAGGED' if ctx_result.density_check.is_flagged else 'OK'}]")
    
    # NEW: Multi-sector consistency
    logger.info(f"  Multi-sector:    {ctx_result.multi_sector.n_sectors_consistent}/"
          f"{ctx_result.multi_sector.n_sectors_available} consistent  "
          f"[boost: ×{ctx_result.multi_sector.confidence_boost:.1f}]")
    
    logger.info(f"  S_S =            {ctx_result.s_stellar:.4f}")

    # ── Phase 4: Vitality Matrix + Orbital Mechanics ──────────────────────────
    logger.info("[PHASE 4] Compositing Vitality Matrix ...")
    matrix = VitalityMatrix(tic_id=tic_id, planet_order=planet_order)

    # Periodicity score: 0 if gate failed
    s_p      = bls.s_periodicity if bls.passed_detection_gate() else 0.0
    proof_p  = bls.proof

    matrix.ingest_scores(
        s_periodicity = s_p,       proof_p = proof_p,
        s_depth       = depth_result.s_depth, proof_d = depth_result.proof,
        s_limb        = limb_result.s_limb,   proof_l = limb_result.proof,
        s_stellar     = ctx_result.s_stellar,  proof_s = ctx_result.proof,
        flags_p       = bls.flags,
        flags_d       = depth_result.flags,
        flags_l       = limb_result.flags,
        flags_s       = ctx_result.flags,
    )

    matrix.compute_orbital_mechanics(
        period_days          = bls.period_best,
        transit_depth        = bls.transit_depth,
        stellar_mass_solar   = ctx_result.metadata.stellar_mass_solar,
        stellar_teff         = ctx_result.metadata.stellar_teff_k,
        stellar_radius_solar = ctx_result.metadata.stellar_radius_solar,
    )

    # ── Phase 3.5: Hard Physical Filters (pre-CVS rejection gates) ───────────
    logger.info("[PHASE 3.5] Hard Physical Filters ...")
    hard_filter = apply_hard_filters(
        planet_radius_earth  = matrix.orbital.planet_radius_earth,
        transit_depth        = bls.transit_depth,
        secondary_snr        = ctx_result.secondary.snr_at_half_phase,
        density_deviation    = ctx_result.density_check.fractional_deviation,
        even_odd_sigma       = eo_result.delta_sigma if not np.isnan(eo_result.delta_sigma) else 0.0,
        is_v_shape           = ie_result.is_v_shape,
        transit_duration_hrs = bls.transit_duration * 24.0,
        period_days          = bls.period_best,
    )
    if not hard_filter.passed:
        logger.warning(f"  HARD REJECT: {hard_filter.rejection}")
        logger.warning(f"  Proof: {hard_filter.proof}")
    elif hard_filter.cvs_penalty < 1.0:
        logger.info(f"  Penalties applied: {'; '.join(hard_filter.flags)}")
    else:
        logger.info("  All hard physical filters PASSED")

    # Apply hard filter to CVS engine
    matrix.cvs_engine.apply_hard_filter(hard_filter)

    # ── Phase 3.6: Critical-FP Veto Gate ─────────────────────────────────────
    # Any critically-flagged FP indicator forces CVS below the ambiguous
    # threshold (FALSE POSITIVE) regardless of component scores. This is
    # single-channel: the same indicators already lower S_S via the Stellar
    # Context score; the veto only prevents a PERFECT other score from
    # overriding a critical stellar/binary red flag.
    veto_reasons = []
    if not np.isnan(ctx_result.centroid.centroid_shift_sigma) and \
            ctx_result.centroid.centroid_shift_sigma > 5.0:
        veto_reasons.append(
            f"CENTROID_CRITICAL | σ={ctx_result.centroid.centroid_shift_sigma:.2f} > 5.0"
        )
    if ctx_result.secondary.is_flagged:
        veto_reasons.append(
            f"SECONDARY_ECLIPSE | SNR={ctx_result.secondary.snr_at_half_phase:.2f}"
        )
    if eo_result.is_eb_flag:
        veto_reasons.append(
            f"EVEN_ODD_EB | Welch t={eo_result.t_stat:.2f}, p={eo_result.p_value:.4f}"
        )
    for reason in veto_reasons:
        logger.warning(f"  VETO: {reason}")
        matrix.cvs_engine.apply_veto(reason)
    if veto_reasons:
        logger.warning(f"  {len(veto_reasons)} critical FP veto(s) → CVS forced to FALSE POSITIVE")

    cvs = matrix.cvs_engine.compute()
    logger.info("")
    logger.info("  +---------------------------------------------+")
    logger.info(f"  |  CVS = {cvs:.4f}   ->   {matrix.cvs_engine.verdict:<28}|")
    logger.info("  +---------------------------------------------+")

    # ── Phase 5: Report ───────────────────────────────────────────────────────
    logger.info("[PHASE 5] Emitting Truthimatics Discovery Card ...")
    reporter = TruthimaticsReport(matrix)
    reporter.attach_ingestion(lc_product)
    reporter.attach_bls(bls)
    reporter.attach_audits(eo_result, depth_result, limb_result)
    reporter.attach_context(ctx_result)

    card = reporter.emit()
    
    # ── Phase 4.5: C99 Sovereign Engine (--engine c99) ──────────────────────
    card["c99_sovereign_card"] = None
    if engine == "c99":
        logger.info("[PHASE 4.5] C99 Sovereign Engine (bin/zspace_card) ...")
        try:
            from c99_bridge import run_c99_sovereign

            candidate = {
                "period_days": float(bls.period_best),
                "transit_depth": float(bls.transit_depth),
                "transit_duration_hrs": float(bls.transit_duration * 24.0),
                "t0_days": float(bls.t0),
                "stellar_mass_solar": float(getattr(ctx_result.metadata, "stellar_mass_solar", 1.0) or 1.0),
                "stellar_radius_solar": float(getattr(ctx_result.metadata, "stellar_radius_solar", 1.0) or 1.0),
                "stellar_teff_k": float(getattr(ctx_result.metadata, "stellar_teff_k", 5772.0) or 5772.0),
                "stellar_logg": float(getattr(ctx_result.metadata, "stellar_logg", 4.44) or 4.44),
                "planet_radius_earth": float(matrix.orbital.planet_radius_earth),
                "bls_snr": float(bls.snr),
                "bls_fap": float(bls.fap),
                "even_odd_delta_sigma": float(eo_result.delta_sigma) if not np.isnan(eo_result.delta_sigma) else 0.0,
                "shape_ratio": float(limb_result.shape_ratio),
                "secondary_snr": float(ctx_result.secondary.snr_at_half_phase),
                "secondary_depth_ratio": 0.0,
                "alias_secondary_ratio": 0.0,
                "coherent_evidence": 0,
                "centroid_sigma": float(ctx_result.centroid.centroid_shift_sigma)
                if not np.isnan(ctx_result.centroid.centroid_shift_sigma) else 0.0,
                "limb_dark_u1": 0.45,
                "limb_dark_u2": 0.15,
                "s_periodicity": float(s_p),
                "s_depth": float(depth_result.s_depth),
                "s_limb": float(limb_result.s_limb),
                "s_stellar": float(ctx_result.s_stellar),
            }
            sov = run_c99_sovereign(candidate, time=time_f, flux=flux_f)
            card["c99_sovereign_card"] = sov
            logger.info(f"  C99 sovereign verdict: {sov.get('sovereign_verdict', '?')}  "
                        f"(kepler a={sov['section_1_kepler']['a_au']:.4f} AU, "
                        f"P_tr={sov['section_4_probability']['P_tr']:.4f}, "
                        f"FP {sov['section_5_fp_ruling']['n_pass']}/"
                        f"{sov['section_5_fp_ruling']['n_tests']})")
        except Exception as e:
            logger.warning(f"  C99 engine unavailable, keeping Python verdict: {e}")
            card["c99_sovereign_card"] = {"error": str(e)}
    
    # Inject V2.0 FP filter results into card
    card["fp_filters_v2"] = {
        "ingress_egress_test": {
            "ingress_fraction": ie_result.ingress_fraction,
            "flat_fraction": ie_result.flat_fraction,
            "is_v_shape": ie_result.is_v_shape,
            "fp_risk": ie_result.fp_risk,
            "proof": ie_result.proof,
            "flags": ie_result.flags,
        },
        "density_constraint": {
            "a_rs_transit": ctx_result.density_check.a_rs_transit,
            "a_rs_catalog": ctx_result.density_check.a_rs_catalog,
            "deviation_percent": round(ctx_result.density_check.fractional_deviation * 100, 2),
            "is_flagged": ctx_result.density_check.is_flagged,
            "proof": ctx_result.density_check.proof,
        },
        "multi_sector_consistency": {
            "n_sectors_available": ctx_result.multi_sector.n_sectors_available,
            "n_sectors_consistent": ctx_result.multi_sector.n_sectors_consistent,
            "sectors_checked": ctx_result.multi_sector.sectors_checked,
            "confidence_boost": ctx_result.multi_sector.confidence_boost,
            "is_consistent": ctx_result.multi_sector.is_consistent,
            "proof": ctx_result.multi_sector.proof,
        },
        "centroid_method": ctx_result.centroid.method,
        "hard_physical_filters": {
            "passed": hard_filter.passed,
            "rejection": hard_filter.rejection,
            "cvs_penalty": hard_filter.cvs_penalty,
            "flags": hard_filter.flags,
            "proof": hard_filter.proof,
        },
    }
    
    if mcmc_result is not None:
        card["fp_filters_v2"]["mcmc_validation"] = {
            "n_walkers": mcmc_result.n_walkers,
            "n_steps": mcmc_result.n_steps,
            "acceptance_fraction": mcmc_result.acceptance_fraction,
            "is_gaussian": mcmc_result.is_gaussian,
            "skewness_max": mcmc_result.skewness_max,
            "kurtosis_max": mcmc_result.kurtosis_max,
            "proof": mcmc_result.proof,
            "flags": mcmc_result.flags,
        }

    # Save to disk
    filename = f"discovery_card_{card['zspace_id']}.json"
    reporter.save(card, filename)
    logger.info(f"  ZSpace ID:  {card['zspace_id']}")
    logger.info(f"  Verdict:    {card['verdict']}")
    logger.info(f"  Saved ->     {filename}")

    return card


# ─────────────────────────────────────────────────────────────────────────────
# Synthetic self-test
# ─────────────────────────────────────────────────────────────────────────────

def run_synthetic_test(engine: str = "python") -> dict:
    """
    Generate a synthetic transit, inject it into the pipeline, and verify
    that the engine recovers the correct period and emits a valid Discovery Card.
    """
    logger.info("="*70)
    logger.info("SYNTHETIC SELF-TEST  |  Injected Neptune-class transit")
    logger.info("="*70)

    TRUE_PERIOD  = 3.7       # days
    TRUE_DEPTH   = 0.009     # ~0.9%
    TRUE_DUR_HRS = 2.0
    TRUE_T0      = 1.2

    time, flux = generate_synthetic_transit(
        n_days       = 27.0,
        cadence_min  = 2.0,
        period_days  = TRUE_PERIOD,
        depth        = TRUE_DEPTH,
        duration_hrs = TRUE_DUR_HRS,
        t0_offset    = TRUE_T0,
        noise_ppm    = 300.0,
    )
    logger.info(f"Injected: P={TRUE_PERIOD} d, depth={TRUE_DEPTH*1e6:.0f} ppm, dur={TRUE_DUR_HRS} h")
    logger.info(f"Data:     {len(time)} cadences over {time[-1]:.1f} days  (2-min cadence, 300 ppm noise)")

    if engine == "c99":
        from c99_bridge import run_c99_flatten
        good = np.isfinite(flux)
        time_c = time[good]
        flux_c = flux[good]
        n_drop = int((~good).sum())
        med = float(np.median(flux_c))
        flux_norm = flux_c / med
        flux_flat = np.asarray(run_c99_flatten(time_c, flux_norm, TRUE_PERIOD))
        lc = _NS(
            tic_id="SYNTHETIC",
            sector=0,
            time=time_c,
            flux_flat=flux_flat,
            flux_norm=flux_norm,
            n_points_cleaned=len(time_c),
            n_dropped_quality=0,
            n_dropped_sigma=n_drop,
            fits_source="local_array_c99",
        )
    else:
        from zspace_engine.ingestion import LightCurveIngester
        lc = LightCurveIngester.from_arrays(
            tic_id  = "SYNTHETIC",
            time    = time,
            flux    = flux,
            quality = None,
            period_hint_days = TRUE_PERIOD,
        )

    card = run_pipeline(
        tic_id        = "SYNTHETIC",
        lc_product    = lc,
        planet_order  = 1,
        fetch_tic_meta = False,   # no MAST for synthetic
        engine        = engine,
    )

    # Verify period recovery
    recovered_period = card["bls_detection"]["period_days"]
    period_err_pct   = abs(recovered_period - TRUE_PERIOD) / TRUE_PERIOD * 100.0
    logger.info("Period recovery check:")
    logger.info(f"  True:      {TRUE_PERIOD:.5f} d")
    logger.info(f"  Recovered: {recovered_period:.5f} d")
    logger.info(f"  Error:     {period_err_pct:.3f}%  {'PASS' if period_err_pct < 5 else 'FAIL'}")

    return card


# ─────────────────────────────────────────────────────────────────────────────
# Real TIC pipeline
# ─────────────────────────────────────────────────────────────────────────────

def run_real_tic(
    tic_id: str,
    run_mcmc: bool = False,
    use_tpf_centroids: bool = False,
    check_multi_sector: bool = False,
    engine: str = "python",
) -> dict:
    """
    Download and process a real TIC target from MAST.
    """
    logger.info(f"Downloading TIC {tic_id} from MAST ...")
    from zspace_engine.ingestion import LightCurveIngester
    ingester   = LightCurveIngester(tic_id=tic_id, mission="TESS", exptime="short", use_cache=True)
    lc_product = ingester.process()
    logger.info(f"Cleaned:  {lc_product.n_points_cleaned} cadences  "
          f"(dropped {lc_product.n_dropped_quality} quality + {lc_product.n_dropped_sigma} sigma)")
    logger.info(f"Source:   {lc_product.fits_source}")

    return run_pipeline(
        tic_id        = tic_id,
        lc_product    = lc_product,
        planet_order  = 1,
        fetch_tic_meta = True,
        run_mcmc       = run_mcmc,
        use_tpf_centroids = use_tpf_centroids,
        check_multi_sector = check_multi_sector,
        engine         = engine,
    )


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_sector_range(sectors_str: str) -> list[int]:
    """
    Parse sector range string into list of sector numbers.
    
    Supports:
    - Single sector: "42" -> [42]
    - Range: "1-50" -> [1, 2, 3, ..., 50]
    - Comma-separated: "1,3,5" -> [1, 3, 5]
    - Mixed: "1,3,5-10,15" -> [1, 3, 5, 6, 7, 8, 9, 10, 15]
    """
    sectors = set()
    parts = sectors_str.split(',')
    
    for part in parts:
        part = part.strip()
        if '-' in part:
            try:
                start, end = part.split('-')
                start = int(start.strip())
                end = int(end.strip())
                if start > end:
                    raise ValueError(f"Invalid range: {part} (start > end)")
                if start < 1:
                    raise ValueError(f"Sector numbers must be >= 1, got {start}")
                sectors.update(range(start, end + 1))
            except ValueError as e:
                if "invalid literal" in str(e):
                    raise ValueError(f"Invalid range format: {part}. Expected format: 'start-end'")
                raise
        else:
            try:
                sector = int(part)
                if sector < 1:
                    raise ValueError(f"Sector numbers must be >= 1, got {sector}")
                sectors.add(sector)
            except ValueError:
                raise ValueError(f"Invalid sector number: {part}")
    
    return sorted(sectors)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Axiom-ZSpace Deterministic Exoplanet Detection Engine (v1.0)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run synthetic self-test (no internet required)
  python run_pipeline.py --synthetic
  
  # Process a single TIC target
  python run_pipeline.py --tic 260128333
  
  # Process with all FP filters enabled
  python run_pipeline.py --tic 260128333 --mcmc --tpf-centroids --multi-sector
  
  # Process entire TESS sector
  python run_pipeline.py --sector 42
  
  # Process multiple sectors
  python run_pipeline.py --sectors 1-50
        """,
    )
    parser.add_argument("--tic",       type=str, default=None,  
                        help="TIC ID to process (e.g., 260128333)")
    parser.add_argument("--synthetic", action="store_true",     
                        help="Run synthetic self-test with injected transit")
    parser.add_argument("--sector",    type=int, default=None,  
                        help="Process entire TESS sector")
    parser.add_argument("--sectors",   type=str, default=None,  
                        help="Process multiple TESS sectors (e.g., '1-50')")
    parser.add_argument("--max-targets", type=int, default=None,
                        help="Max targets per sector (for testing)")
    parser.add_argument("--output",    type=str, default=None,   
                        help="Base output directory")
    
    # V2.0 FP filter flags
    parser.add_argument("--mcmc", action="store_true",
                        help="Enable MCMC posterior validation (slower, more rigorous)")
    parser.add_argument("--tpf-centroids", action="store_true",
                        help="Enable TPF pixel-level centroid analysis")
    parser.add_argument("--multi-sector", action="store_true",
                        help="Enable multi-sector consistency checking")
    parser.add_argument("--engine", type=str, choices=["python", "c99"], default="python",
                        help="Proof engine: 'python' (default) or 'c99' "
                             "(C99-Version/bin/zspace_card, Purce-generated kernels)")
    
    args = parser.parse_args()

    if not args.tic and not args.synthetic and not args.sector and not args.sectors:
        parser.print_help()
        logger.info("No target specified.  Running synthetic self-test by default.")
        args.synthetic = True

    # Determine output directory
    if args.output is None:
        if args.sector or args.sectors:
            args.output = CONFIG.get("output", {}).get("base_directory", "axiom_output")
        else:
            args.output = "."
    
    Path(args.output).mkdir(parents=True, exist_ok=True)

    cards = []

    # Handle multiple sectors processing
    if args.sectors:
        try:
            from zspace_engine.sector_processor import SectorProcessor
            
            sector_list = parse_sector_range(args.sectors)
            
            logger.info("="*70)
            logger.info(f"MULTI-SECTOR PROCESSING  |  {len(sector_list)} sectors: {sector_list[0]}-{sector_list[-1]}")
            logger.info("="*70)
            
            total_stats = {
                'sectors_processed': 0,
                'sectors_failed': 0,
                'total_targets': 0,
                'total_processed': 0,
                'total_new_discoveries': 0,
                'total_known_planets': 0,
                'total_false_positives': 0,
                'total_failed': 0,
                'sector_summaries': []
            }
            
            for sector_num in sector_list:
                try:
                    logger.info(f"\nPROCESSING SECTOR {sector_num} ({total_stats['sectors_processed']+1}/{len(sector_list)})")
                    
                    processor = SectorProcessor(
                        sector=sector_num,
                        output_dir=args.output,
                        config=CONFIG,
                        max_targets=args.max_targets
                    )
                    summary = processor.process_sector()
                    
                    total_stats['sectors_processed'] += 1
                    total_stats['total_targets'] += summary['total_targets']
                    total_stats['total_processed'] += summary['processed']
                    total_stats['total_new_discoveries'] += summary['new_discoveries']
                    total_stats['total_known_planets'] += summary['known_planets']
                    total_stats['total_false_positives'] += summary['false_positives']
                    total_stats['total_failed'] += summary['failed']
                    total_stats['sector_summaries'].append({
                        'sector': sector_num,
                        'new_discoveries': summary['new_discoveries'],
                        'known_planets': summary['known_planets'],
                        'false_positives': summary['false_positives']
                    })
                    
                except Exception as e:
                    total_stats['sectors_failed'] += 1
                    logger.error(f"Sector {sector_num} failed: {e}", exc_info=True)
                    continue
            
            # Save summary
            overall_summary_path = Path(args.output) / "multi_sector_summary.json"
            with open(overall_summary_path, 'w') as f:
                json.dump(total_stats, f, indent=2)
            
            logger.info(f"\nTotal new discoveries: {total_stats['total_new_discoveries']}")
            return
            
        except ValueError as e:
            logger.error(f"Invalid sector range: {e}")
            sys.exit(1)

    # Handle single sector processing
    if args.sector:
        try:
            from zspace_engine.sector_processor import SectorProcessor
            
            processor = SectorProcessor(
                sector=args.sector,
                output_dir=args.output,
                config=CONFIG,
                max_targets=args.max_targets
            )
            summary = processor.process_sector()
            
            logger.info(f"Sector {args.sector}: {summary['new_discoveries']} discoveries")
            return
            
        except Exception:
            logger.error(f"Sector {args.sector} failed:", exc_info=True)
            sys.exit(1)

    if args.synthetic:
        try:
            card = run_synthetic_test(engine=args.engine)
            cards.append(card)
        except Exception:
            logger.error("Synthetic test failed:", exc_info=True)
            sys.exit(1)

    if args.tic:
        try:
            card = run_real_tic(
                args.tic,
                run_mcmc=args.mcmc,
                use_tpf_centroids=args.tpf_centroids,
                check_multi_sector=args.multi_sector,
                engine=args.engine,
            )
            cards.append(card)
        except Exception:
            logger.error(f"Pipeline failed for TIC {args.tic}:", exc_info=True)
            sys.exit(1)

    logger.info("="*70)
    logger.info(f"PIPELINE COMPLETE  |  {len(cards)} Discovery Card(s) emitted")
    logger.info("="*70)


if __name__ == "__main__":
    main()
