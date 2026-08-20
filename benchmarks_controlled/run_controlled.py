#!/usr/bin/env python3
"""
run_controlled.py  ·  Honest Controlled Evaluation Runner
===========================================================
Runs the FULL local Axiom-ZSpace pipeline (ingestion → BLS → auditors → CVS
→ sovereign validator) on a controlled synthetic dataset, then computes
honest detection/validation metrics:

  * Recall@period   : % of injected planets whose signal is found AND folded
                      at the correct ephemeris (≤5% period error).
  * Detection recall: % detected at ANY period (before ephem check).
  * FPR (false-positive rate): % of contamination targets wrongly certified
                      as NEW_DISCOVERY by the sovereign validator.
  * Precision / F1 over the certified-discovery decision at period level.

The validator archive query is stubbed to OFFLINE deterministically (no
network, no ground-truth leakage): the sovereign verdict then reflects ONLY
the pipeline's intrinsic physics/QA gating.

Usage:
  python benchmarks_controlled/run_controlled.py --true 80 --false 80 \
      --out benchmarks_controlled/runs/RUNNAME
"""

import argparse
import json
import logging
import math
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from zspace_engine.logging_config import setup_logging, get_logger
from zspace_engine.ingestion import LightCurveIngester
from zspace_engine.detectors import BLSDetector
from zspace_engine.auditors import TransitAuditor
from zspace_engine.validator import AxiomValidator
from zspace_engine import thresholds as T
from zspace_engine.config import fap_threshold, cvs_planet_threshold, snr_threshold

setup_logging()
logger = get_logger(__name__)

CERTIFIED = ("NEW_DISCOVERY", "OFFLINE_NEW_DISCOVERY")

# Positive-evidence override for the strict power-FAP (FP-2) firewall. When
# ENABLED, a candidate with strong repeated-observation evidence can override
# a marginal power-FAP. Measured on this dataset this is UNSAFE: pure-noise
# folds reach SNR 6-8 and dipF≈0.94 (every depth-based scalar fails), so the
# override certifies noise (FPR explosion). The strict power-FAP is therefore
# the honest firewall. Controlled via config/production.yaml -> thresholds.
ENABLE_COHERENT_OVERRIDE = bool(T.threshold("coherent_override_enabled"))
COHERENT_OVERRIDE_MIN_SNR = float(T.threshold("coherent_min_snr"))
COHERENT_OVERRIDE_MIN_TRANSITS = int(T.threshold("coherent_min_transits"))
COHERENT_OVERRIDE_MIN_DIP_FRACTION = float(T.threshold("coherent_min_dip_fraction"))


# ─────────────────────────────────────────────────────────────────────────────
# Offline validator stub: deterministic OFFLINE, no ground-truth leak.
# ─────────────────────────────────────────────────────────────────────────────

class _OfflineQuerier:
    def query(self, tic_id):
        return ([], [], "CONTROLLED_OFFLINE")


def _patch_validator_offline(validator: AxiomValidator) -> None:
    validator._querier = _OfflineQuerier()


# ─────────────────────────────────────────────────────────────────────────────
# Pipeline re-run for one target (local arrays, blind search)
# ─────────────────────────────────────────────────────────────────────────────

def validate_candidate(
    tic_id: str,
    cand,
    time_arr: np.ndarray,
    flux_arr: np.ndarray,
    stellar: dict,
    run_dir: Path,
    return_unvalidated: bool = False,
    engine: str = "python",
):
    """Run auditors + sovereign validator for one BLS candidate.

    Returns a marker tuple:
      (candidate, result_dict_or_None)
    When `return_unvalidated`, the candidate is returned unchanged (used to
    detect whether a candidate would certify if it failed the FP-2 FAP gate).
    """
    auditor = TransitAuditor(run_mcmc=False)
    eo = auditor.even_odd_test(time_arr, flux_arr, cand.period_best, cand.t0, cand.transit_duration)
    dc = auditor.depth_consistency_score(time_arr, flux_arr, cand.period_best, cand.t0, cand.transit_duration, eo)
    limb = auditor.limb_shape_score(
        period=cand.period_best, duration=cand.transit_duration,
        transit_depth=cand.transit_depth,
        time=time_arr, flux=flux_arr, t0=cand.t0,
    )
    sec = auditor.secondary_eclipse_test(
        time_arr, flux_arr, cand.period_best, cand.t0, cand.transit_duration,
    )

    # Positive-evidence override for the power-FAP null test. The strict
    # power-FAP (red-noise conservative) is the FPR firewall, but a genuine but
    # SHALLOW transit sits in the red-noise tail and over-fires it. A candidate
    # with overwhelming REPEATED-OBSERVATION evidence — >=3 independent transits
    # whose measured depths are consistent (low CV), coherent matched-filter
    # SNR, and no phase-0.5 eclipse — is certified on that positive evidence,
    # not on a weakened null test. Single events (1 transit) and ephemeral
    # noise folds fail this and stay rejected by the power-FAP.
    n_occ = int(len(dc.depths)) if dc.depths is not None else 0
    # Dip-recurrence: the fraction of independent transits that actually show a
    # dip. A genuine periodic transit dips on nearly every epoch; ephemeral
    # noise dips recur on ~50% of epochs randomly. This is a robust
    # positive-evidence discriminator even at very low SNR (where depth CV is
    # noise-inflated and unusable).
    _deps = np.asarray(dc.depths, dtype=float).ravel() if dc.depths is not None else np.array([])
    dip_fraction = float(np.mean(_deps > 0.0)) if _deps.size > 0 else 0.0
    if ENABLE_COHERENT_OVERRIDE:
        coherent_evidence = int(
            n_occ >= 3
            and dip_fraction >= 0.6
            and cand.snr >= 6.5
            and sec.secondary_ratio < 0.30
        )
    else:
        coherent_evidence = 0

    # Alias-secondary check (FP-5c): a grazing EB whose BLS period lands on a
    # sub-harmonic (found = true_period / 2) folds BOTH eclipses onto phase 0,
    # hiding the phase-0.5 eclipse from the FP-5/FP-5b windows — the candidate
    # looks like a clean planet. Re-fold at 2× and 3× the candidate period: at
    # the true (integer-multiple) period the EB's secondary reappears as an
    # ASYMMETRIC eclipse (0.2 < ratio < 0.9) with high SNR, while a REAL planet
    # folds to SYMMETRIC equal-depth transits (ratio ~ 1.0) or a sub-threshold
    # SNR. The old FP-5c over-vetoed because it flagged ANY phase-0.5 dip at the
    # doubled period; a planet ALSO shows a 0.5 dip (it is a genuine second
    # transit) but with ratio -> 1.0.
    ALIAS_SEC_MIN_SNR = float(T.threshold("fp5c_alias_min_snr"))
    _alias_ratio = 0.0
    for _mult in (2.0, 3.0):
        _r = auditor.secondary_eclipse_test(
            time_arr, flux_arr, cand.period_best * _mult, cand.t0, cand.transit_duration,
        )
        if 0.20 < _r.secondary_ratio < 0.90 and _r.secondary_snr >= ALIAS_SEC_MIN_SNR:
            _alias_ratio = max(_alias_ratio, _r.secondary_ratio)
    _alias_secondary = _alias_ratio

    s_p = cand.s_periodicity if cand.passed_detection_gate() else 0.0
    w_p, w_d, w_l, w_s = (
        float(T.threshold("cvs_w_periodicity")), float(T.threshold("cvs_w_depth")),
        float(T.threshold("cvs_w_limb")), float(T.threshold("cvs_w_secondary")),
    )
    cvs = (
        w_p * s_p +
        w_d * dc.s_depth +
        w_l * limb.s_limb +
        w_s * 0.5
    ) / (w_p + w_d + w_l + w_s)
    cvs_verdict = "PLANET CANDIDATE" if cvs > cvs_planet_threshold() else "UNCERTAIN"

    from zspace_engine.constants import R_EARTH_SOLAR
    st_rad = stellar.get("st_rad", 0.6) or 0.6
    bls_rp = st_rad * math.sqrt(abs(cand.transit_depth)) / R_EARTH_SOLAR

    validator = AxiomValidator(output_dir=str(run_dir / "validation"), verbose=False)
    _patch_validator_offline(validator)

    eo_sigma = eo.delta_sigma if not np.isnan(eo.delta_sigma) else 0.0
    shape_r = limb.shape_ratio if not np.isnan(limb.shape_ratio) else 1.0

    vresult = validator.validate(
        tic_id=tic_id,
        period_days=cand.period_best,
        transit_depth=cand.transit_depth,
        transit_duration_hrs=cand.transit_duration * 24,
        t0_btjd=cand.t0,
        stellar_mass_solar=stellar.get("st_mass", 0.6) or 0.6,
        stellar_radius_solar=st_rad,
        stellar_teff_k=stellar.get("st_teff", 3900.0) or 3900.0,
        stellar_logg=stellar.get("st_logg", 4.66) or 4.66,
        planet_radius_earth=bls_rp,
        cvs_score=cvs,
        cvs_verdict=cvs_verdict,
        bls_snr=cand.snr,
        bls_fap=cand.fap_power,           # strict power-spectrum FAP → FP-2 firewall
        even_odd_delta_sigma=eo_sigma,
        shape_ratio=shape_r,
        secondary_snr=sec.secondary_snr,
        secondary_depth_ratio=sec.secondary_ratio,
        alias_secondary_ratio=_alias_secondary,
        coherent_evidence=coherent_evidence,
        time=time_arr,
        flux=flux_arr,
    )

    result = {
        "validation_status": vresult.status,
        "cvs": cvs,
        "cvs_verdict": cvs_verdict,
        "snr": cand.snr,
        "fap_power": cand.fap_power,
        "fap_snr": cand.fap_snr,
        "secondary_ratio": sec.secondary_ratio,
        "secondary_snr": sec.secondary_snr,
        "coherent_evidence": coherent_evidence,
        "n_transits": n_occ,
        "dip_fraction": round(dip_fraction, 3),
        "alias_secondary_ratio": round(_alias_secondary, 3),
    }
    if vresult.output_file and Path(vresult.output_file).exists():
        try:
            card = json.loads(Path(vresult.output_file).read_text())
            result["sovereign_verdict"] = card.get("sovereign_verdict")
        except Exception:
            pass

    if engine == "c99":
        from c99_bridge import run_c99_sovereign
        cand_dict = {
            "period_days": cand.period_best,
            "transit_depth": cand.transit_depth,
            "transit_duration_hrs": cand.transit_duration * 24,
            "t0_days": cand.t0,
            "stellar_mass_solar": stellar.get("st_mass", 0.6) or 0.6,
            "stellar_radius_solar": stellar.get("st_rad", 0.6) or 0.6,
            "stellar_teff_k": stellar.get("st_teff", 3900.0) or 3900.0,
            "stellar_logg": stellar.get("st_logg", 4.66) or 4.66,
            "planet_radius_earth": bls_rp,
            "bls_snr": cand.snr,
            "bls_fap": cand.fap_power,
            "even_odd_delta_sigma": eo_sigma,
            "shape_ratio": shape_r,
            "secondary_snr": sec.secondary_snr,
            "secondary_depth_ratio": sec.secondary_ratio,
            "alias_secondary_ratio": _alias_secondary,
            "coherent_evidence": coherent_evidence,
            "centroid_sigma": 0.0,
            "limb_dark_u1": 0.45,
            "limb_dark_u2": 0.15,
            "s_periodicity": s_p,
            "s_depth": float(dc.s_depth),
            "s_limb": float(limb.s_limb),
            "s_stellar": 0.5,
        }
        try:
            c99_card = run_c99_sovereign(cand_dict, time_arr, flux_arr)
            sv = c99_card.get("sovereign_verdict")
            result["sovereign_verdict"] = sv
            result["sovereign_verdict_c99"] = sv
            result["validation_status"] = (
                "OFFLINE_NEW_DISCOVERY"
                if sv in ("SOVEREIGN_PASS", "CONDITIONAL_PASS") else "FALSE_POSITIVE"
            )
            result["c99_error"] = None
        except Exception as e:
            result["sovereign_verdict_c99"] = None
            result["c99_error"] = str(e)
    return cand, result


def evaluate_target(target, run_dir: Path, engine: str = "python") -> Dict:
    """Run the full pipeline on a SyntheticTarget. Blind: no period prior."""
    tic_id = target.tic_id
    res = {
        "kind": target.kind, "subkind": target.subkind, "tic_id": tic_id,
        "label_period": target.label_period, "detected": False,
        "period_found": None, "period_error_pct": None, "snr": None,
        "cvs": None, "sovereign_verdict": None, "validation_status": None,
        "ladder_rank": None, "n_candidates_tested": 0,
        "sovereign_verdict_c99": None, "c99_error": None,
        "error": None, "processing_time_sec": 0.0,
        "white_ppm": target.meta.get("white_ppm"),
        "target_snr": target.meta.get("target_snr"),
        "injected_depth": target.injected_depth,
    }

    start = time.time()
    try:
        lc = LightCurveIngester.from_arrays(
            tic_id=tic_id,
            time=target.time,
            flux=target.flux,
            quality=target.quality,
            period_hint_days=None,
        )
        time_arr = lc.time
        flux_arr = lc.flux_flat

        baseline_days = float(time_arr[-1] - time_arr[0])
        max_search = min(T.threshold("detector_period_max"), baseline_days / 3.0)

        # Weak-target recovery: a fine frequency ladder on the raw (unbinned)
        # light curve. Binning to 30k destroyed the shallow true peaks.
        detector = BLSDetector(
            period_min=T.threshold("detector_period_min"), period_max=max(0.6, max_search),
            snr_threshold=snr_threshold(), fap_threshold=fap_threshold(),
            frequency_factor=float(T.threshold("detector_frequency_factor")),
        )
        bls = detector.run(time=time_arr, flux=flux_arr)   # global max (fallback)
        ladder = detector.top_candidates(
            time=time_arr, flux=flux_arr, k=int(T.threshold("ladder_k")),
            min_relative_snr=float(T.threshold("ladder_min_relative_snr")),
        )

        # Detection is decided from the CANDIDATE LADDER, not the global max:
        # for weak targets the unconditional maximum is often a low-SNR noise
        # spike while the true transit sits lower in the ladder at high SNR.
        if not ladder:
            ladder = [bls] if bls else []
        elif not any(abs(math.log(c.period_best / bls.period_best)) < 0.05 for c in ladder):
            ladder = [bls] + ladder

        if not ladder or not any(c.snr > float(T.threshold("no_detection_snr_floor")) for c in ladder):
            best_snr = max((c.snr for c in ladder), default=0.0)
            res["validation_status"] = "NO_DETECTION"
            res["error"] = f"No BLS detection (best ladder SNR={best_snr:.1f})"
            return res

        res["detected"] = True
        res["snr"] = max(c.snr for c in ladder)

        certified = None
        first_status = None
        n_tested = 0
        from zspace_engine.ephemeris import EphemerisResolver
        resolver = EphemerisResolver()
        for cand in ladder:
            n_tested += 1
            _, cand_res = validate_candidate(tic_id, cand, time_arr, flux_arr, target.stellar, run_dir, engine=engine)
            if first_status is None:
                first_status = cand_res["validation_status"]
            if cand_res["validation_status"] in CERTIFIED and certified is None:
                certified = cand

                # ── Alias/epoch resolution (Stage A): a sub-harmonic alias
                # (P_true/2, P_true/3) certifies cleanly in ITS OWN fold but is
                # the same physical signal. Test the 2P/3P folds; if the
                # candidate is an integer alias, RE-VALIDATE the fundamental
                # period with the full sovereign validator and adopt the true
                # ephemeris only if that second validation certifies.
                from zspace_engine.ephemeris import EphemerisResolver
                ephem = resolver.resolve(
                    time_arr, flux_arr, cand.period_best, cand.t0,
                    period_min=0.5, period_max=max(0.6, max_search),
                )
                adoption = None
                if ephem.multiple > 1:
                    resolved_cand = None
                    try:
                        # Transit duration is period-invariant physically, but the
                        # alias-fold fit is biased. Scale by (P_phys/P_ali)^(1/3)
                        # (Kepler: T_dur ∝ P^(1/3) for a fixed host/geometry).
                        dur_scale = (ephem.physical_period / cand.period_best) ** (1.0 / 3.0)
                        resolved_cand = BLSDetector(
                            period_min=0.5, period_max=max(0.6, max_search),
                            snr_threshold=snr_threshold(), fap_threshold=fap_threshold(),
                            frequency_factor=float(T.threshold("detector_frequency_factor")),
                        ).run_at_period(
                            time=time_arr, flux=flux_arr,
                            target_period=ephem.physical_period,
                            duration_days=cand.transit_duration * dur_scale,
                        )
                    except Exception as e:
                        logger.warning(f"{tic_id}: resolution re-BLS failed: {e}")
                    if resolved_cand is not None:
                        _, rcand_res = validate_candidate(
                            tic_id, resolved_cand, time_arr, flux_arr, target.stellar, run_dir,
                            engine=engine,
                        )
                        # Adopt ONLY when the fundamental re-validates; otherwise
                        # keep the original verdict (zero regression).
                        if rcand_res["validation_status"] in CERTIFIED:
                            adoption = (resolved_cand, rcand_res)

                if adoption is not None:
                    resolved_cand, rcand_res = adoption
                    certified = resolved_cand
                    resolution_applied = True
                    res.update({
                        "period_found": resolved_cand.period_best,
                        "snr": resolved_cand.snr,
                        "cvs": rcand_res["cvs"],
                        "cvs_verdict": rcand_res["cvs_verdict"],
                        "sovereign_verdict_c99": rcand_res.get("sovereign_verdict_c99"),
                        "c99_error": rcand_res.get("c99_error"),
                        "validation_status": rcand_res["validation_status"],
                        "sovereign_verdict": rcand_res.get("sovereign_verdict"),
                        "ladder_rank": n_tested,
                        "ephemeris_resolution": ephem.classifier,
                        "ephemeris_physical_period": round(ephem.physical_period, 6),
                        "ephemeris_pattern": ephem.pattern,
                        "ephemeris_evidence": ephem.evidence,
                    })
                else:
                    res.update({
                        "period_found": cand.period_best,
                        "snr": cand.snr,
                        "cvs": cand_res["cvs"],
                        "cvs_verdict": cand_res["cvs_verdict"],
                        "sovereign_verdict_c99": cand_res.get("sovereign_verdict_c99"),
                        "c99_error": cand_res.get("c99_error"),
                        "validation_status": cand_res["validation_status"],
                        "sovereign_verdict": cand_res.get("sovereign_verdict"),
                        "ladder_rank": n_tested,
                        "ephemeris_resolution": ephem.classifier,
                        "ephemeris_physical_period": round(ephem.physical_period, 6),
                        "ephemeris_pattern": ephem.pattern,
                        "ephemeris_evidence": ephem.evidence,
                    })
                if target.label_period:
                    res["period_error_pct"] = abs(res["period_found"] - target.label_period) / target.label_period * 100.0
                break

        res["n_candidates_tested"] = n_tested
        if engine == "c99":
            # Nothing certified: report the first (highest-power) candidate's
            # sovereign outcome for honest FALSE_POSITIVE accounting.
            res["sovereign_verdict_c99"] = cand_res.get("sovereign_verdict_c99")
            res["c99_error"] = cand_res.get("c99_error")
        if certified is None:
            if target.label_period:
                res["period_error_pct"] = abs(bls.period_best - target.label_period) / target.label_period * 100.0
            res["validation_status"] = first_status or "UNCERTAIN"
            res["period_found"] = bls.period_best
            res["sovereign_verdict"] = None
            res["error"] = "No ladder candidate certified by sovereign validator"

    except Exception as e:
        res["error"] = str(e)
        logger.error(f"Error evaluating {tic_id}: {e}")

    finally:
        res["processing_time_sec"] = time.time() - start
    return res


# ─────────────────────────────────────────────────────────────────────────────
# True-target tiers
# ─────────────────────────────────────────────────────────────────────────────

TRUE_TIERS = [  # (weight, snr_lo, snr_hi)
    ("hard",  5.5,  8.0),
    ("medium", 8.0, 14.0),
    ("easy", 14.0, 30.0),
]

FALSE_PLAN = [
    ("eb",            8),
    ("grazing_eb",    8),
    ("rotation",     16),
    ("single_event", 16),
    ("noise",        32),
]


def build_true_set(n: int, seed: int = 20260814) -> List:
    import benchmarks_controlled.synthetic as S
    rng = np.random.default_rng(seed)
    targets = []
    # log-uniform periods 0.6..13.4
    lo, hi = math.log(0.6), math.log(13.4)
    i = 0
    while len(targets) < n:
        period = float(np.exp(rng.uniform(lo, hi)))
        tier = TRUE_TIERS[i % len(TRUE_TIERS)][1:]
        snr = float(rng.uniform(tier[0], tier[1]))
        targets.append(S.generate_true_planet(len(targets), period, snr, seed=seed))
        i += 1
    return targets


def build_false_set(n: int, seed: int = 20260814) -> List:
    import benchmarks_controlled.synthetic as S
    rng = np.random.default_rng(seed + 99)
    targets = []
    # build per plan proportional to n
    total = sum(sz for _, sz in FALSE_PLAN)
    counts = {k: max(1, int(round(n * sz / total))) for k, sz in FALSE_PLAN}
    idx = {k: 0 for k, _ in FALSE_PLAN}
    lo, hi = math.log(0.6), math.log(13.4)
    while len(targets) < n:
        for kind, _ in FALSE_PLAN:
            if idx[kind] >= counts[kind]:
                continue
            k = idx[kind]
            period = float(np.exp(rng.uniform(lo, hi)))
            if kind == "eb":
                t = S.generate_false_eb(k, period, depth=float(rng.uniform(0.005, 0.03)), seed=seed)
            elif kind == "grazing_eb":
                t = S.generate_false_grazing_eb(k, period, depth=float(rng.uniform(0.005, 0.02)), seed=seed)
            elif kind == "rotation":
                t = S.generate_false_rotation(k, period, amplitude=float(rng.uniform(0.001, 0.008)), seed=seed)
            elif kind == "single_event":
                t = S.generate_false_single_event(k, depth=float(rng.uniform(0.005, 0.05)), seed=seed)
            else:
                t = S.generate_false_noise(k, white_ppm=float(rng.uniform(150.0, 1500.0)), seed=seed)
            targets.append(t)
            idx[kind] += 1
            if len(targets) >= n:
                break
    return targets


# ─────────────────────────────────────────────────────────────────────────────
# Metrics
# ─────────────────────────────────────────────────────────────────────────────

def compute_metrics(true_results: List[Dict], false_results: List[Dict]) -> Dict:
    CERTIFIED = ("NEW_DISCOVERY", "OFFLINE_NEW_DISCOVERY")

    tp_period = 0        # true target: detected + correct period + certified
    tp_any = 0           # true target: detected + certified (any period)
    fn = 0
    fp_contamination = 0 # false set: wrongly certified as planet
    fp_wrong_ephem = 0   # true set: certified at a WRONG period
    tn = 0
    det_any_true = 0

    for r in true_results:
        if r["detected"]:
            det_any_true += 1
        certified = r["validation_status"] in CERTIFIED
        ephem_ok = r["period_error_pct"] is not None and r["period_error_pct"] <= 5.0
        if certified and ephem_ok:
            tp_period += 1
            tp_any += 1
        elif certified:
            tp_any += 1     # signal certified but at WRONG ephemeris
            fp_wrong_ephem += 1
        else:
            fn += 1

    for r in false_results:
        certified = r["validation_status"] in CERTIFIED
        if certified:
            fp_contamination += 1
        else:
            tn += 1

    fp = fp_contamination + fp_wrong_ephem
    n_true = len(true_results)
    n_false = len(false_results)
    recall_period = tp_period / n_true if n_true else 0.0
    recall_any = tp_any / n_true if n_true else 0.0
    fpr_contamination = fp_contamination / n_false if n_false else 0.0
    precision = tp_period / (tp_period + fp) if (tp_period + fp) else 0.0
    f1 = 2 * precision * recall_period / (precision + recall_period) if (precision + recall_period) else 0.0

    return {
        "n_true": n_true, "n_false": n_false,
        "tp_period": tp_period, "tp_any": tp_any,
        "fp_contamination": fp_contamination, "fp_wrong_ephem": fp_wrong_ephem,
        "fp": fp, "fn": fn, "tn": tn,
        "detected_any_true": det_any_true,
        "recall_period": round(recall_period, 4),
        "recall_any": round(recall_any, 4),
        "fpr": round(fpr_contamination, 4),
        "precision": round(precision, 4),
        "f1": round(f1, 4),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Report
# ─────────────────────────────────────────────────────────────────────────────

def write_report(run_dir: Path, true_res, false_res, metrics, started_at: str) -> None:
    report = f"""# Axiom-ZSpace Controlled-Honesty Evaluation

**Generated:** {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}
**Run:** {run_dir.name}
**Started:** {started_at}
**Method:** full local pipeline, BLIND search (no period prior, no ground-truth
expected period), validator archive query stubbed OFFLINE so the sovereign
verdict reflects ONLY intrinsic physics/QA gating.

## Key Metrics

| Metric | Value |
|--------|-------|
| True (injected) planets | {metrics['n_true']} |
| Contamination targets | {metrics['n_false']} |
| **Recall@correct period (TP)** | {metrics['recall_period']*100:.1f}%  ({metrics['tp_period']}/{metrics['n_true']}) |
| **Detection recall (any period)** | {metrics['recall_any']*100:.1f}%  ({metrics['tp_any']}/{metrics['n_true']}) |
| Detected-at-any-period (raw BLS) | {metrics['detected_any_true']}/{metrics['n_true']} |
| **FPR (contamination certified as planet)** | {metrics['fpr']*100:.2f}%  ({metrics['fp_contamination']} FP / {metrics['n_false']} ) |
| Wrong-ephemeris certs on TRUE set | {metrics['fp_wrong_ephem']} (certified at a period outside the 5% window) |
| **Precision (certified & correct-period)** | {metrics['precision']*100:.1f}% |
| **F1 (period-level)** | {metrics['f1']:.3f} |
| Confusion | TP={metrics['tp_period']} FP={metrics['fp']} FN={metrics['fn']} TN={metrics['tn']} |

## Recall by injected SNR

"""
    from collections import defaultdict
    buckets = defaultdict(lambda: [0, 0])
    for r in true_res:
        ts = r.get("target_snr")
        ok = (r["validation_status"] in ("NEW_DISCOVERY", "OFFLINE_NEW_DISCOVERY"))
        perr = r.get("period_error_pct")
        epoch_ok = perr is not None and perr <= 5.0
        if ts is None:
            continue
        b = "5.5-8" if ts < 8 else ("8-14" if ts < 14 else "14-30")
        buckets[b][1] += 1
        if ok and epoch_ok:
            buckets[b][0] += 1
    report += "| SNR tier | recovered/n | recall |\n|---|---|---|\n"
    for b in ("5.5-8", "8-14", "14-30"):
        rec, tot = buckets.get(b, (0, 0))
        report += f"| {b} | {rec}/{tot} | {rec/tot*100 if tot else 0:.1f}% |\n"
    report += "\n## Per-subkind contamination results\n\n| kind | count | certified-FP |\n|---|---|---|\n"
    from collections import Counter
    fb = defaultdict(list)
    for r in false_res:
        fb[r["subkind"]].append(r)
    for kind, arr in sorted(fb.items()):
        cert = sum(1 for r in arr if r["validation_status"] in ("NEW_DISCOVERY", "OFFLINE_NEW_DISCOVERY"))
        report += f"| {kind} | {len(arr)} | {cert} |\n"
    report += f"""

## Files

- `results_true.json` — per-target results for injected planets
- `results_false.json` — per-target results for contamination
- `validation/` — sovereign validation cards written by the pipeline
- `EVALUATION_REPORT.md` — this report

**Note on honesty:** raw "detected" flags exclude nothing; wrong-ephemeris
certifications on the true set are counted separately from contamination
false-certifications (the FPR row above is contamination-only).
"""
    (run_dir / "EVALUATION_REPORT.md").write_text(report, encoding="utf-8")


# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Controlled honesty evaluation")
    parser.add_argument("--true", type=int, default=80)
    parser.add_argument("--false", type=int, default=80)
    parser.add_argument("--out", default=None, help="run dir (optional)")
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--engine", choices=["python", "c99"], default="python",
                        help="sovereign engine: python reference or c99 binary (default: python)")
    args = parser.parse_args()

    run_dir = Path(args.out) if args.out else (
        Path(__file__).resolve().parent / "runs" /
        datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    )
    run_dir.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    logger.info(f"Controlled evaluation -> {run_dir}")

    true_targets = build_true_set(args.true, seed=args.seed)
    false_targets = build_false_set(args.false, seed=args.seed)
    logger.info(f"Built {len(true_targets)} true + {len(false_targets)} false targets")

    true_results, false_results = [], []
    for i, t in enumerate(true_targets, 1):
        logger.info(f"[true {i}/{len(true_targets)}] {t.tic_id} P={t.label_period:.3f} snr={t.meta.get('target_snr'):.1f}")
        true_results.append(evaluate_target(t, run_dir, engine=args.engine))
    for i, t in enumerate(false_targets, 1):
        logger.info(f"[false {i}/{len(false_targets)}] {t.tic_id} kind={t.subkind}")
        false_results.append(evaluate_target(t, run_dir, engine=args.engine))

    (run_dir / "results_true.json").write_text(
        json.dumps(true_results, indent=2), encoding="utf-8")
    (run_dir / "results_false.json").write_text(
        json.dumps(false_results, indent=2), encoding="utf-8")

    metrics = compute_metrics(true_results, false_results)
    write_report(run_dir, true_results, false_results, metrics, started)

    print("\n" + "=" * 60)
    print("CONTROLLED HONESTY EVALUATION")
    print("=" * 60)
    print(f"Recall@correct period : {metrics['recall_period']*100:.1f}%  ({metrics['tp_period']}/{metrics['n_true']})")
    print(f"Detection recall       : {metrics['recall_any']*100:.1f}%  ({metrics['tp_any']}/{metrics['n_true']})")
    print(f"FPR (contamination)    : {metrics['fpr']*100:.2f}%  ({metrics['fp_contamination']} FP / {metrics['n_false']})")
    print(f"Precision (period)     : {metrics['precision']*100:.1f}%")
    print(f"F1                     : {metrics['f1']:.3f}")
    print(f"Report: {run_dir / 'EVALUATION_REPORT.md'}")


if __name__ == "__main__":
    main()