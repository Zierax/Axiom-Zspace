"""
thresholds.py — Central tunable threshold catalog for the Axiom-ZSpace
detection/validation pipeline.

Every tunable number that affects detection, vetting, or certification lives
here. Values come from ``config/production.yaml`` (section ``thresholds``),
with per-profile sets so operators can choose their FPR/recall trade-off
without touching code:

    conservative — zero-false-positive priority (measured FPR 0/80 fixed-seed;
                    4.25% across the 800-target BIG400 suite)
    balanced     — default; same measured values as above
    sensitive    — recall priority, accepts some FPR (EXPERIMENTAL, requires
                   a full benchmark re-run before use)

Usage
-----
    from zspace_engine import thresholds as T
    T.threshold("fp7_density_band")      # current active-profile value
    T.active_profile()                   # e.g. "balanced"
    T.summary()                          # dump table
"""

from __future__ import annotations

from typing import Any, Dict, List

from zspace_engine.config import load_config

PROFILES: tuple = ("conservative", "balanced", "sensitive")
ACTIVE_KEY = "profile"


# ── Built-in defaults (used only when the YAML omits a key) ────────────────
# Values marked [M] are directly measured on the controlled benchmark
# (100 true + 80 false synthetic targets, seed 20260814).

PROFILE_DEFAULTS: Dict[str, Dict[str, Any]] = {
    # ── conservative: zero-FPR priority ──────────────────────────────────
    "conservative": {
        "fp1_snr_min": 5.5,                  # [M] BLS detection SNR floor
        "fp2_fap_max": 0.05,                 # [M] power-spectrum FAP ceiling
        "fp3_eo_sigma_max": 3.0,             # even/odd depth Δσ ceiling
        "fp4_shape_min": 0.4,                # [M] U/V shape ratio floor
        "fp5_secondary_snr_max": 3.0,        # phase-0.5 eclipse SNR ceiling
        "fp5b_secondary_ratio_max": 0.30,    # eclipse depth ratio ceiling
        "fp5c_alias_band": [0.20, 0.90],     # asymmetric alias-eclipse band
        "fp5c_alias_min_snr": 15.0,          # alias eclipse SNR floor
        "fp7_density_band": [0.2, 5.0],      # [M] transit/TIC density band
        "fp8_impact_max": 0.9,               # grazing impact-parameter cap
        "fp10_min_transits": 2,
        "verdict_max_fail_pass": 2,          # ≤2 fails → SOVEREIGN_PASS
        "verdict_max_fail_conditional": 3,   # ≤3 fails → CONDITIONAL_PASS
        "density_conflict_low": 0.5,
        "density_conflict_high": 2.0,
        "density_conflict_snr_min": 10.0,
        "ladder_k": 20,
        "ladder_min_relative_snr": 0.05,
        "detector_frequency_factor": 20.0,
        "detector_period_min": 0.5,
        "detector_period_max": 13.5,
        "coherent_override_enabled": False,  # [M] power-FAP firewall stays strict
        "coherent_min_snr": 6.5,
        "coherent_min_transits": 3,
        "coherent_min_dip_fraction": 0.6,
        "ev_odd_sigma_eb": 3.0,
        "ev_odd_pvalue_eb": 0.01,
        "ingress_fraction_vshape": 0.45,
        "secondary_ratio_veto": 0.30,
        "context_secondary_snr": 3.0,
        "context_density_mismatch": 0.20,
        "cvs_w_periodicity": 0.97,
        "cvs_w_depth": 0.83,
        "cvs_w_limb": 0.61,
        "cvs_w_secondary": 0.31,
        "no_detection_snr_floor": 5.5,
    },
    # ── balanced: default profile ────────────────────────────────────────
    "balanced": {
        "fp1_snr_min": 5.5,
        "fp2_fap_max": 0.05,
        "fp3_eo_sigma_max": 3.0,
        "fp4_shape_min": 0.4,
        "fp5_secondary_snr_max": 3.0,
        "fp5b_secondary_ratio_max": 0.30,
        "fp5c_alias_band": [0.20, 0.90],
        "fp5c_alias_min_snr": 15.0,
        "fp7_density_band": [0.2, 5.0],
        "fp8_impact_max": 0.9,
        "fp10_min_transits": 2,
        "verdict_max_fail_pass": 2,
        "verdict_max_fail_conditional": 3,
        "density_conflict_low": 0.5,
        "density_conflict_high": 2.0,
        "density_conflict_snr_min": 10.0,
        "ladder_k": 20,
        "ladder_min_relative_snr": 0.05,
        "detector_frequency_factor": 20.0,
        "detector_period_min": 0.5,
        "detector_period_max": 13.5,
        "coherent_override_enabled": False,
        "coherent_min_snr": 6.5,
        "coherent_min_transits": 3,
        "coherent_min_dip_fraction": 0.6,
        "ev_odd_sigma_eb": 3.0,
        "ev_odd_pvalue_eb": 0.01,
        "ingress_fraction_vshape": 0.45,
        "secondary_ratio_veto": 0.30,
        "context_secondary_snr": 3.0,
        "context_density_mismatch": 0.20,
        "cvs_w_periodicity": 0.97,
        "cvs_w_depth": 0.83,
        "cvs_w_limb": 0.61,
        "cvs_w_secondary": 0.31,
        "no_detection_snr_floor": 5.5,
    },
    # ── sensitive: recall priority — EXPERIMENTAL, NOT MEASURED ──────────
    # Every deviation from balanced requires a full benchmark re-run
    # (100 true + 80 false) before it may be adopted.
    "sensitive": {
        "fp1_snr_min": 5.0,
        "fp2_fap_max": 0.10,
        "fp3_eo_sigma_max": 3.0,
        "fp4_shape_min": 0.3,
        "fp5_secondary_snr_max": 3.0,
        "fp5b_secondary_ratio_max": 0.35,
        "fp5c_alias_band": [0.20, 0.90],
        "fp5c_alias_min_snr": 15.0,
        "fp7_density_band": [0.15, 6.0],
        "fp8_impact_max": 0.95,
        "fp10_min_transits": 2,
        "verdict_max_fail_pass": 3,
        "verdict_max_fail_conditional": 4,
        "density_conflict_low": 0.4,
        "density_conflict_high": 2.5,
        "density_conflict_snr_min": 10.0,
        "ladder_k": 20,
        "ladder_min_relative_snr": 0.05,
        "detector_frequency_factor": 20.0,
        "detector_period_min": 0.5,
        "detector_period_max": 13.5,
        "coherent_override_enabled": True,
        "coherent_min_snr": 6.5,
        "coherent_min_transits": 3,
        "coherent_min_dip_fraction": 0.55,
        "ev_odd_sigma_eb": 3.0,
        "ev_odd_pvalue_eb": 0.01,
        "ingress_fraction_vshape": 0.45,
        "secondary_ratio_veto": 0.30,
        "context_secondary_snr": 3.0,
        "context_density_mismatch": 0.25,
        "cvs_w_periodicity": 0.97,
        "cvs_w_depth": 0.83,
        "cvs_w_limb": 0.61,
        "cvs_w_secondary": 0.31,
        "no_detection_snr_floor": 5.0,
    },
}


# ── Catalog: human-readable metadata for every threshold ──────────────────
# direction: "gt" (value must be ABOVE), "lt" (BELOW), "range" (inside),
#            "flag" (boolean switch)
CATALOG: Dict[str, Dict[str, Any]] = {
    "fp1_snr_min": {
        "name": "FP-1 BLS SNR minimum",
        "unit": "σ", "direction": "gt", "weight": "critical",
        "purpose": "BLS signal-to-noise ratio floor: below this the periodic power peak is not statistically distinguishable from the noise.",
        "evidence": "Measured: truths with SNR ≥ 5.5 at the correct period were all detectable; false noise-only candidates with SNR 5.0-5.5 exist but were vetoed by other gates.",
        "pros_tight": "Fewer noise peaks enter the ladder (lower FPR).",
        "cons_tight": "Loses the weakest genuine transits (recall drop); several label-SNR≈5.5 truths sit right at the edge.",
        "pros_loose": "Recovers shallow planet transits.",
        "cons_loose": "Pure-noise folds measured at SNR up to 8 → FPR risk unless later gates are strong (see coherent_override notes).",
        "fpr_risk_loose": "MEDIUM",
    },
    "fp2_fap_max": {
        "name": "FP-2 power-spectrum FAP ceiling",
        "unit": "probability", "direction": "lt", "weight": "critical",
        "purpose": "False-Alarm-Probability of the BLS power peak under a self-calibrating exponential-tail noise model (red-noise conservative). This is the FPR firewall.",
        "evidence": "Measured: saturates at 0.5-1.0 for ALL candidates in red-noise-dominated spectra (truths AND noise) — it cannot separate; certification therefore depends on the coherent override or the raw SNR. 8 high-SNR truths with FAP≈1 died on this gate alone.",
        "pros_tight": "Strict firewall; nothing certifies on a weak null.",
        "cons_tight": "Frequent false rejection of genuine shallow-but-coherent signals (the reason the override exists).",
        "pros_loose": "Accepts more signals without touching FPR if the coherent path stays off (override OFF measured 0/80 fixed-seed; 4.25% across the 800-target BIG400 suite).",
        "cons_loose": "With override ON, pure-noise folds pass FAP 0.3-1.0 and certify → measured contamination FPR 12.5-62.5% (PROBE_FPR68 series).",
        "fpr_risk_loose": "LOW with override OFF, HIGH with override ON",
    },
    "fp3_eo_sigma_max": {
        "name": "FP-3 Even/Odd depth Δσ ceiling",
        "unit": "σ", "direction": "lt", "weight": "critical",
        "purpose": "Even and odd transit depths must be statistically identical; a large Δσ indicates an EB whose primary/secondary alternate.",
        "evidence": "Measured: truths sit at Δσ 0.05-2.0; grazing EBs in the false set reach Δσ 4-30 with SNR>100 → clean separation.",
        "pros_tight": "Strongest single EB killer (catches deep EBs even when secondary windows miss).",
        "cons_tight": "None measured within the observed distributions; do not tighten below ~2.5 or noise can trigger.",
        "pros_loose": "No measured benefit.",
        "cons_loose": "Deep EBs with Δσ 4-8 would certify → immediate FPR.",
        "fpr_risk_loose": "HIGH",
    },
    "fp4_shape_min": {
        "name": "FP-4 Shape Ratio (U vs V) floor",
        "unit": "ratio", "direction": "gt", "weight": "major",
        "purpose": "Ratio of wing residual power vs centre: U-shaped (flat-bottom) planetary transits score >1; V-shaped (EB/grazing) score <1.",
        "evidence": "Measured: correct-period truths with per-cadence noise score 0.48-0.99 (!!) — the old >1.0 threshold rejected 4/8 recovered truths on pure measurement noise. False set: only 2/80 top candidates score <0.4, both independently vetoed by critical gates.",
        "pros_tight": "Rejects V-shaped grazing EBs.",
        "cons_tight": "Rejects grazing planets and noise-realized shallow transits (measured 0.48-0.99).",
        "pros_loose": "Recovers grazing/near-grazing and shallow signals.",
        "cons_loose": "Admits some V-shaped false dips; safe ONLY while density (FP-7) or even/odd (FP-3) critical gates still separate them (measured: yes).",
        "fpr_risk_loose": "LOW at 0.4 (measured 0/80 fixed-seed; 4.25% BIG400), grows below ~0.3",
    },
    "fp5_secondary_snr_max": {
        "name": "FP-5 secondary-eclipse SNR ceiling",
        "unit": "σ", "direction": "lt", "weight": "critical",
        "purpose": "No significant flux dip at phase 0.5 (a real companion would show one).",
        "evidence": "Measured: EBs show secondary SNR 5-90; planets ≤2.2.",
        "pros_tight": "Kills EBs whose secondary is measured.",
        "cons_tight": "At very low primary depth the phase-0.5 window catches noise; sub-2.5 thresholds can fire randomly on faint folds.",
        "pros_loose": "Tolerates marginal secondary windows.",
        "cons_loose": "Direct EB admission risk.",
        "fpr_risk_loose": "HIGH",
    },
    "fp5b_secondary_ratio_max": {
        "name": "FP-5b secondary depth-ratio ceiling",
        "unit": "fraction", "direction": "lt", "weight": "critical",
        "purpose": "Phase-0.5 eclipse depth must be <30% of the primary (companion is not self-luminous).",
        "evidence": "Measured: planets ≤0.12; EBs 0.4-1.0. One recovered truth (SYN100076) sat at 0.33 — the 0.30 edge rejected it on noise.",
        "pros_tight": "Clean EB discriminator.",
        "cons_tight": "Borderline true values (0.30-0.35) get rejected on noise.",
        "pros_loose": "Recovers borderline planets.",
        "cons_loose": "Grazing EBs with compact secondaries (ratio 0.2-0.4) admitted; must stay below ~0.4.",
        "fpr_risk_loose": "MEDIUM",
    },
    "fp5c_alias_band": {
        "name": "FP-5c alias-eclipse ratio band",
        "unit": "fraction", "direction": "range", "weight": "critical",
        "purpose": "Fold at 2×/3× the candidate period: an EB alias shows an ASYMMETRIC 0.2-0.9-ratio eclipse at SNR≥15; a planet folds to symmetric transits (ratio≈1.0).",
        "evidence": "Measured: catches grazing-EB sub-harmonics that faked clean planet folds (false-set rows 0.75/0.94/1.16 density in-band are rejected by FP-5c).",
        "pros_tight": "The only gate that sees hidden EB secondaries folded onto phase 0.",
        "cons_tight": "Requires SNR≥15 to avoid noise flags; at huge EB SNRs (100+) it fires reliably.",
        "pros_loose": "None.",
        "cons_loose": "EB aliases certify as clean planets — measured FPR source.",
        "fpr_risk_loose": "HIGH",
    },
    "fp7_density_band": {
        "name": "FP-7 stellar-density ratio band",
        "unit": "ratio", "direction": "range", "weight": "critical",
        "purpose": "Stellar density implied by transit depth+duration must match the TIC-stellar density (ratio≈1); an EB's half-duration geometry violates it.",
        "evidence": "Measured: recovered truths sit at 0.38-4.0 (noise around 1.0); false-set top candidates sit at ≤0.2 (n=73/80) or ≥9.4 (grazing EBs). In [0.2,0.4)∪(2.5,5] there are ZERO false top candidates and 6 truths → the widened band is FPR-free.",
        "pros_tight": "Rejects density-divergent EBs.",
        "cons_tight": "Rejects grazing and noise-realized planets (0.4-2.5 produced 5/8 measured false rejections).",
        "pros_loose": "Recovers grazing/near-grazing planets.",
        "cons_loose": "Beyond [0.15, 6.0] grazing-EB densities (9-590) enter → FPR grows; keep 5.0 hard cap.",
        "fpr_risk_loose": "LOW to 5.0 (measured 0/80 fixed-seed; 4.25% BIG400), HIGH beyond ~6",
    },
    "fp8_impact_max": {
        "name": "FP-8 impact-parameter cap",
        "unit": "b", "direction": "lt", "weight": "moderate",
        "purpose": "Transit must not be grazing (b<0.9); V-shaped grazing geometry looks like an EB.",
        "evidence": "Measured: one recovered truth (SYN100020) at b=0.95 — as a moderate-weight fail it does NOT trip the circuit breaker, so the cap can stay 0.9.",
        "pros_tight": "Flares the V-grazing alarm.",
        "cons_tight": "Rejects the most grazing planets (b 0.9-1.0).",
        "pros_loose": "Admits grazing planets.",
        "cons_loose": "Grazing EBs also live at b→1; loose caps need the density cap ≤6 to stay safe.",
        "fpr_risk_loose": "MEDIUM (coupled with fp7)",
    },
    "fp10_min_transits": {
        "name": "FP-10 minimum independent transits",
        "unit": "count", "direction": "gt", "weight": "critical",
        "purpose": "At least 2 transits must be observed; single-event dips are never certifiable.",
        "evidence": "Measured: single-event false targets all fail here; every true target (P≤13.5 d over ~1050 d) has ≫2.",
        "pros_tight": "Kills single-event contamination.",
        "cons_tight": "None in the synthetic band (P≤13.5).",
        "pros_loose": "None.",
        "cons_loose": "Single transients certify — unacceptable.",
        "fpr_risk_loose": "ABSOLUTE (design invariant)",
    },
    "verdict_max_fail_pass": {
        "name": "Verdict: max fails for SOVEREIGN_PASS",
        "unit": "count", "direction": "lt", "weight": "meta",
        "purpose": "All critical gates must pass; up to N non-critical fails tolerated.",
        "evidence": "Measured: recovered truths have 0-1 fails after the gate widening.",
        "pros_tight": "Fewer marginal candidates certified.",
        "cons_tight": "Moderate-weight noisy measures (FP-4/FP-8) can veto otherwise clean signals.",
        "pros_loose": "Tolerates noisy minor gates.",
        "cons_loose": "More noise candidates pass; only with strong critical gates.",
        "fpr_risk_loose": "LOW (critical breaker dominates)",
    },
    "coherent_override_enabled": {
        "name": "FP-2 coherent multi-transit override",
        "unit": "flag", "direction": "flag", "weight": "meta",
        "purpose": "Let overwhelming repeated-observation evidence (≥3 consistent transits + SNR + no secondary) override the saturated power-FAP.",
        "evidence": "Measured: ON caused FPR explosion (contamination FPR up to 62.5% on 8-target probes — PROBE_FPR68 series); OFF measured 0/80 fixed-seed and 4.25% across the 800-target BIG400 suite. The dip-fraction discriminator at 0.6 is a coin-flip for real shallow transits (0.5-0.55 measured) — reliable only above per-transit SNR≈3.",
        "pros_tight(OFF)": "FPR stays 0/80 on the fixed-seed suite; 4.25% across BIG400.",
        "cons_tight(OFF)": "FAP-saturated truths stay rejected unless other gates pass (SYN100037).",
        "pros_loose(ON)": "Rescues FAP-saturated coherent signals.",
        "cons_loose(ON)": "Measured FPR explosion at dip-fraction 0.6; needs statistical (binomial) per-epoch test before re-enabling.",
        "fpr_risk_loose": "HIGH (measured)",
    },
    "coherent_min_snr": {
        "name": "Override minimum matched SNR", "unit": "σ", "direction": "gt", "weight": "meta",
        "purpose": "Floor for the coherent path. Measured noise folds reach 6-8σ → 6.5 is the expected value of the noise, not a separation point.",
        "evidence": "Measured: noise folds at SNR 6-8 (PROBE_FPR68).",
        "pros_tight": "Harder to exploit by noise.",
        "cons_tight": "Weak-but-real signals can't override.",
        "pros_loose": "More overrides.",
        "cons_loose": "FPR with override ON and SNR floor ≤6.5 is already measured-explosive.",
        "fpr_risk_loose": "HIGH (with override ON)",
    },
    "coherent_min_dip_fraction": {
        "name": "Override dip-recurrence floor",
        "unit": "fraction", "direction": "gt", "weight": "meta",
        "purpose": "Fraction of epochs showing a dip; genuine periodicity ~1.0, ephemeral noise ~0.5.",
        "evidence": "Measured: noise folds show dipF≈0.9+ while passing FAP (probes certify 12.5-62.5% of contamination) — the metric is useless at SNR<3 per transit (sign flips); truths at per-transit SNR 1-2 measure 0.5-0.55, below the 0.6 floor.",
        "pros_tight": "Feels safe at 0.6+.",
        "cons_tight": "Rejects exactly the shallow truths it was meant to rescue.",
        "pros_loose": "Nothing without a binomial significance test.",
        "cons_loose": "Noise folds show dipF≈0.9+ → they certify → measured FPR explosion. Keep the floor, keep override OFF.",
        "fpr_risk_loose": "HIGH, do not loosen while override is ON",
    },
    "ladder_k": {
        "name": "Ladder depth (candidates tested)",
        "unit": "count", "direction": "gt", "weight": "meta",
        "purpose": "How many period candidates are validated per target; recall grows with k, CPU scales with k.",
        "evidence": "Measured: correct-period candidate is the #1 ladder entry for 7/8 recovered truths; sub-harmonics occupy ranks 1-3 on weaker targets.",
        "pros_tight": "Faster.",
        "cons_tight": "Misses weaker true periods ranked lower.",
        "pros_loose": "Better weak-target recall (alias + true candidates both tested).",
        "cons_loose": "CPU cost; FPR unchanged (each candidate still must pass gates).",
        "fpr_risk_loose": "NONE",
    },
    "detector_frequency_factor": {
        "name": "BLS frequency oversampling", "unit": "x" , "direction": "gt", "weight": "meta",
        "purpose": "Frequencies per trial per bin; 20× recovers shallow peaks that 1× misses (prior 30k-bin runs destroyed weak signals).",
        "evidence": "Measured: weak truths (SNR≈6-8) only appeared at factor 20 on unbinned data.",
        "pros_tight": "Faster.",
        "cons_tight": "Aliasing/missed weak peaks.",
        "pros_loose": "Weak-signal recall.",
        "cons_loose": "CPU cost, more ladder noise entries (gated downstream).",
        "fpr_risk_loose": "NONE",
    },
    "ev_odd_sigma_eb": {
        "name": "Auditor Even/Odd EB Δσ flag", "unit": "σ", "direction": "gt", "weight": "auditor",
        "purpose": "EB flag inside the auditors (paired with p<0.01).",
        "evidence": "Aligned with FP-3 (3.0).",
        "pros_tight": "Flags EBs earlier.",
        "cons_tight": "False flags at noisy depths.",
        "pros_loose": "None.",
        "cons_loose": "EBs pass auditor stage.",
        "fpr_risk_loose": "MEDIUM",
    },
    "ingress_fraction_vshape": {
        "name": "Ingress-fraction V-shape flag", "unit": "fraction", "direction": "gt", "weight": "auditor",
        "purpose": "Ingress/duration >0.45 ⇒ V-shaped ⇒ EB-like flag in auditor scoring.",
        "evidence": "Complements FP-4; grazing planets can exceed 0.45 — penalty only, not hard veto.",
        "pros_tight": "Stronger V-shape separation.",
        "cons_tight": "Penalizes grazing planets in CVS.",
        "pros_loose": "Fairer to grazing planets.",
        "cons_loose": "EB-like shapes less penalized.",
        "fpr_risk_loose": "LOW (CVS is soft)",
    },
    "context_density_mismatch": {
        "name": "Context a/R★ mismatch tolerance", "unit": "fraction", "direction": "lt", "weight": "context",
        "purpose": "Catalog vs transit a/R★ fractional deviation above 20% flags density mismatch in context checks.",
        "evidence": "Soft flag; feeds CVS penalties not gates.",
        "pros_tight": "Flags suspicious densities.",
        "cons_tight": "Noise excursions flag clean planets.",
        "pros_loose": "None measured.",
        "cons_loose": "Density-divergent EBs less penalized.",
        "fpr_risk_loose": "LOW",
    },
    "cvs_w_periodicity": {
        "name": "CVS weight: periodicity (S_P)", "unit": "w", "direction": "meta", "weight": "meta",
        "purpose": "Composite Vitality Score mixing weight; the bench derives CVS = (0.97·S_P + 0.83·S_depth + 0.61·S_limb + 0.31·0.5)/Σw.",
        "evidence": "CVS ≥0.80 (core.THRESHOLD_PLANET, mirrored in config detection.cvs_planet_threshold) labels PLANET CANDIDATE; all certified candidates in the benchmark have CVS≥0.8.",
        "pros_tight": "Periodicity dominates classification.",
        "cons_tight": "Depth/shape quality underweighted.",
        "pros_loose": "More balanced score.",
        "cons_loose": "Score drifts with pipeline changes.",
        "fpr_risk_loose": "LOW (informational)",
    },
    "no_detection_snr_floor": {
        "name": "NO_DETECTION floor", "unit": "σ", "direction": "gt", "weight": "meta",
        "purpose": "If no ladder candidate exceeds this SNR the target is declared NO_DETECTION (undetected) instead of certified/FP.",
        "evidence": "Aligned with fp1_snr_min; governs reported recall, not FPR.",
        "pros_tight": "Honest accounting.",
        "cons_tight": "Marks weak truths undetected even when recoverable.",
        "pros_loose": "More candidates enter validation.",
        "cons_loose": "None beyond CPU.",
        "fpr_risk_loose": "NONE",
    },
}

_ACTIVE_CACHE: str = ACTIVE_KEY
_PROFILE_CACHE: Dict[str, Dict[str, Any]] = {}


def _merge_profile(profile: str, raw: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(PROFILE_DEFAULTS[profile])
    for k, v in (raw or {}).items():
        if v is not None:
            out[k] = v
    return out


def _load() -> None:
    global _ACTIVE_CACHE, _PROFILE_CACHE
    if _PROFILE_CACHE:
        return
    cfg = load_config()
    th = cfg.get("thresholds") or {}
    active = str(th.get(ACTIVE_KEY, "balanced")).lower()
    if active not in PROFILES:
        active = "balanced"
    _ACTIVE_CACHE = active
    for p in PROFILES:
        _PROFILE_CACHE[p] = _merge_profile(p, th.get(p))


def active_profile() -> str:
    _load()
    return _ACTIVE_CACHE


def set_profile(name: str) -> None:
    global _ACTIVE_CACHE
    _load()
    if name not in PROFILES:
        raise ValueError(f"unknown profile {name!r}; choose from {PROFILES}")
    _ACTIVE_CACHE = name


def threshold(key: str, profile: str | None = None) -> Any:
    _load()
    p = profile or _ACTIVE_CACHE
    return _PROFILE_CACHE[p][key]


def profile_values(profile: str | None = None) -> Dict[str, Any]:
    _load()
    return dict(_PROFILE_CACHE[profile or _ACTIVE_CACHE])


def catalog() -> Dict[str, Dict[str, Any]]:
    return CATALOG


def summary() -> str:
    _load()
    lines = [f"Active profile: {_ACTIVE_CACHE}"]
    lines.append(f"{'key':<32}{'value':<20}{'unit':<10}{'weight':<10}name")
    for k, meta in CATALOG.items():
        v = _PROFILE_CACHE[_ACTIVE_CACHE].get(k)
        lines.append(f"{k:<32}{str(v):<20}{meta.get('unit',''):<10}{meta.get('weight',''):<10}{meta['name']}")
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse
    import sys

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    parser = argparse.ArgumentParser(description="Threshold catalog CLI")
    parser.add_argument("--show", action="store_true", help="print active profile table")
    parser.add_argument("--profile", default=None, help="view a specific profile table")
    parser.add_argument("--set", default=None, help="persist an active profile to config YAML")
    args = parser.parse_args()

    if args.set:
        _load()
        if args.set not in PROFILES:
            raise SystemExit(f"unknown profile {args.set!r}; choose from {PROFILES}")
        import yaml
        from zspace_engine.config import config_path
        path = config_path()
        if path is None:
            raise SystemExit("no config file found (set AXIOM_CONFIG or create config/production.yaml)")
        with open(path, "r", encoding="utf-8") as fh:
            raw = yaml.safe_load(fh) or {}
        raw.setdefault("thresholds", {})["profile"] = args.set
        with open(path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(raw, fh, sort_keys=False, allow_unicode=True)
        print(f"active profile -> {args.set} in {path}")

    if args.profile:
        if args.profile not in PROFILES:
            raise SystemExit(f"unknown profile {args.profile!r}")
        _load()
        vals = _PROFILE_CACHE[args.profile]
        for k, meta in CATALOG.items():
            print(f"{k:<32}{str(vals.get(k)):<20}{meta['name']}")
    if args.show or (not args.set and not args.profile):
        print(summary())