"""
core.py  ·  The Deterministic Kernel
=====================================
Implements the Composite Vitality Score (CVS) as the sole arbiter of
candidate planetary status.  No probabilistic black-boxes; every score
is derived from measurable photometric and astrometric quantities.

Physics contract
----------------
CVS = Σ(w_x · S_x) / Σ(w_x)

Component weights (hard-coded per Truthimatics V3.0 spec):
  w_P = 0.97   Periodicity
  w_δ = 0.83   Depth Consistency
  w_τ = 0.61   Limb Shape (Mandel-Agol discriminator)
  w_S = 0.31   Stellar Context

Score bands:
  CVS ≥ 0.80  → PLANET CANDIDATE
  CVS ≥ 0.55  → LIKELY PLANET CANDIDATE
  CVS ≥ 0.35  → AMBIGUOUS / REQUIRES FOLLOW-UP
  CVS  < 0.35 → FALSE POSITIVE
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Tuple, List


# ─────────────────────────────────────────────────────────────────────────────
# Hard-coded weight constants  (Truthimatics V3.0 §4.1)
# ─────────────────────────────────────────────────────────────────────────────
W_PERIODICITY: float = 0.97
W_DEPTH: float       = 0.83
W_LIMB: float        = 0.61
W_STELLAR: float     = 0.31

WEIGHT_REGISTRY: Dict[str, float] = {
    "periodicity": W_PERIODICITY,
    "depth":       W_DEPTH,
    "limb":        W_LIMB,
    "stellar":     W_STELLAR,
}

# CVS decision thresholds
THRESHOLD_PLANET     = 0.80
THRESHOLD_LIKELY     = 0.55
THRESHOLD_AMBIGUOUS  = 0.35


# ─────────────────────────────────────────────────────────────────────────────
# Orbital Mechanics — pure Newtonian / Equilibrium Temperature
# ─────────────────────────────────────────────────────────────────────────────

# Import IAU 2015 physical constants
from .constants import (
    G_SI,
    M_SUN,
    R_SUN,
    L_SUN,
    SIGMA_SB,
    AU,
    R_EARTH_SOLAR,
)


def semi_major_axis_au(period_days: float, stellar_mass_solar: float) -> Tuple[float, str]:
    """
    Kepler's Third Law  →  a = (G·M·T²/4π²)^(1/3)

    Returns
    -------
    a_au  : semi-major axis in AU
    proof : human-readable logic string
    """
    if period_days <= 0:
        raise ValueError(f"Period must be > 0; got {period_days}")
    if stellar_mass_solar <= 0:
        raise ValueError(f"Stellar mass must be > 0; got {stellar_mass_solar}")

    T_sec = period_days * 86400.0          # days → seconds
    M_kg  = stellar_mass_solar * M_SUN     # solar masses → kg

    a_m = (G_SI * M_kg * T_sec**2 / (4.0 * math.pi**2)) ** (1.0 / 3.0)
    a_au_val = a_m / AU

    proof = (
        f"Kepler III (Newton/Kepler law, IAU 2015 constants): "
        f"a = (G·M·T²/4π²)^(1/3) | "
        f"T={period_days:.4f} d, M={stellar_mass_solar:.4f} M☉ "
        f"→ a={a_au_val:.5f} AU"
    )
    return a_au_val, proof


def equilibrium_temperature_k(
    stellar_teff: float,
    stellar_radius_solar: float,
    a_au: float,
    albedo: float = 0.30,
) -> Tuple[float, str]:
    """
    T_eq = T_eff · (R_★ / 2a)^(1/2) · (1 - A_B)^(1/4)

    Parameters
    ----------
    stellar_teff          : stellar effective temperature (K)
    stellar_radius_solar  : stellar radius in solar radii
    a_au                  : semi-major axis in AU
    albedo                : Bond albedo (default 0.30 per spec)
    """
    if a_au <= 0:
        raise ValueError("Semi-major axis must be > 0")

    R_m = stellar_radius_solar * R_SUN
    a_m = a_au * AU

    T_eq = stellar_teff * math.sqrt(R_m / (2.0 * a_m)) * (1.0 - albedo) ** 0.25

    flags = []
    if not math.isfinite(T_eq) or T_eq > 1.0e5:
        flags.append(
            f"FLAG: PHYSICALLY_UNLIKELY_EQ_T (T_eq={T_eq:.2f} K non-finite or >1e5 K)"
        )

    proof = (
        f"T_eq (equilibrium re-radiation law; IAU 2015 constants): "
        f"T_eff·√(R★/2a)·(1-A)^0.25 | "
        f"T_eff={stellar_teff:.1f} K, R★={stellar_radius_solar:.3f} R☉, "
        f"a={a_au:.5f} AU, A={albedo:.2f} "
        f"→ T_eq={T_eq:.2f} K"
    )
    if flags:
        proof += " | " + " | ".join(flags)
    return T_eq, proof


def planet_radius_earth(
    transit_depth: float,
    stellar_radius_solar: float,
    stellar_radius_err_solar: float = 0.0,
    transit_depth_err: float = 0.0,
) -> Tuple[float, str]:
    """
    From transit depth: δ = (R_p / R_★)²  →  R_p = R_★ · √δ

    NOTE: assumes a flat-square (geometric) transit depth with no limb
    darkening correction. For uniformly bright disks this is exact; with
    limb darkening the depth-vs-(R_p/R★) relation differs at the ~1% level,
    so this is reported as the geometric estimate.

    Returns planet radius in Earth radii.
    """
    if transit_depth < 0 or transit_depth > 1:
        raise ValueError(f"Transit depth must be in [0, 1]; got {transit_depth}")
    if stellar_radius_solar <= 0:
        raise ValueError(f"Stellar radius must be > 0; got {stellar_radius_solar}")

    rp_solar = stellar_radius_solar * math.sqrt(transit_depth)
    rp_earth = rp_solar / R_EARTH_SOLAR

    # Error propagation (independent errors, geometric depth):
    #   σ_Rp/Rp = sqrt( (σ_R★/R★)² + (0.5·σ_δ/δ)² )
    if transit_depth > 0:
        sigma_rp_rel = math.sqrt(
            (stellar_radius_err_solar / stellar_radius_solar) ** 2
            + (0.5 * transit_depth_err / transit_depth) ** 2
        )
        rp_earth_err = rp_earth * sigma_rp_rel
        err_str = f" ± {rp_earth_err:.3f} R⊕"
    else:
        rp_earth_err = float("nan")
        err_str = " ± NaN (δ=0)"

    proof = (
        f"R_p (geometric transit-depth law; IAU 2015 constants): "
        f"R★·√δ | δ={transit_depth:.6f}, R★={stellar_radius_solar:.3f} R☉ "
        f"→ R_p={rp_earth:.3f} R⊕{err_str}"
    )
    return rp_earth, proof


# ─────────────────────────────────────────────────────────────────────────────
# Hard Physical Filters (pre-CVS rejection gates)
# ─────────────────────────────────────────────────────────────────────────────

# Thresholds for hard physical filters
MAX_PLANET_RADIUS_EARTH    = 25.0    # Nothing larger than ~2.2 Jupiter radii is a planet
SECONDARY_SNR_EB_THRESHOLD = 15.0    # Secondary eclipse SNR above this -> definite EB hard reject
SECONDARY_SNR_PENALTY_THRESHOLD = 5.0  # Secondary eclipse SNR above this -> CVS penalty
DENSITY_DEVIATION_HARD_MAX = 5.0     # 500% density deviation -> impossible planet
EVEN_ODD_EB_HARD_THRESHOLD = 5.0     # Even/odd sigma above this -> definite EB
MAX_TRANSIT_DEPTH          = 0.05    # 5% depth -> implied radius too large for planet
MAX_DURATION_FRACTION      = 0.25    # Transit longer than 25% of period is unphysical


@dataclass
class HardFilterResult:
    """Result of hard physical filter checks."""
    passed:      bool
    rejection:   str           # empty string if passed
    flags:       List[str] = field(default_factory=list)
    cvs_penalty: float = 1.0  # multiplier applied to CVS (1.0 = no penalty)
    proof:       str   = ""


def apply_hard_filters(
    planet_radius_earth:   float,
    transit_depth:         float,
    secondary_snr:         float = 0.0,
    density_deviation:     float = 0.0,    # fractional (e.g. 0.7 = 70%)
    even_odd_sigma:        float = 0.0,
    is_v_shape:            bool  = False,
    transit_duration_hrs:  float = 0.0,
    period_days:           float = 1.0,
) -> HardFilterResult:
    """
    Apply hard physical constraint filters BEFORE CVS computation.

    These are non-negotiable physics checks. If a candidate fails any
    of these, it is immediately classified as FALSE_POSITIVE regardless
    of what the CVS score would have been.

    Returns HardFilterResult with passed=False and a rejection reason
    if any hard constraint is violated.
    """
    flags = []
    cvs_penalty = 1.0

    # Stage 1: Geometric Filter -- radius cap
    if planet_radius_earth > MAX_PLANET_RADIUS_EARTH:
        return HardFilterResult(
            passed=False,
            rejection=f"TOO_LARGE: R_p={planet_radius_earth:.1f} R_earth > {MAX_PLANET_RADIUS_EARTH} R_earth",
            flags=["HARD_REJECT_RADIUS"],
            proof=f"Planet radius {planet_radius_earth:.1f} R_earth exceeds physical maximum "
                  f"({MAX_PLANET_RADIUS_EARTH} R_earth). No planet can be this large.",
        )

    # Stage 1b: Depth sanity check
    if transit_depth > MAX_TRANSIT_DEPTH:
        return HardFilterResult(
            passed=False,
            rejection=f"TOO_DEEP: depth={transit_depth*100:.2f}% > {MAX_TRANSIT_DEPTH*100:.0f}%",
            flags=["HARD_REJECT_DEPTH"],
            proof=f"Transit depth {transit_depth*100:.2f}% exceeds {MAX_TRANSIT_DEPTH*100:.0f}%. "
                  f"Implies a companion far too large to be a planet.",
        )

    # Stage 1c: Duration sanity check
    if period_days > 0 and transit_duration_hrs > 0:
        duration_fraction = (transit_duration_hrs / 24.0) / period_days
        if duration_fraction > MAX_DURATION_FRACTION:
            return HardFilterResult(
                passed=False,
                rejection=f"DURATION_TOO_LONG: {duration_fraction*100:.1f}% of period",
                flags=["HARD_REJECT_DURATION"],
                proof=f"Transit duration is {duration_fraction*100:.1f}% of orbital period. "
                      f"Physically impossible for a planetary transit.",
            )

    # Stage 2: Secondary eclipse -> definite EB (only at extreme SNR)
    if secondary_snr > SECONDARY_SNR_EB_THRESHOLD:
        return HardFilterResult(
            passed=False,
            rejection=f"ECLIPSING_BINARY: secondary_snr={secondary_snr:.1f} > {SECONDARY_SNR_EB_THRESHOLD}",
            flags=["HARD_REJECT_SECONDARY_ECLIPSE"],
            proof=f"Secondary eclipse detected at SNR={secondary_snr:.1f} (threshold {SECONDARY_SNR_EB_THRESHOLD}). "
                  f"Planets reflect too little light for this SNR. This is an eclipsing binary.",
        )

    # Moderate secondary eclipse -> penalty (could be hot Jupiter thermal emission
    # or artifact from wrong BLS period folding other transits at phase 0.5)
    if secondary_snr > SECONDARY_SNR_PENALTY_THRESHOLD:
        penalty = max(0.3, 1.0 - (secondary_snr - SECONDARY_SNR_PENALTY_THRESHOLD) / 20.0)
        cvs_penalty *= penalty
        flags.append(f"SECONDARY_PENALTY: snr={secondary_snr:.1f} -> CVS x{penalty:.2f}")

    # Stage 2b: Even/odd -> definite EB
    if even_odd_sigma > EVEN_ODD_EB_HARD_THRESHOLD:
        return HardFilterResult(
            passed=False,
            rejection=f"ECLIPSING_BINARY: even_odd_sigma={even_odd_sigma:.1f} > {EVEN_ODD_EB_HARD_THRESHOLD}",
            flags=["HARD_REJECT_EVEN_ODD"],
            proof=f"Even/odd transit depth difference at {even_odd_sigma:.1f} sigma "
                  f"(threshold {EVEN_ODD_EB_HARD_THRESHOLD}). Definite eclipsing binary.",
        )

    # Stage 2c: Combined V-shape + large radius -> EB candidate
    if is_v_shape and planet_radius_earth > 15.0:
        return HardFilterResult(
            passed=False,
            rejection=f"EB_CANDIDATE: V-shape + R_p={planet_radius_earth:.1f} R_earth",
            flags=["HARD_REJECT_VSHAPE_LARGE"],
            proof=f"V-shaped transit with large implied radius ({planet_radius_earth:.1f} R_earth). "
                  f"Strong indicator of eclipsing binary or blended EB.",
        )

    # Stage 3: Soft penalties (don't reject, but penalize CVS)
    # NOTE: single-channel design — density mismatch and V-shape are scored
    # once, in the Stellar Context score (S_S, context.py) and the limb-shape
    # audit (S_τ, auditors.py) respectively. Only a PHYSICALLY IMPOSSIBLE
    # density deviation (>500%) hard-rejects here; milder deviations are
    # handled by the S_S component so no double penalty is applied.
    if density_deviation > DENSITY_DEVIATION_HARD_MAX:
        return HardFilterResult(
            passed=False,
            rejection=f"DENSITY_IMPOSSIBLE: deviation={density_deviation*100:.0f}% > {DENSITY_DEVIATION_HARD_MAX*100:.0f}%",
            flags=["HARD_REJECT_DENSITY_IMPOSSIBLE"],
            proof=f"a/R★ deviation {density_deviation*100:.0f}% exceeds the physically "
                  f"impossible threshold ({DENSITY_DEVIATION_HARD_MAX*100:.0f}%). No planet "
                  f"can produce this stellar-density mismatch.",
        )

    proof_parts = []
    if flags:
        proof_parts = [f"Hard filters applied: {'; '.join(flags)}"]
    else:
        proof_parts = ["All hard physical filters PASSED"]

    return HardFilterResult(
        passed=True,
        rejection="",
        flags=flags,
        cvs_penalty=cvs_penalty,
        proof="; ".join(proof_parts),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Score containers
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ComponentScore:
    """A single component score with its logic proof."""
    name:   str
    weight: float
    value:  float          # normalised 0 → 1
    proof:  str
    flags:  List[str] = field(default_factory=list)

    def weighted(self) -> float:
        return self.weight * self.value

    def validate(self) -> None:
        if not (0.0 <= self.value <= 1.0):
            raise ValueError(
                f"Score '{self.name}' value {self.value:.4f} is outside [0, 1]. "
                "Physical constraint violated."
            )
        if self.weight not in WEIGHT_REGISTRY.values():
            raise ValueError(
                f"Score '{self.name}' weight {self.weight} is not a registered Truthimatics weight."
            )


@dataclass
class OrbitalMechanics:
    """Derived orbital parameters with provenance strings."""
    period_days:           float = 0.0
    semi_major_axis_au:    float = 0.0
    equilibrium_temp_k:    float = 0.0
    planet_radius_earth:   float = 0.0
    transit_depth:         float = 0.0
    albedo:                float = 0.30
    proof_semi_major:      str   = ""
    proof_teq:             str   = ""
    proof_radius:          str   = ""


# ─────────────────────────────────────────────────────────────────────────────
# Composite Vitality Score
# ─────────────────────────────────────────────────────────────────────────────

class CompositeVitalityScore:
    """
    Computes and stores the Composite Vitality Score (CVS).

    Usage
    -----
    cvs = CompositeVitalityScore()
    cvs.register(ComponentScore(...))
    result = cvs.compute()
    """

    def __init__(self) -> None:
        self._components: List[ComponentScore] = []
        self._cvs: float | None = None
        self._verdict: str = "UNCOMPUTED"
        self._proof_chain: List[str] = []
        self._vetoes: List[str] = []

    def register(self, score: ComponentScore) -> None:
        score.validate()
        self._components.append(score)

    def apply_veto(self, reason: str) -> None:
        """Register a critical-FP veto. Any veto caps CVS at the ambiguous
        threshold so the candidate cannot be classified as a planet, no
        matter how strong the periodicity/depth scores are."""
        self._vetoes.append(reason)

    def apply_hard_filter(self, hard_filter: 'HardFilterResult') -> None:
        """Store a hard filter result for application during compute()."""
        self._hard_filter = hard_filter

    def compute(self) -> float:
        if not self._components:
            raise RuntimeError("No component scores registered; cannot compute CVS.")

        total_weight     = sum(c.weight  for c in self._components)
        weighted_sum     = sum(c.weighted() for c in self._components)
        raw_cvs          = weighted_sum / total_weight

        # Build logic proof chain
        self._proof_chain = []
        for c in self._components:
            self._proof_chain.append(
                f"  [{c.name}] w={c.weight} * S={c.value:.4f} = {c.weighted():.4f} | {c.proof}"
            )
        self._proof_chain.append(
            f"  CVS_raw = {weighted_sum:.4f} / {total_weight:.4f} = {raw_cvs:.4f}"
        )

        # Apply hard filter penalty if present
        hard_filter = getattr(self, '_hard_filter', None)
        if hard_filter and not hard_filter.passed:
            # Hard rejection -- override CVS to 0
            self._cvs = 0.0
            self._verdict = "FALSE POSITIVE"
            self._proof_chain.append(
                f"  HARD_REJECT: {hard_filter.rejection} | CVS forced to 0.0"
            )
            return self._cvs

        if hard_filter and hard_filter.cvs_penalty < 1.0:
            self._cvs = raw_cvs * hard_filter.cvs_penalty
            self._proof_chain.append(
                f"  HARD_FILTER_PENALTY: CVS = {raw_cvs:.4f} x {hard_filter.cvs_penalty:.2f} = {self._cvs:.4f}"
            )
            if hard_filter.flags:
                for flag in hard_filter.flags:
                    self._proof_chain.append(f"  FLAG: {flag}")
        else:
            self._cvs = raw_cvs

        # Critical-FP veto gate: any veto forces CVS below the ambiguous
        # threshold so the candidate is classified as FALSE POSITIVE no
        # matter how strong the periodicity/depth scores are.
        if self._vetoes:
            cvs_cap = THRESHOLD_AMBIGUOUS - 1.0e-6  # below 0.35 → FALSE POSITIVE
            if self._cvs > cvs_cap:
                self._cvs = cvs_cap
            for reason in self._vetoes:
                self._proof_chain.append(f"  VETO: {reason} → CVS capped at {cvs_cap:.4f}")

        self._verdict = self._classify(self._cvs)
        return self._cvs

    @staticmethod
    def _classify(cvs: float) -> str:
        if cvs >= THRESHOLD_PLANET:
            return "PLANET CANDIDATE"
        elif cvs >= THRESHOLD_LIKELY:
            return "LIKELY PLANET CANDIDATE"
        elif cvs >= THRESHOLD_AMBIGUOUS:
            return "AMBIGUOUS / REQUIRES FOLLOW-UP"
        else:
            return "FALSE POSITIVE"

    @property
    def score(self) -> float:
        if self._cvs is None:
            raise RuntimeError("CVS not yet computed. Call .compute() first.")
        return self._cvs

    @property
    def verdict(self) -> str:
        return self._verdict

    @property
    def proof_chain(self) -> List[str]:
        return self._proof_chain

    @property
    def components(self) -> List[ComponentScore]:
        return self._components

    def summary_dict(self) -> dict:
        return {
            "cvs":     round(self._cvs, 6) if self._cvs is not None else None,
            "verdict": self._verdict,
            "components": [
                {
                    "name":     c.name,
                    "weight":   c.weight,
                    "score":    round(c.value, 6),
                    "weighted": round(c.weighted(), 6),
                    "proof":    c.proof,
                    "flags":    c.flags,
                }
                for c in self._components
            ],
            "proof_chain": self._proof_chain,
        }


# ─────────────────────────────────────────────────────────────────────────────
# VitalityMatrix — Orchestrator shim used by the full pipeline
# ─────────────────────────────────────────────────────────────────────────────

class VitalityMatrix:
    """
    High-level orchestrator.  Accepts raw component score values
    (produced by detectors, auditors, context modules) and returns
    a fully computed CVS with orbital mechanics.

    This is the single point of truth for the pipeline's decision.
    """

    def __init__(
        self,
        tic_id: str,
        planet_order: int = 1,
    ) -> None:
        self.tic_id       = str(tic_id)
        self.planet_order = planet_order
        self.zspace_id    = f"ZS-T-{self.tic_id}-{self.planet_order:02d}"
        self.cvs_engine   = CompositeVitalityScore()
        self.orbital      = OrbitalMechanics()

    def ingest_scores(
        self,
        s_periodicity: float, proof_p: str,
        s_depth:       float, proof_d: str,
        s_limb:        float, proof_l: str,
        s_stellar:     float, proof_s: str,
        flags_p: List[str] | None = None,
        flags_d: List[str] | None = None,
        flags_l: List[str] | None = None,
        flags_s: List[str] | None = None,
    ) -> None:
        self.cvs_engine.register(ComponentScore("periodicity", W_PERIODICITY, s_periodicity, proof_p, flags_p or []))
        self.cvs_engine.register(ComponentScore("depth",       W_DEPTH,       s_depth,       proof_d, flags_d or []))
        self.cvs_engine.register(ComponentScore("limb",        W_LIMB,        s_limb,        proof_l, flags_l or []))
        self.cvs_engine.register(ComponentScore("stellar",     W_STELLAR,     s_stellar,     proof_s, flags_s or []))

    def compute_orbital_mechanics(
        self,
        period_days: float,
        transit_depth: float,
        stellar_mass_solar: float,
        stellar_teff: float,
        stellar_radius_solar: float,
    ) -> OrbitalMechanics:
        a_au, proof_a = semi_major_axis_au(period_days, stellar_mass_solar)
        t_eq, proof_t = equilibrium_temperature_k(stellar_teff, stellar_radius_solar, a_au)
        r_p, proof_r  = planet_radius_earth(transit_depth, stellar_radius_solar)

        self.orbital = OrbitalMechanics(
            period_days=period_days,
            semi_major_axis_au=a_au,
            equilibrium_temp_k=t_eq,
            planet_radius_earth=r_p,
            transit_depth=transit_depth,
            albedo=0.30,
            proof_semi_major=proof_a,
            proof_teq=proof_t,
            proof_radius=proof_r,
        )
        return self.orbital

    def finalize(self) -> dict:
        cvs = self.cvs_engine.compute()
        return {
            "zspace_id":        self.zspace_id,
            "tic_id":           self.tic_id,
            "planet_order":     self.planet_order,
            "cvs":              self.cvs_engine.summary_dict(),
            "orbital_mechanics": {
                "period_days":         round(self.orbital.period_days, 6),
                "semi_major_axis_au":  round(self.orbital.semi_major_axis_au, 6),
                "equilibrium_temp_k":  round(self.orbital.equilibrium_temp_k, 2),
                "planet_radius_earth": round(self.orbital.planet_radius_earth, 4),
                "transit_depth_ppm":   round(self.orbital.transit_depth * 1e6, 2),
                "albedo":              self.orbital.albedo,
                "proof_semi_major":    self.orbital.proof_semi_major,
                "proof_teq":           self.orbital.proof_teq,
                "proof_radius":        self.orbital.proof_radius,
            },
        }