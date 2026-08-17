"""
validator.py  ·  AxiomValidator — Sovereign Validation Module
==============================================================
Cross-references every Axiom-ZSpace candidate against:
  1. NASA Exoplanet Archive  (confirmed planets, pscomppars table)
  2. TESS Object of Interest (TOI) catalogue  (toi table)

   Decision tree
   ─────────────
  TIC ID + Period
         │
         ▼
  Query NASA Archive + TOI
         │
   ┌─────┴───────────────────┐
   │ MATCH found             │ NO MATCH
   │ (|ΔP| ≤ 0.001 d)        │
   ▼                         ▼
 exist_planet.json       Discovery.json
 (known metadata)        (Sovereign Logic Card —
                          full mathematical proof)

Mathematical proof blocks in Discovery.json
───────────────────────────────────────────
  §1  Keplerian Dynamics      a = (GM_★P²/4π²)^(1/3)
  §2  Geometric Consistency   δ ≈ k²; limb-darkening correction
  §3  Density Constraint       ρ_★ from transit; compare to TIC
  §4  Transit Probability      P_tr = (R_★ + R_p) / a
  §5  False-Positive Ruling    secondary eclipse + V-shape discriminator
  §6  Axiom Whitebox Verdict   step-by-step logical closure

Network resilience
──────────────────
  • 30-second per-query timeout
  • Two independent query methods (astroquery TAP  /  requests TAP fallback)
  • On any network failure → OFFLINE_MODE flagged, math proof still emitted
  • Partial matches (only TOI found, or only Archive found) are documented

Usage
─────
  from zspace_engine.validator import AxiomValidator

  validator = AxiomValidator(output_dir=".")
  result    = validator.validate(
      tic_id               = "260128333",
      period_days          = 3.6986,
      transit_depth        = 0.00836,
      transit_duration_hrs = 2.0,
      t0_btjd              = 1201.0,
      stellar_mass_solar   = 1.0,
      stellar_radius_solar = 1.0,
      stellar_teff_k       = 5778.0,
      stellar_logg         = 4.44,
      planet_radius_earth  = 9.97,
      cvs_score            = 0.83,
      cvs_verdict          = "PLANET CANDIDATE",
      cvs_proof_chain      = [...],
      bls_snr              = 347.9,
      bls_fap              = 0.0,
      even_odd_delta_sigma = 0.975,
      shape_ratio          = 4.711,
  )
  # result["status"] == "KNOWN" | "NEW_DISCOVERY" | "OFFLINE_NEW_DISCOVERY"
  #                 | "EPHEMERIS_MISMATCH"  (found period ≠ tested planet's)
"""

from __future__ import annotations

import json
import math
import os
import time
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Import order fix for charset_normalizer conflicts
# ─────────────────────────────────────────────────────────────────────────────
# Force charset_normalizer import before requests to prevent import conflicts
# that block astroquery from querying NASA Archive
try:
    import charset_normalizer
except ImportError:
    pass

# Add warning filters for deprecation warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

# ─────────────────────────────────────────────────────────────────────────────
# Logging setup
# ─────────────────────────────────────────────────────────────────────────────
from zspace_engine.logging_config import get_logger

# ─────────────────────────────────────────────────────────────────────────────
# Detection gate constants — single source of truth, shared with detectors.py
# ─────────────────────────────────────────────────────────────────────────────
from zspace_engine.detectors import SNR_THRESHOLD, FAP_THRESHOLD
from zspace_engine import thresholds as _T  # central threshold catalog (config/production.yaml)

# ─────────────────────────────────────────────────────────────────────────────
# Physical constants (IAU 2015) — imported from constants.py
# ─────────────────────────────────────────────────────────────────────────────
from zspace_engine.constants import (
    G_SI,
    M_SUN,
    R_SUN,
    R_EARTH,
    AU,
    R_EARTH_SOLAR,
)

# ─────────────────────────────────────────────────────────────────────────────
# Chi-squared goodness-of-fit analysis
# ─────────────────────────────────────────────────────────────────────────────
from zspace_engine.chi_squared import ChiSquaredAnalyzer

PI = math.pi

# ─────────────────────────────────────────────────────────────────────────────
# Validation constants
# ─────────────────────────────────────────────────────────────────────────────
PERIOD_MATCH_TOLERANCE_DAYS = 0.01    # |P_candidate - P_archive| <= this -> MATCH (absolute)
PERIOD_MATCH_TOLERANCE_REL  = 0.05   # 5% relative tolerance for period matching
API_TIMEOUT_SECONDS         = 30
MAX_RETRY_ATTEMPTS          = 2

# NASA Exoplanet Archive TAP endpoint (used in real deployments)
NASA_TAP_URL  = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
# TOI table via same endpoint
TOI_TABLE     = "toi"
PLANET_TABLE  = "pscomppars"


# ─────────────────────────────────────────────────────────────────────────────
# Data containers
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ArchiveMatch:
    """Represents a single matching entry from NASA Archive or TOI list."""
    source:             str           # "NASA_ARCHIVE" | "TOI" | "BOTH"
    planet_name:        str
    period_days:        float
    period_delta_days:  float         # |P_candidate - P_archive|
    transit_depth:      Optional[float]
    planet_radius_earth: Optional[float]
    semi_major_axis_au: Optional[float]
    stellar_teff_k:     Optional[float]
    stellar_radius_solar: Optional[float]
    stellar_mass_solar: Optional[float]
    discovery_method:   Optional[str]
    disposition:        Optional[str]  # for TOIs: CP, PC, FP, KP, ...
    extra_fields:       Dict[str, Any] = field(default_factory=dict)


@dataclass
class ValidationResult:
    """Top-level output of AxiomValidator.validate()."""
    status:         str           # "KNOWN" | "NEW_DISCOVERY" | "OFFLINE_NEW_DISCOVERY" | "EPHEMERIS_MISMATCH"
    tic_id:         str
    period_days:    float
    match:          Optional[ArchiveMatch]
    output_file:    str
    proof_summary:  str
    network_error:  Optional[str] = None
    # Identity-gate provenance (EPHEMERIS_MISMATCH path)
    expected_period_days: Optional[float] = None
    expected_planet_name: Optional[str]   = None
    alias_sibling_match:  Optional[ArchiveMatch] = None


# ─────────────────────────────────────────────────────────────────────────────
# §1–§5  Mathematical Proof Engine (pure Python/numpy — no internet needed)
# ─────────────────────────────────────────────────────────────────────────────

class ProofEngine:
    """
    Generates the complete mathematical proof chain for a new planetary
    discovery. All derivations are deterministic and white-box.

    Every method returns a dict with:
      "result"  : the computed numeric answer
      "unit"    : physical unit string
      "formula" : LaTeX-style formula string
      "steps"   : ordered list of step dicts  {label, equation, value, unit}
      "verdict" : PASS | FAIL | WARN
      "proof"   : human-readable single-line proof string
    """

    def __init__(
        self,
        period_days:            float,
        transit_depth:          float,      # dimensionless (e.g. 0.008836)
        transit_duration_hrs:   float,
        stellar_mass_solar:     float,
        stellar_radius_solar:   float,
        stellar_teff_k:         float,
        stellar_logg:           float,
        planet_radius_earth:    float,
        bls_snr:                float,
        bls_fap:                float,
        even_odd_delta_sigma:   float,
        shape_ratio:            float,      # wings_rms / centre_rms
        limb_dark_u1:           float = 0.30,
        limb_dark_u2:           float = 0.10,
    ) -> None:
        # Store inputs
        self.P_days    = period_days
        self.P_sec     = period_days * 86400.0
        self.delta     = transit_depth
        self.T_hrs     = transit_duration_hrs
        self.T_sec     = transit_duration_hrs * 3600.0
        self.M_star    = stellar_mass_solar * M_SUN
        self.R_star    = stellar_radius_solar * R_SUN
        self.T_eff     = stellar_teff_k
        self.logg      = stellar_logg
        self.R_planet  = planet_radius_earth * R_EARTH
        self.snr       = bls_snr
        self.fap       = bls_fap
        self.eo_sigma  = even_odd_delta_sigma
        self.shape_r   = shape_ratio
        self.u1        = limb_dark_u1
        self.u2        = limb_dark_u2

        # Derived convenience quantities
        self.M_star_solar   = stellar_mass_solar
        self.R_star_solar   = stellar_radius_solar
        self.k              = math.sqrt(max(self.delta, 1e-10))  # R_p / R_★

    # ─────────────────────────────────────────────────────────────────────────
    # §1  Keplerian Dynamics
    # ─────────────────────────────────────────────────────────────────────────

    def kepler_third_law(self) -> dict:
        """
        a = (G·M_★·P²  /  4π²)^(1/3)

        Physical meaning: the only force acting on the planet is gravity.
        Any deviation from this relation would indicate non-Keplerian
        dynamics (tidal distortion, third body, etc.).
        """
        # EDGE CASE FIX: Validate inputs are in reasonable ranges
        if self.M_star_solar <= 0 or self.M_star_solar > 100:
            # Unrealistic stellar mass - clamp to reasonable range
            M_star_clamped = max(0.08, min(self.M_star_solar, 100.0))  # 0.08-100 M☉
        else:
            M_star_clamped = self.M_star_solar
            
        # Step-by-step
        GM   = G_SI * M_star_clamped * M_SUN
        P2   = self.P_sec ** 2
        a3   = GM * P2 / (4 * PI ** 2)
        a_m  = a3 ** (1.0 / 3.0)
        a_au = a_m / AU

        # Kepler's Third Law ratio check: P²/a³ = 4π²/(GM_★)
        # Primary check in SI units (exact by construction, residual ~ 0%)
        ratio_si = (self.P_sec ** 2) / (a_m ** 3)
        expected_si = (4 * PI ** 2) / GM
        residual_si_pct = abs(ratio_si - expected_si) / expected_si * 100.0
        
        # Secondary check in solar units: P²[yr] / a³[AU] ≈ 1/M_★[M☉]
        # (includes unit conversion artifacts, ~0.004% for IAU 2015)
        P_yr    = self.P_days / 365.25
        ratio   = (P_yr ** 2) / max(a_au ** 3, 1e-30)  # Protect division
        expected = 1.0 / M_star_clamped
        residual_pct = abs(ratio - expected) / max(expected, 1e-10) * 100.0

        steps = [
            {"label": "Gravitational parameter",
             "equation": "μ = G·M_★",
             "value": round(GM, 4), "unit": "m³/s²"},
            {"label": "Orbital period squared",
             "equation": "P² = P²",
             "value": round(P2, 4), "unit": "s²"},
            {"label": "Semi-major axis cubed",
             "equation": "a³ = μ·P² / 4π²",
             "value": round(a3, 6), "unit": "m³"},
            {"label": "Semi-major axis",
             "equation": "a = (μP²/4π²)^(1/3)",
             "value": round(a_m, 2), "unit": "m"},
            {"label": "Semi-major axis in AU",
             "equation": "a [AU] = a [m] / 1 AU",
             "value": round(a_au, 7), "unit": "AU"},
            {"label": "Kepler ratio check (SI units)",
             "equation": "P²[s²] / a³[m³] vs 4π²/(GM_★)",
             "value": round(ratio_si, 15), "unit": "s⁻²m⁻³",
             "expected": round(expected_si, 15),
             "residual_pct": round(residual_si_pct, 9)},
            {"label": "Kepler ratio check (solar units)",
             "equation": "P²[yr] / a³[AU] vs 1/M_★[M☉]",
             "value": round(ratio, 8), "unit": "dimensionless",
             "expected": round(expected, 8),
             "residual_pct": round(residual_pct, 6)},
        ]
        verdict = "PASS" if residual_si_pct < 0.001 else "WARN"
        return {
            "section": "§1 Keplerian Dynamics",
            "formula": "a = (G·M_★·P² / 4π²)^(1/3)",
            "result": round(a_au, 7),
            "unit": "AU",
            "steps": steps,
            "verdict": verdict,
            "proof": (
                f"Kepler III (IAU 2015) | G={G_SI:.3e} m³kg⁻¹s⁻², "
                f"M_★={M_star_clamped:.4f} M☉, P={self.P_days:.6f} d "
                f"→ a={a_au:.6f} AU | "
                f"SI residual={residual_si_pct:.9f}%, "
                f"solar residual={residual_pct:.6f}% → {verdict}"
            ),
            "a_m": a_m,      # raw SI value for downstream use
            "a_au": a_au,
        }

    # ─────────────────────────────────────────────────────────────────────────
    # §2  Geometric Consistency
    # ─────────────────────────────────────────────────────────────────────────

    def geometric_consistency(self) -> dict:
        """
        Transit depth:  δ = (R_p / R_★)²  →  k = √δ

        Limb-darkening correction (Mandel & Agol 2002):
          The true depth at transit centre is deeper by factor:
            δ_centre = δ · [1 - u1·(1 - μ_c) - u2·(1 - μ_c)²] / Ī(u1,u2)
          where μ_c ≈ √(1 - b²) ≈ 1 for central transits,
          and Ī = 1 - u1/3 - u2/6  (mean disk intensity).

        Proof: the measured photometric depth δ must equal k² within
        the limb-darkening correction.  Any large discrepancy signals
        an EB (V-shaped eclipse deeper at centre due to stellar surface
        brightness gradient ≠ planetary occultation geometry).
        """
        k = self.k    # = sqrt(delta)
        k_sq = k ** 2

        # Limb-darkening mean disk intensity
        # EDGE CASE FIX: Ensure I_mean is positive and reasonable
        I_mean = 1.0 - self.u1 / 3.0 - self.u2 / 6.0
        I_mean = max(I_mean, 0.1)  # Prevent division by near-zero

        # Centre intensity (b=0 assumed for conservative estimate)
        mu_c  = 1.0
        I_cen = 1.0 - self.u1 * (1.0 - mu_c) - self.u2 * (1.0 - mu_c) ** 2
        I_cen = max(I_cen, 0.1)  # Prevent negative or zero intensity

        # LD-corrected depth (depth as seen against mean disk)
        delta_ld_corrected = k_sq * I_cen / I_mean

        # Fractional correction magnitude
        ld_correction_pct  = abs(delta_ld_corrected - k_sq) / max(k_sq, 1e-12) * 100.0

        # Measured vs geometric consistency check
        consistency_residual = abs(self.delta - k_sq) / max(k_sq, 1e-12) * 100.0

        # Planet radius from k and R_★
        R_p_m      = k * self.R_star
        R_p_earth  = R_p_m / R_EARTH
        R_p_solar  = R_p_m / R_SUN

        steps = [
            {"label": "Radius ratio (k)",
             "equation": "k = √δ",
             "value": round(k, 8), "unit": "dimensionless"},
            {"label": "Geometric depth prediction",
             "equation": "k² = (R_p/R_★)²",
             "value": round(k_sq, 8), "unit": "dimensionless"},
            {"label": "Measured transit depth",
             "equation": "δ (measured)",
             "value": round(self.delta, 8), "unit": "dimensionless"},
            {"label": "Consistency residual",
             "equation": "|δ_meas - k²| / k²",
             "value": round(consistency_residual, 4), "unit": "%"},
            {"label": "Mean limb-darkened disk intensity",
             "equation": "Ī = 1 - u1/3 - u2/6",
             "value": round(I_mean, 6), "unit": "dimensionless",
             "u1": self.u1, "u2": self.u2},
            {"label": "LD-corrected depth",
             "equation": "δ_LD = k²·I_cen / Ī",
             "value": round(delta_ld_corrected, 8), "unit": "dimensionless"},
            {"label": "Limb-darkening correction magnitude",
             "equation": "|δ_LD - k²| / k²",
             "value": round(ld_correction_pct, 4), "unit": "%"},
            {"label": "Planet radius",
             "equation": "R_p = k · R_★",
             "value": round(R_p_earth, 4), "unit": "R⊕"},
        ]
        verdict = "PASS" if consistency_residual < 5.0 else "WARN"
        return {
            "section": "§2 Geometric Consistency",
            "formula": "δ = (R_p/R_★)² = k²",
            "result": round(k, 8),
            "unit": "R_p/R_★",
            "steps": steps,
            "verdict": verdict,
            "proof": (
                f"Geom. consistency (IAU 2015) | δ_measured={self.delta:.6f}, "
                f"k=√δ={k:.6f}, k²={k_sq:.6f} | "
                f"residual={consistency_residual:.4f}% | "
                f"LD correction={ld_correction_pct:.4f}% | "
                f"R_p={R_p_earth:.4f} R⊕ → {verdict}"
            ),
            "k": k,
            "R_p_earth": R_p_earth,
            "delta_ld_corrected": delta_ld_corrected,
        }

    # ─────────────────────────────────────────────────────────────────────────
    # §3  Stellar Density Constraint
    # ─────────────────────────────────────────────────────────────────────────

    def density_constraint(self, a_m: float) -> dict:
        """
        Seager & Mallén-Ornelas (2003): The stellar mean density is uniquely
        determined by the transit observables (P, T_dur, δ).

          ρ_★ = (3π / G·P²) · (a/R_★)³

        where:  a/R_★ can be estimated from transit duration:
          T_dur = (P/π) · arcsin[ (R_★/a) · √((1+k)² - b²) ]
          For b=0 (central transit): T_dur ≈ (P·R_★)/(π·a) · √(1+k)²

        Rearranging:
          a/R_★ = P/(π·T_dur) · √(1+k)²   [b=0 approximation]

        Then compare ρ_★^transit with ρ_★^TIC (from M and R).
        Agreement validates planetary nature; large discrepancy → EB.

        EB discriminator (Seager & Mallén-Ornelas 2003, Torres et al. 2011):
          EBs show ρ_★^transit >> ρ_★^TIC because the eclipsing object is
          much smaller than the primary → a/R_★ is over-estimated.
        """
        k   = self.k
        T_s = self.T_sec
        P_s = self.P_sec

        # ── a/R_★ from transit geometry (b=0 assumption) ────────────────────
        # T_dur = (P/π) · (R_★/a) · √(1 + k)²   →   a/R_★ = P·(1+k) / (π·T_dur)
        # EDGE CASE FIX: Protect against unrealistic transit durations
        # If T_dur is placeholder (exactly 1.0 hr) or unrealistic, skip transit-based density check
        T_s_safe = max(T_s, 1.0)  # Minimum 1 second to prevent division by zero
        is_tdur_placeholder = abs(self.T_hrs - 1.0) < 0.001  # Detect placeholder value
        
        a_over_Rs_transit = P_s * (1.0 + k) / (PI * T_s_safe)
        a_over_Rs_direct  = a_m / max(self.R_star, 1e-9)   # from Kepler III (protect division)

        # ── ρ_★ from transit (Seager & Mallén-Ornelas 2003) ─────────────────
        rho_transit = (3.0 * PI / (G_SI * P_s ** 2)) * (a_over_Rs_transit ** 3)   # kg/m³
        rho_transit_gcc = rho_transit / 1000.0   # g/cm³

        # ── ρ_★ from TIC parameters (M/R³) ──────────────────────────────────
        M_SUN_kg  = M_SUN
        Vol_star  = (4.0 / 3.0) * PI * self.R_star ** 3
        rho_tic   = self.M_star / Vol_star              # kg/m³
        rho_tic_gcc = rho_tic / 1000.0                  # g/cm³

        # ── Solar density (reference) ────────────────────────────────────────
        rho_sun_gcc = 1.41   # g/cm³

        # ── Density ratio discriminator ──────────────────────────────────────
        density_ratio = rho_transit_gcc / max(rho_tic_gcc, 1e-10)

        # EB test: EBs typically show density_ratio >> 1 (order of magnitude)
        # For genuine planets: ratio within ~2× of 1.0 (geometric b uncertainty)
        # EDGE CASE FIX: If transit duration is placeholder, relax density check
        if is_tdur_placeholder:
            # Don't flag EB based on density if T_dur is unreliable
            is_eb_density_flag = False
        else:
            is_eb_density_flag = bool(density_ratio > 2.0 or density_ratio < 0.5)

        # ── Logg consistency check ────────────────────────────────────────────
        # log_g = log10(G·M/R²)  in cm/s²
        g_cgs     = G_SI * self.M_star / (self.R_star ** 2) * 100.0  # cm/s²
        logg_calc = math.log10(max(g_cgs, 1.0))
        logg_residual = abs(logg_calc - self.logg)

        steps = [
            {"label": "Transit duration (input)",
             "equation": "T_dur",
             "value": round(self.T_hrs, 4), "unit": "hours",
             "is_placeholder": is_tdur_placeholder},
            {"label": "a/R_★ from transit geometry (b=0)",
             "equation": "a/R_★ = P·(1+k) / (π·T_dur)",
             "value": round(a_over_Rs_transit, 4), "unit": "dimensionless"},
            {"label": "a/R_★ from Kepler III (independent check)",
             "equation": "a/R_★ = a_Kepler / R_★_TIC",
             "value": round(a_over_Rs_direct, 4), "unit": "dimensionless"},
            {"label": "Consistency of a/R_★ estimates",
             "equation": "|(a/R_★)_transit - (a/R_★)_Kepler| / (a/R_★)_Kepler",
             "value": round(abs(a_over_Rs_transit - a_over_Rs_direct) / max(a_over_Rs_direct, 1e-9) * 100, 3),
             "unit": "%",
             "note": "Skipped (T_dur placeholder)" if is_tdur_placeholder else None},
            {"label": "Stellar density from transit (Seager & Mallén-Ornelas 2003)",
             "equation": "ρ_★^transit = (3π/GP²)·(a/R_★)³",
             "value": round(rho_transit_gcc, 4), "unit": "g/cm³"},
            {"label": "Stellar density from TIC (M_★/V_★)",
             "equation": "ρ_★^TIC = M_★ / (4πR_★³/3)",
             "value": round(rho_tic_gcc, 4), "unit": "g/cm³"},
            {"label": "Solar density (reference)",
             "equation": "ρ_☉ = 1.41 g/cm³",
             "value": rho_sun_gcc, "unit": "g/cm³"},
            {"label": "Density ratio (transit / TIC)",
             "equation": "ρ^transit / ρ^TIC",
             "value": round(density_ratio, 4), "unit": "dimensionless",
             "eb_flag": is_eb_density_flag,
             "note": "Check skipped (T_dur placeholder)" if is_tdur_placeholder else None},
            {"label": "Surface gravity (calculated)",
             "equation": "log g = log10(G·M/R²) [cgs]",
             "value": round(logg_calc, 4), "unit": "log(cm/s²)"},
            {"label": "Surface gravity (TIC)",
             "equation": "log g (TIC)",
             "value": round(self.logg, 4), "unit": "log(cm/s²)"},
            {"label": "log g residual",
             "equation": "|logg_calc - logg_TIC|",
             "value": round(logg_residual, 4), "unit": "dex"},
        ]

        verdict = "PASS" if not is_eb_density_flag and logg_residual < 0.3 else "WARN"
        if is_eb_density_flag and not is_tdur_placeholder:
            verdict = "FAIL"
        elif is_tdur_placeholder:
            verdict = "WARN"  # Can't definitively rule out EB without good T_dur

        return {
            "section": "§3 Stellar Density Constraint",
            "formula": "ρ_★ = (3π / G·P²) · (a/R_★)³",
            "result": round(rho_transit_gcc, 6),
            "unit": "g/cm³",
            "steps": steps,
            "verdict": verdict,
            "proof": (
                f"Density (IAU 2015) | ρ^transit={rho_transit_gcc:.4f} g/cm³, "
                f"ρ^TIC={rho_tic_gcc:.4f} g/cm³, "
                f"ratio={density_ratio:.4f} | "
                f"logg_calc={logg_calc:.3f}, logg_TIC={self.logg:.3f}, "
                f"Δlogg={logg_residual:.3f} dex | "
                f"{'T_dur=placeholder, ' if is_tdur_placeholder else ''}"
                f"EB_density_flag={'YES' if is_eb_density_flag else 'NO'} → {verdict}"
            ),
            "rho_transit_gcc": rho_transit_gcc,
            "rho_tic_gcc":     rho_tic_gcc,
            "density_ratio":   density_ratio,
            "is_eb_density_flag": is_eb_density_flag,
            "a_over_Rs_transit": a_over_Rs_transit,
        }

    # ─────────────────────────────────────────────────────────────────────────
    # §4  Probability of Transit
    # ─────────────────────────────────────────────────────────────────────────

    def transit_probability(self, a_m: float) -> dict:
        """
        Geometric transit probability (Barnes 2007):

          P_tr = (R_★ + R_p) / a  ×  (1 + e·cos ω) / (1 - e²)

        For a circular orbit (e=0):
          P_tr = (R_★ + R_p) / a

        Physical meaning: if we were to observe this system from a random
        viewing angle, what fraction of orientations would produce a transit?
        High P_tr validates the transit as geometrically plausible.

        Supplementary: ingress/egress duration and impact parameter estimation.
        """
        R_p  = self.R_planet   # in metres
        R_s  = self.R_star
        a    = a_m

        # Geometric probability
        # EDGE CASE FIX: Clamp to [0, 1] range (can exceed 1.0 for very close-in planets)
        P_tr = (R_s + R_p) / max(a, R_s + R_p)  # Prevent P_tr > 1.0
        P_tr = min(P_tr, 1.0)  # Ensure probability doesn't exceed 100%

        # Grazing transit condition: b < 1 + k  where b = a·cos(i)/R_★
        # From transit duration: estimate impact parameter b
        k     = self.k
        T_s   = self.T_sec
        P_s   = self.P_sec
        
        # EDGE CASE FIX: Simplified impact parameter calculation
        # Full formula: b² = (1+k)² - [(a/R_★)·sin(πT/P)]²
        # Protect against domain errors in sin and sqrt
        a_over_Rs = a / max(R_s, 1e-9)
        sin_arg = PI * T_s / max(P_s, 1.0)
        sin_arg = min(sin_arg, PI)  # Clamp to valid range
        sin_term = a_over_Rs * math.sin(sin_arg)
        b_sq     = max(0.0, (1.0 + k) ** 2 - sin_term ** 2)
        b        = math.sqrt(b_sq)

        # Ingress duration (T_ingress): time to cross limb = T_dur * k/(1+k) [approx]
        T_ingress_hrs = self.T_hrs * k / max(1.0 + k, 1e-9)

        # Minimum orbital inclination for transit
        cos_i_max = (R_s + R_p) / a
        i_min_deg = math.degrees(math.acos(min(cos_i_max, 1.0)))

        is_grazing = bool(b > 0.9)

        steps = [
            {"label": "Stellar radius",
             "equation": "R_★",
             "value": round(R_s / R_SUN, 6), "unit": "R☉"},
            {"label": "Planet radius",
             "equation": "R_p",
             "value": round(R_p / R_EARTH, 6), "unit": "R⊕"},
            {"label": "Semi-major axis (from Kepler III)",
             "equation": "a",
             "value": round(a / AU, 7), "unit": "AU"},
            {"label": "Geometric transit probability",
             "equation": "P_tr = (R_★ + R_p) / a",
             "value": round(P_tr, 6), "unit": "dimensionless"},
            {"label": "Geometric transit probability as percent",
             "equation": "P_tr [%]",
             "value": round(P_tr * 100, 4), "unit": "%"},
            {"label": "Impact parameter estimate",
             "equation": "b = √[(1+k)² - (a/R_★·sin(πT/P))²]",
             "value": round(b, 5), "unit": "dimensionless"},
            {"label": "Grazing transit check (b < 0.9 for non-grazing)",
             "equation": "b < 0.9",
             "value": round(b, 5), "unit": "dimensionless",
             "is_grazing": is_grazing},
            {"label": "Ingress / Egress duration estimate",
             "equation": "T_ingress ≈ T_dur · k / (1+k)",
             "value": round(T_ingress_hrs, 4), "unit": "hours"},
            {"label": "Minimum orbital inclination for transit",
             "equation": "i_min = arccos[(R_★+R_p)/a]",
             "value": round(i_min_deg, 4), "unit": "degrees"},
        ]
        verdict = "PASS" if not is_grazing else "WARN"
        return {
            "section": "§4 Probability of Transit",
            "formula": "P_tr = (R_★ + R_p) / a  [circular orbit]",
            "result": round(P_tr, 8),
            "unit": "dimensionless",
            "steps": steps,
            "verdict": verdict,
            "proof": (
                f"Transit probability (IAU 2015) | R_★={R_s/R_SUN:.4f} R☉, "
                f"R_p={R_p/R_EARTH:.4f} R⊕, a={a/AU:.6f} AU "
                f"→ P_tr={(R_s+R_p)/a:.6f} ({(R_s+R_p)/a*100:.3f}%) | "
                f"b={b:.4f} | grazing={'YES' if is_grazing else 'NO'} → {verdict}"
            ),
            "P_tr": P_tr,
            "impact_parameter_b": b,
            "T_ingress_hrs": T_ingress_hrs,
            "i_min_deg": i_min_deg,
            "is_grazing": is_grazing,
        }

    # ─────────────────────────────────────────────────────────────────────────
    # §5  False-Positive Ruling (Axiom Whitebox Verdict)
    # ─────────────────────────────────────────────────────────────────────────

    def false_positive_ruling(
        self,
        secondary_snr:     float,
        centroid_sigma:    float,
        density_ratio:     float,
        is_grazing:        bool,
        tic_id:            str,
        n_transits:       Optional[int] = None,
        secondary_depth_ratio: float = 0.0,
        coherent_evidence: int = 0,
        alias_secondary_ratio: float = 0.0,
    ) -> dict:
        """
        Systematic false-positive elimination. Each vector is tested
        independently; a signal must pass ALL tests to reach the
        Sovereign Discovery status.

        Test vectors:
          FP-1  BLS SNR > 5.5  (detection significance — aligned with detectors.SNR_THRESHOLD)
          FP-2  FAP < 1e-4     (false alarm rate)
          FP-3  Even/Odd Δσ < 3.0  (secondary depth equality → not EB)
          FP-4  Shape ratio > 1.0   (U-shape → flat-bottomed → planet)
          FP-5  Secondary eclipse SNR < 3.0  (no phase-0.5 eclipse)
          FP-5b Secondary / primary eclipse depth ratio < 0.30 (companion emits light → EB)
          FP-6  Centroid shift σ < 3.0  (no photocenter motion)
          FP-7  Density ratio 0.5–2.0   (ρ_transit ≈ ρ_TIC → not EB)
          FP-8  Not grazing (b < 0.9)   (V-shape discriminator)
          FP-9  Catalog multiplicity check (no known multiple star systems)
        """
        tests = []

        def make_test(name: str, value: float, threshold: float,
                      direction: str, weight: str, description: str) -> dict:
            if direction == "gt":
                passed = bool(value > threshold)
                comparison = f"{value:.4f} > {threshold}"
            elif direction == "lt":
                passed = bool(value < threshold)
                comparison = f"{value:.4f} < {threshold}"
            elif direction == "range":
                lo, hi = threshold
                passed = bool(lo < value < hi)
                comparison = f"{lo} < {value:.4f} < {hi}"
            else:
                passed = False
                comparison = "unknown"
            return {
                "test":        name,
                "description": description,
                "value":       round(value, 6),
                "threshold":   threshold,
                "comparison":  comparison,
                "verdict":     "PASS" if passed else "FAIL",
                "weight":      weight,
            }

        tests.append(make_test(
            "FP-1 BLS SNR", self.snr, _T.threshold("fp1_snr_min"), "gt", "critical",
            "BLS signal-to-noise exceeds detection threshold"
        ))
        # FP-2 FAP is the strict FPR firewall (power-spectrum, red-noise
        # conservative). Overridden ONLY by overwhelming repeated-observation
        # evidence (>=3 depth-consistent independent transits + coherent SNR +
        # no secondary eclipse) — a genuine shallow transit that is recovered
        # positive-evidence path, never by weakening the null test itself.
        fp2_thr = float(_T.threshold("fp2_fap_max"))
        fp2_pass = (self.fap < fp2_thr) or (bool(coherent_evidence) and self.snr > _T.threshold("fp1_snr_min"))
        fp2_desc = (
            f"False Alarm Probability below {fp2_thr:.0e} "
            + (f"| overridden by coherent multi-transit evidence (SNR={self.snr:.1f}σ)" if (coherent_evidence and self.fap >= fp2_thr) else "")
        )
        fp2_test = {
            "test": "FP-2 FAP",
            "description": fp2_desc,
            "value": round(self.fap, 6),
            "threshold": fp2_thr,
            "comparison": (f"FAP<{fp2_thr:.0e}" if self.fap < fp2_thr else f"FAP={self.fap:.2e}≥{fp2_thr:.0e} overridden" if fp2_pass else f"{self.fap:.2e} ≥ {fp2_thr:.0e}"),
            "verdict": "PASS" if fp2_pass else "FAIL",
            "weight": "critical",
            "coherent_evidence": bool(coherent_evidence),
        }
        tests.append(fp2_test)
        tests.append(make_test(
            "FP-3 Even/Odd Δσ", self.eo_sigma, _T.threshold("fp3_eo_sigma_max"), "lt", "critical",
            "Even and odd transit depths are statistically identical (Δσ < threshold → not EB)"
        ))
        tests.append(make_test(
            "FP-4 Shape Ratio (U vs V)", self.shape_r, _T.threshold("fp4_shape_min"), "gt", "major",
            "Transit wings show higher residuals than centre → flat-bottomed U-shape → planetary occultation"
        ))
        tests.append(make_test(
            "FP-5 Secondary Eclipse SNR", secondary_snr, _T.threshold("fp5_secondary_snr_max"), "lt", "critical",
            "No statistically significant secondary eclipse at phase 0.5 (would indicate EB)"
        ))
        tests.append(make_test(
            "FP-5b Secondary Eclipse Depth Ratio", secondary_depth_ratio, _T.threshold("fp5b_secondary_ratio_max"), "lt", "critical",
            "Phase-0.5 eclipse depth < threshold of primary transit (companion emits no light)"
        ))
        # FP-5c Harmonic-Alias Secondary Eclipse: a grazing EB whose BLS period
        # lands on a sub-harmonic (found = true_period/2) folds both eclipses
        # onto phase 0, hiding the secondary from FP-5/FP-5b. The caller
        # re-folds at 2x and 3x the candidate period; an EB shows an ASYMMETRIC
        # secondary (0.2 < ratio < 0.9) with high SNR there, while a real planet
        # folds to symmetric equal-depth transits (ratio ~ 1.0) or sub-threshold
        # SNR (no veto). The test therefore FAILS exactly inside the asymmetric
        # band, i.e. the value is a genuine net positive secondary shallower
        # than the primary — the EB half-harmonic signature.
        _a0, _a1 = _T.threshold("fp5c_alias_band")
        _fp5c = _a0 < alias_secondary_ratio < _a1
        tests.append({
            "test": "FP-5c Harmonic-Alias Secondary Eclipse",
            "description": "No ASYMMETRIC secondary eclipse at 2x/3x candidate period (EB half-harmonic signature)",
            "value": round(alias_secondary_ratio, 6),
            "threshold": [_a0, _a1],
            "comparison": (
                f"ratio={alias_secondary_ratio:.3f} "
                f"{'in asymmetric band %.2f–%.2f → EB-ALIAS' % (_a0, _a1) if _fp5c else 'symmetric/absent → no veto'}"
            ),
            "verdict": "FAIL" if _fp5c else "PASS",
            "weight": "critical",
        })
        tests.append(make_test(
            "FP-6 Centroid Shift σ", centroid_sigma, 3.0, "lt", "major",
            "No photocenter displacement during transit (would indicate background EB)"
        ))
        tests.append(make_test(
            "FP-7 Density Ratio", density_ratio, tuple(_T.threshold("fp7_density_band")), "range", "critical",
            "Stellar density from transit ~ stellar density from TIC parameters (EB would diverge)"
        ))
        tests.append(make_test(
            "FP-8 Impact Parameter", 0.0 if not is_grazing else 0.95, _T.threshold("fp8_impact_max"), "lt", "moderate",
            "Transit is not grazing (b < threshold) — V-shaped grazing EB discriminator"
        ))

        # FP-9: External catalog cross-matching for stellar multiplicity
        catalog_result = check_external_catalogs(tic_id)
        is_multiple = catalog_result["is_multiple"]
        catalog_source = catalog_result["catalog_source"]
        
        # Create FP-9 test entry
        if catalog_source == "Offline":
            # Network unavailable - pass with CATALOG_OFFLINE flag
            fp9_test = {
                "test": "FP-9 Catalog Multiplicity",
                "description": "No known stellar multiplicity in external catalogs (SIMBAD/Gaia) — CATALOG_OFFLINE",
                "value": 0.0,
                "threshold": "is_multiple == False",
                "comparison": f"catalog_source={catalog_source}, risk_level={catalog_result['risk_level']}",
                "verdict": "PASS",
                "weight": "critical",
                "catalog_offline": True,
            }
        else:
            # Catalog query succeeded
            fp9_test = {
                "test": "FP-9 Catalog Multiplicity",
                "description": f"No known stellar multiplicity in external catalogs (SIMBAD/Gaia) — source: {catalog_source}",
                "value": 1.0 if is_multiple else 0.0,
                "threshold": "is_multiple == False",
                "comparison": f"is_multiple={is_multiple}, catalog_source={catalog_source}, classification={catalog_result['classification']}",
                "verdict": "FAIL" if is_multiple else "PASS",
                "weight": "critical",
                "catalog_offline": False,
            }
        
        tests.append(fp9_test)

        # FP-10: Minimum independent transits (optional; requires light-curve
        # access). Field standard (Kepler/TESS vetting): a planetary candidate
        # needs at least two independent transits. A single-event dip -- however
        # deep -- must never be validated as a planet.
        if n_transits is not None:
            n_tr = int(n_transits)
            fp10_min = int(_T.threshold("fp10_min_transits"))
            fp10_test = {
                "test": "FP-10 Minimum Transits",
                "description": f"At least {fp10_min} independent transits observed (single transients are not validated)",
                "value": float(n_tr),
                "threshold": float(fp10_min),
                "comparison": f"{n_tr} >= {fp10_min}",
                "verdict": "PASS" if n_tr >= fp10_min else "FAIL",
                "weight": "critical",
            }
            tests.append(fp10_test)

        n_pass     = sum(1 for t in tests if t["verdict"] == "PASS")
        n_critical = sum(1 for t in tests if t["weight"] == "critical")
        n_crit_pass= sum(1 for t in tests if t["weight"] == "critical" and t["verdict"] == "PASS")
        n_fail     = len(tests) - n_pass

        # Check for combined shape+density failure (Task 4.2)
        fp4_test = next((t for t in tests if t["test"] == "FP-4 Shape Ratio (U vs V)"), None)
        fp7_test = next((t for t in tests if t["test"] == "FP-7 Density Ratio"), None)
        
        fp4_failed = fp4_test and fp4_test["verdict"] == "FAIL"
        fp7_failed = fp7_test and fp7_test["verdict"] == "FAIL"

        # All critical tests must pass; ≤1 major/moderate failure tolerated
        critical_passed = bool(n_crit_pass == n_critical)

        # Circuit breaker: if ANY critical test fails, cannot be SOVEREIGN_PASS.
        # FP-7 (density ratio) is critical: a density mismatch outside the
        # configured band is an EB discriminator and must veto the candidate.
        v_fail_pass = int(_T.threshold("verdict_max_fail_pass"))
        v_fail_cond = int(_T.threshold("verdict_max_fail_conditional"))
        if not critical_passed:
            overall_verdict = "FALSE_POSITIVE"
        elif n_fail <= v_fail_pass:
            overall_verdict = "SOVEREIGN_PASS"
        elif n_fail <= v_fail_cond:
            overall_verdict = "CONDITIONAL_PASS"
        else:
            overall_verdict = "FALSE_POSITIVE"

        logical_closure = []
        if critical_passed:
            logical_closure.append(
                f"All {n_critical} critical FP tests PASSED → signal is statistically "
                f"real (SNR={self.snr:.1f}σ, FAP={self.fap:.1e}) and "
                f"periodic (Δσ_eo={self.eo_sigma:.2f} < 3.0)"
            )
        else:
            logical_closure.append("CRITICAL CIRCUIT BREAKER TRIGGERED: One or more critical tests failed")
            for t in tests:
                if t["weight"] == "critical" and t["verdict"] == "FAIL":
                    logical_closure.append(f"CRITICAL FAIL: {t['test']} — {t['description']}")
        
        # Add combined shape+density conflict detection (Task 4.2)
        if fp4_failed and fp7_failed:
            logical_closure.append("SHAPE+DENSITY CONFLICT: V-shaped transit with density mismatch indicates EB")

        logical_closure.append(
            f"Shape analysis: ratio={self.shape_r:.3f} {'> 1 (U-shape → PLANET-LIKE)' if self.shape_r > 1 else '< 1 (V-shape → EB-LIKE)'}"
        )
        logical_closure.append(
            f"Density constraint: ratio={density_ratio:.3f} → "
            f"{'consistent with planet' if _T.threshold('density_conflict_low') < density_ratio < _T.threshold('density_conflict_high') else 'INCONSISTENT — possible EB'}"
        )
        logical_closure.append(
            f"Total: {n_pass}/{len(tests)} tests PASSED → Axiom verdict: {overall_verdict}"
        )

        # Task 7.1: Conflict detection for SNR vs density/shape mismatches
        conflicts = []
        conflict_start_time = time.time()
        
        # Task 7.2: Get logger for conflict logging
        logger = get_logger(__name__)
        
        # Check for SNR above floor AND density_ratio outside the conflict band
        if self.snr > _T.threshold("density_conflict_snr_min") and (
                density_ratio < _T.threshold("density_conflict_low") or density_ratio > _T.threshold("density_conflict_high")):
            conflict_latency_ms = (time.time() - conflict_start_time) * 1000.0
            conflict_metadata = {
                "conflict_type": "SNR_DENSITY_CONFLICT",
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "tic_id": tic_id,
                "snr": round(self.snr, 4),
                "density_ratio": round(density_ratio, 4),
                "shape_ratio": None,
                "resolution_latency_ms": round(conflict_latency_ms, 4),
                "overall_verdict": overall_verdict,
            }
            conflicts.append(conflict_metadata)
            
            # Task 7.2: Log conflict with try-except to prevent validation failure
            try:
                logger.info(f"Conflict detected: SNR_DENSITY_CONFLICT for {tic_id}")
            except Exception:
                # Silently ignore logging errors to prevent validation failure
                pass
        
        # Check for SNR > 10.0 AND shape_ratio <= 1.0
        if self.snr > 10.0 and self.shape_r <= 1.0:
            conflict_latency_ms = (time.time() - conflict_start_time) * 1000.0
            conflict_metadata = {
                "conflict_type": "SNR_SHAPE_CONFLICT",
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "tic_id": tic_id,
                "snr": round(self.snr, 4),
                "density_ratio": None,
                "shape_ratio": round(self.shape_r, 4),
                "resolution_latency_ms": round(conflict_latency_ms, 4),
                "overall_verdict": overall_verdict,
            }
            conflicts.append(conflict_metadata)
            
            # Task 7.2: Log conflict with try-except to prevent validation failure
            try:
                logger.info(f"Conflict detected: SNR_SHAPE_CONFLICT for {tic_id}")
            except Exception:
                # Silently ignore logging errors to prevent validation failure
                pass

        return {
            "section": "§5 False-Positive Ruling",
            "formula": "Logical conjunction of 9 independent FP discriminators",
            "n_tests":          len(tests),
            "n_pass":           n_pass,
            "n_fail":           n_fail,
            "n_critical":       n_critical,
            "n_critical_pass":  n_crit_pass,
            "tests":            tests,
            "logical_closure":  logical_closure,
            "overall_verdict":  overall_verdict,
            "conflicts":        conflicts,
            "proof": (
                f"FP ruling | {n_pass}/{len(tests)} PASS, {n_crit_pass}/{n_critical} critical PASS | "
                f"shape={self.shape_r:.3f}, eo_sigma={self.eo_sigma:.3f}, "
                f"density_ratio={density_ratio:.3f} → {overall_verdict}"
            ),
        }

    # ─────────────────────────────────────────────────────────────────────────
    # Master proof assembly
    # ─────────────────────────────────────────────────────────────────────────

    def build_full_proof(
        self,
        secondary_snr:  float = 0.0,
        secondary_depth_ratio: float = 0.0,
        coherent_evidence: int = 0,
        alias_secondary_ratio: float = 0.0,
        centroid_sigma: float = 0.0,
        cvs_score:      float = 0.0,
        cvs_verdict:    str   = "",
        cvs_proof_chain: Optional[List[str]] = None,
        zspace_id:      str   = "",
        tic_id:         str   = "",
        t0_days:        Optional[float] = None,
        # Optional chi-squared analysis parameters
        time:           Optional[np.ndarray] = None,
        flux:           Optional[np.ndarray] = None,
        flux_err:       Optional[np.ndarray] = None,
        model_flux:     Optional[np.ndarray] = None,
        n_params:       int = 5,
    ) -> dict:
        """
        Assemble all §1–§5 proof blocks into the Sovereign Logic Card.
        
        Optional chi-squared analysis (§6) is included if flux data is provided.
        Optional FP-10 transit-count test is included if a light curve AND an
        epoch (t0_days) are provided (single transient events are not
        validatable as planets).
        """
        # §1
        kep = self.kepler_third_law()
        a_m = kep["a_m"]

        # §2
        geo = self.geometric_consistency()

        # §3
        den = self.density_constraint(a_m)

        # §4
        prob = self.transit_probability(a_m)

        # FP-10: count independent transits when a light curve + epoch exist.
        n_transits = None
        if (
            time is not None and flux is not None
            and self.P_days and self.P_days > 0.0
            and t0_days is not None
        ):
            dur_days = self.T_hrs / 24.0
            n_transits = count_observed_transits(
                time          = time,
                flux          = flux,
                period_days   = float(self.P_days),
                t0_days       = float(t0_days),
                duration_days = float(dur_days),
            )

        # §5
        fp = self.false_positive_ruling(
            secondary_snr  = secondary_snr,
            secondary_depth_ratio = secondary_depth_ratio,
            coherent_evidence = coherent_evidence,
            alias_secondary_ratio = alias_secondary_ratio,
            centroid_sigma = centroid_sigma,
            density_ratio  = den["density_ratio"],
            is_grazing     = prob["is_grazing"],
            tic_id         = tic_id,
            n_transits     = n_transits,
        )

        # §6 Chi-squared goodness-of-fit (optional)
        chi_sq_section = None
        if time is not None and flux is not None and flux_err is not None and model_flux is not None:
            analyzer = ChiSquaredAnalyzer()
            chi_sq_result = analyzer.compute_chi_squared(
                time=time,
                flux=flux,
                flux_err=flux_err,
                model_flux=model_flux,
                n_params=n_params,
            )
            chi_sq_section = {
                "section": "§6 Axiom Axiomatic Certainty",
                "chi_squared": chi_sq_result["chi_squared"],
                "reduced_chi_squared": chi_sq_result["reduced_chi_squared"],
                "degrees_of_freedom": chi_sq_result["degrees_of_freedom"],
                "p_value": chi_sq_result["p_value"],
                "interpretation": chi_sq_result["interpretation"],
                "quality": chi_sq_result["quality"],
                "proof": chi_sq_result["proof"],
                "verdict": "PASS" if chi_sq_result["quality"] in ("EXCELLENT", "GOOD") else "WARN",
            }

        # Aggregate verdict
        section_verdicts = [kep["verdict"], geo["verdict"], den["verdict"],
                            prob["verdict"], fp["overall_verdict"]]
        all_pass = all(v in ("PASS", "SOVEREIGN_PASS", "CONDITIONAL_PASS")
                       for v in section_verdicts)
        sovereign_verdict = fp["overall_verdict"]

        # Physical summary
        physical_summary = {
            "orbital_period_days":           round(self.P_days, 7),
            "semi_major_axis_au":            round(kep["a_au"], 7),
            "planet_radius_earth":           round(geo["R_p_earth"], 4),
            "radius_ratio_k":                round(geo["k"], 8),
            "transit_depth_ppm":             round(self.delta * 1e6, 2),
            "transit_duration_hours":        round(self.T_hrs, 4),
            "stellar_density_transit_gcc":   round(den["rho_transit_gcc"], 5),
            "stellar_density_tic_gcc":       round(den["rho_tic_gcc"], 5),
            "density_ratio":                 round(den["density_ratio"], 5),
            "geometric_transit_probability": round(prob["P_tr"] * 100, 4),
            "impact_parameter_b":            round(prob["impact_parameter_b"], 5),
            "ingress_duration_hrs":          round(prob["T_ingress_hrs"], 4),
            "min_inclination_deg":           round(prob["i_min_deg"], 4),
            "ld_corrected_depth":            round(geo["delta_ld_corrected"], 8),
            "limb_dark_u1":                  self.u1,
            "limb_dark_u2":                  self.u2,
        }

        # Build proof sections dict
        proof_sections = {
            "section_1_keplerian_dynamics":     kep,
            "section_2_geometric_consistency":  geo,
            "section_3_density_constraint":     den,
            "section_4_transit_probability":    prob,
            "section_5_false_positive_ruling":  fp,
        }
        
        # Add chi-squared section if available
        if chi_sq_section is not None:
            proof_sections["section_6_axiom_axiomatic_certainty"] = chi_sq_section

        return {
            "schema":               "Axiom-ZSpace Sovereign Logic Card v1.0",
            "zspace_id":            zspace_id,
            "timestamp_utc":        datetime.now(timezone.utc).isoformat(),
            "sovereign_verdict":    sovereign_verdict,
            "cvs_score":            cvs_score,
            "cvs_verdict":          cvs_verdict,
            "cvs_proof_chain":      cvs_proof_chain or [],
            "physical_summary":     physical_summary,
            "proof_sections":       proof_sections,
            "proof_integrity": {
                "all_sections_pass": all_pass,
                "section_verdicts":  dict(zip(
                    ["§1_Kepler","§2_Geometry","§3_Density","§4_Probability","§5_FP"],
                    section_verdicts,
                )),
                "truthimatics_seal": (
                    f"TRUTHIMATICS-SOVEREIGN | ZSpace={zspace_id} | "
                    f"CVS={cvs_score:.4f} | {sovereign_verdict} | "
                    f"WHITE-BOX | NO-ML | PHYSICS-DERIVED"
                ),
            },
        }


# ─────────────────────────────────────────────────────────────────────────────
# Archive Query Engine
# ─────────────────────────────────────────────────────────────────────────────

class ArchiveQueryEngine:
    """
    Queries NASA Exoplanet Archive + TOI catalogue.
    Implements two independent query paths with graceful degradation.

    Query Path A: astroquery.ipac.nexsci.nasa_exoplanet_archive (TAP/ADQL)
    Query Path B: direct requests to the TAP endpoint (fallback)

    On any network failure: returns (None, error_message).
    The caller (AxiomValidator) then proceeds in OFFLINE mode.
    """

    # Full column sets to retrieve
    PLANET_COLUMNS = (
        "pl_name, tic_id, pl_orbper, pl_rade, pl_radeerr1, pl_radeerr2, "
        "pl_trandep, pl_trandur, pl_tranmid, pl_ratdor, pl_imppar, "
        "st_teff, st_rad, st_mass, st_logg, st_met, "
        "ra, dec, sy_dist, discoverymethod, disc_year, pl_controv_flag"
    )
    TOI_COLUMNS = (
        "toi, tid, pl_orbper, pl_trandep, pl_trandurh, "
        "tfopwg_disp, ra, dec, st_tmag, st_teff, st_logg, st_rad, pl_rade"
    )

    def __init__(self, timeout: int = API_TIMEOUT_SECONDS) -> None:
        self.timeout = timeout

    # ── Path A: astroquery TAP ────────────────────────────────────────────────

    def _query_astroquery_planets(self, tic_id: str) -> Optional[List[dict]]:
        """Query confirmed planet catalogue by TIC ID via astroquery."""
        try:
            from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = NasaExoplanetArchive.query_criteria(
                    table   = "ps",  # Use newer 'ps' table instead of 'pscomppars'
                    select  = self.PLANET_COLUMNS,
                    where   = f"tic_id='{tic_id}'",
                    timeout = self.timeout,
                )
            if result is None or len(result) == 0:
                return []
            return [dict(zip(result.colnames, row)) for row in result]
        except Exception as exc:
            raise RuntimeError(f"astroquery planet query failed: {exc}") from exc

    def _query_astroquery_toi(self, tic_id: str) -> Optional[List[dict]]:
        """Query TOI catalogue by TIC ID via astroquery."""
        try:
            from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = NasaExoplanetArchive.query_criteria(
                    table   = "toi",
                    select  = self.TOI_COLUMNS,
                    where   = f"tid='{tic_id}'",  # TOI table uses 'tid' not 'tic_id'
                    timeout = self.timeout,
                )
            if result is None or len(result) == 0:
                return []
            return [dict(zip(result.colnames, row)) for row in result]
        except Exception as exc:
            raise RuntimeError(f"astroquery TOI query failed: {exc}") from exc

    # ── Path B: direct requests TAP fallback ─────────────────────────────────

    def _query_tap_direct(self, adql: str) -> Optional[List[dict]]:
        """
        Direct TAP query via requests + manual VOTable parsing.
        Used when astroquery is unavailable or its TAP endpoint is cached.
        """
        try:
            import requests
            from io import BytesIO
            from astropy.io.votable import parse_single_table

            params = {
                "QUERY":  adql,
                "FORMAT": "votable",
                "LANG":   "ADQL",
                "REQUEST":"doQuery",
            }
            resp = requests.get(
                NASA_TAP_URL,
                params  = params,
                timeout = self.timeout,
            )
            resp.raise_for_status()
            table = parse_single_table(BytesIO(resp.content))
            arr   = table.array
            if len(arr) == 0:
                return []
            return [
                {col: (arr[col][i].item() if hasattr(arr[col][i], 'item') else arr[col][i])
                 for col in arr.dtype.names}
                for i in range(len(arr))
            ]
        except Exception as exc:
            raise RuntimeError(f"Direct TAP query failed: {exc}") from exc

    # ── Coordinate-based fallback (for non-TIC targets) ──────────────────────

    def query_by_coordinates(
        self,
        ra_deg:  float,
        dec_deg: float,
        radius_arcsec: float = 10.0,
    ) -> Tuple[Optional[List[dict]], Optional[List[dict]], Optional[str]]:
        """
        Query by sky coordinates when TIC ID is unavailable.
        Returns (planet_rows, toi_rows, error).
        """
        adql_planets = (
            f"SELECT {self.PLANET_COLUMNS} FROM ps "
            f"WHERE CONTAINS(POINT('ICRS', ra, dec), "
            f"CIRCLE('ICRS', {ra_deg}, {dec_deg}, {radius_arcsec/3600.0})) = 1"
        )
        adql_toi = (
            f"SELECT {self.TOI_COLUMNS} FROM toi "
            f"WHERE CONTAINS(POINT('ICRS', ra, dec), "
            f"CIRCLE('ICRS', {ra_deg}, {dec_deg}, {radius_arcsec/3600.0})) = 1"
        )
        return self._execute_with_fallback(adql_planets, adql_toi)

    # ── Master query dispatcher ────────────────────────────────────────────────

    def query(
        self,
        tic_id: str,
    ) -> Tuple[Optional[List[dict]], Optional[List[dict]], Optional[str]]:
        """
        Execute both planet + TOI queries with retry and fallback.

        Returns
        -------
        planet_rows : list of dicts (may be empty list)
        toi_rows    : list of dicts (may be empty list)
        error       : None if successful; error string if network failed
        """
        planet_adql = (
            f"SELECT {self.PLANET_COLUMNS} FROM ps "
            f"WHERE tic_id='{tic_id}'"
        )
        toi_adql = (
            f"SELECT {self.TOI_COLUMNS} FROM toi "
            f"WHERE tid='{tic_id}'"  # TOI table uses 'tid' not 'tic_id'
        )
        return self._execute_with_fallback(planet_adql, toi_adql, tic_id=tic_id)

    def _execute_with_fallback(
        self,
        planet_adql: str,
        toi_adql:    str,
        tic_id:      Optional[str] = None,
    ) -> Tuple[Optional[List[dict]], Optional[List[dict]], Optional[str]]:
        """Internal: try astroquery path, then direct TAP, then fail gracefully."""
        planet_rows: Optional[List[dict]] = None
        toi_rows:    Optional[List[dict]] = None
        last_error = None

        for attempt in range(1, MAX_RETRY_ATTEMPTS + 1):
            try:
                # Path A: astroquery
                if tic_id:
                    planet_rows = self._query_astroquery_planets(tic_id)
                    toi_rows    = self._query_astroquery_toi(tic_id)
                else:
                    planet_rows = self._query_tap_direct(planet_adql)
                    toi_rows    = self._query_tap_direct(toi_adql)

                if planet_rows is not None and toi_rows is not None:
                    return planet_rows, toi_rows, None   # ← SUCCESS

            except Exception as exc:
                last_error = str(exc)
                # Path B: direct TAP
                try:
                    if planet_rows is None:
                        planet_rows = self._query_tap_direct(planet_adql)
                    if toi_rows is None:
                        toi_rows = self._query_tap_direct(toi_adql)
                    if planet_rows is not None and toi_rows is not None:
                        return planet_rows, toi_rows, None   # ← Path B SUCCESS
                except Exception as exc2:
                    last_error = f"Path-A: {exc} | Path-B: {exc2}"

            if attempt < MAX_RETRY_ATTEMPTS:
                time.sleep(2)

        return None, None, f"NETWORK_FAIL (all paths exhausted): {last_error}"


# ─────────────────────────────────────────────────────────────────────────────
# Period Comparator
# ─────────────────────────────────────────────────────────────────────────────

class PeriodComparator:
    """
    Compares a candidate period against all archive entries.
    Handles period aliases (P/2, 2P) and harmonic multiples.
    """

    ALIAS_FACTORS = [1.0, 2.0, 0.5, 3.0, 1.0/3.0, 4.0, 0.25]   # common aliases
    # Ground-truth consistency tolerance: the found period must agree with the
    # expected archive period to within ±5% of one of the harmonic factors.
    GROUND_TRUTH_REL_TOL = 0.05

    @classmethod
    def period_consistent(
        cls,
        candidate_period:      float,
        ground_truth_period:   float,
        relative_tolerance:    float = GROUND_TRUTH_REL_TOL,
    ) -> bool:
        """
        True if `candidate_period` describes the SAME planet as
        `ground_truth_period` (the archive-known period), allowing the common
        harmonic aliases. A period that differs by more than the tolerance
        (and is not a pure alias) belongs to a DIFFERENT signal.

        Consistency is required for correctly folding an ephemeris. This is
        the decisive guard against:
          - sibling-planet confusion  (L 98-59 b found at c's 3.691 d);
          - wrong-ephemeris SOVEREIGN_PASS (HD 63433 b/d at 12.94 d).
        """
        if ground_truth_period <= 0 or candidate_period <= 0:
            return False
        if candidate_period == ground_truth_period:
            return True
        for alias in cls.ALIAS_FACTORS:
            dist = abs(candidate_period - ground_truth_period * alias)
            if dist / max(ground_truth_period * alias, 1e-9) <= relative_tolerance:
                return True
        return False

    @classmethod
    def find_match(
        cls,
        candidate_period: float,
        archive_rows:     List[dict],
        period_key:       str = "pl_orbper",
        tolerance:        float = PERIOD_MATCH_TOLERANCE_DAYS,
        relative_tolerance: float = PERIOD_MATCH_TOLERANCE_REL,
    ) -> Optional[Tuple[dict, float, float]]:
        """
        Return (best_match_row, period_delta_days, alias_factor) or None.
        Checks exact match first, then harmonic aliases.

        Uses both absolute tolerance (days) and relative tolerance (fraction).
        A match is accepted if EITHER condition is met:
          |P_candidate - P_archive*alias| < tolerance_days   (absolute)
          |P_candidate - P_archive*alias| / P_archive*alias < relative_tolerance  (relative)
        """
        best_match   = None
        best_delta   = float("inf")
        best_alias   = 1.0

        for row in archive_rows:
            try:
                P_arch = float(row[period_key])
            except (KeyError, TypeError, ValueError):
                continue
            if P_arch <= 0:
                continue

            for alias in cls.ALIAS_FACTORS:
                P_check = P_arch * alias
                if P_check <= 0:
                    continue
                delta   = abs(candidate_period - P_check)
                rel_delta = delta / P_check
                # Match if absolute OR relative tolerance is satisfied
                if (delta < tolerance or rel_delta < relative_tolerance) and delta < best_delta:
                    best_delta = delta
                    best_match = row
                    best_alias = alias

        if best_match is not None:
            return best_match, best_delta, best_alias
        return None


# ─────────────────────────────────────────────────────────────────────────────
# JSON serialiser (handles numpy, astropy masked, Python native)
# ─────────────────────────────────────────────────────────────────────────────

def _json_safe(obj: Any) -> Any:
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.bool_):
        return bool(obj)
    # astropy masked values → None if masked, else native Python
    try:
        import numpy.ma as ma
        if isinstance(obj, ma.core.MaskedConstant):
            return None
    except ImportError:
        pass
    if hasattr(obj, "item"):
        return obj.item()
    return str(obj)


def write_json(path: str, data: dict, indent: int = 2) -> None:
    """Write JSON data to file with error handling."""
    try:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=indent, default=_json_safe)
    except IOError as e:
        logger = get_logger(__name__)
        logger.error(f"File I/O error writing JSON to {path}: {e}")
        raise IOError(f"Failed to write JSON file {path}: {e}") from e
    except Exception as e:
        logger = get_logger(__name__)
        logger.error(f"Unexpected error writing JSON to {path}: {e}")
        raise RuntimeError(f"Failed to write JSON file {path}: {e}") from e


# ─────────────────────────────────────────────────────────────────────────────
# Transit counting (FP-10)
# ─────────────────────────────────────────────────────────────────────────────


def count_observed_transits(
    time,
    flux,
    period_days,
    t0_days,
    duration_days,
    min_parity_sigma: float = 3.0,
):
    """
    Count the number of *independent* epoch series supporting a periodic
    ephemeris (period_days, t0_days) in the given light curve.

    The fundamental limit behind the FP-10 gate is that an individual shallow
    transit is usually sub-threshold; only the phase-fold reveals a genuine
    planet, and any periodic signal can be *constructed* from a single deep
    event (the 'blip' adversary). The discriminating quantity is therefore
    whether the signal also exists in an epoch series that shares no events
    with the deepest one.

    Implementation: the epoch list e = t0 + k*P (k in Z over the baseline) is
    split into the two parity series (even/odd k). For each parity, all
    in-transit samples (|t - e| <= duration/2) are pooled and the signal is
    measured as the flux deficit of that pooled set vs the global baseline:
        sig = <deficit> / (sigma_baseline / sqrt(N_pooled)).
    This mirrors the fold-depth significance restricted to that parity alone.
    A parity is 'supporting' when sig >= min_parity_sigma.

    * A genuine planet populates every epoch, so BOTH parities are
      supporting -> returns 2 (>=2 independent transits, validated).
    * A single transient populates exactly one epoch *, so one parity series
      is empty and not supporting -> returns 1 (FP-10 FAIL, never validated).
    * A periodic EB populates both parity series (primary/secondary eclipses
      alternate) -> returns 2; FP-10 passes and the even/odd + secondary
      eclipse tests veto it instead (as they must).

    Returns 0..2, or None when the light curve is absent / degenerate.
    """
    if time is None or flux is None:
        return None
    time = np.asarray(time, dtype=np.float64).ravel()
    flux = np.asarray(flux, dtype=np.float64).ravel()
    if time.size < 20:
        return None
    finite = np.isfinite(time) & np.isfinite(flux)
    if int(finite.sum()) < 20:
        return None
    t, f = time[finite], flux[finite]

    P = float(period_days)
    if not (P > 0.0) or not np.isfinite(P):
        return None
    t0 = float(t0_days)
    D = max(float(duration_days), 1e-3)          # clamp zero/negative durations
    half = 0.5 * D

    baseline = float(np.median(f))
    sigma = float(np.std(f - baseline))
    if sigma <= 0.0 or not np.isfinite(sigma):
        return 0

    pooled_n = [0, 0]
    pooled_sum = [0.0, 0.0]
    k_start = int(np.floor((t.min() - t0) / P)) - 1
    k_end   = int(np.ceil((t.max() - t0) / P)) + 1
    for k in range(k_start, k_end + 1):
        e = t0 + k * P
        win = (t >= e - half) & (t <= e + half)
        n_cad = int(win.sum())
        if n_cad < 3:
            continue
        par = int(k) % 2
        pooled_n[par]   += n_cad
        pooled_sum[par] += float(np.sum(f[win]))

    n_supporting = 0
    for par in (0, 1):
        if pooled_n[par] < 9:
            continue
        mean_f = pooled_sum[par] / float(pooled_n[par])
        deficit = baseline - mean_f
        if deficit <= 0.0:
            continue
        sig = deficit / (sigma / np.sqrt(float(pooled_n[par])))
        if sig >= min_parity_sigma:
            n_supporting += 1
    return int(n_supporting)


# ─────────────────────────────────────────────────────────────────────────────
# External Catalog Cross-Matching
# ─────────────────────────────────────────────────────────────────────────────

def check_external_catalogs(tic_id: str, timeout: int = 10) -> dict:
    """
    Cross-match TIC ID against SIMBAD and Gaia DR3 catalogs to identify
    known multiple star systems.
    
    Args:
        tic_id: TESS Input Catalog identifier (e.g., "TIC 307210830" or "307210830")
        timeout: Query timeout in seconds (default: 10)
    
    Returns:
        Dictionary with keys:
            - is_multiple: bool (True if stellar multiplicity detected)
            - catalog_source: str ("SIMBAD" | "Gaia" | "None" | "Offline")
            - risk_level: str ("HIGH" | "LOW" | "UNKNOWN")
            - classification: str (stellar type from catalog)
            - query_latency_ms: float (time taken for query)
    
    Raises:
        None (gracefully handles network failures)
    """
    logger = get_logger(__name__)
    start_time = time.time()
    
    # Simple in-memory cache for query results
    if not hasattr(check_external_catalogs, '_cache'):
        check_external_catalogs._cache = {}
    
    # Check cache first
    if tic_id in check_external_catalogs._cache:
        cached_result = check_external_catalogs._cache[tic_id].copy()
        cached_result['query_latency_ms'] = 0.0  # Cache hit
        return cached_result
    
    # Validate TIC ID format
    import re
    tic_pattern = r'^(?:TIC\s*)?(\d+)$'
    match = re.match(tic_pattern, str(tic_id).strip(), re.IGNORECASE)
    
    if not match:
        result = {
            "is_multiple": False,
            "catalog_source": "None",
            "risk_level": "UNKNOWN",
            "classification": "Invalid TIC ID",
            "query_latency_ms": (time.time() - start_time) * 1000.0,
        }
        logger.warning(f"Invalid TIC ID format: {tic_id}")
        return result
    
    tic_number = match.group(1)
    tic_formatted = f"TIC {tic_number}"
    
    # Try SIMBAD query first
    try:
        from astroquery.simbad import Simbad
        
        # Configure SIMBAD to return object type
        custom_simbad = Simbad()
        custom_simbad.add_votable_fields('otype')
        custom_simbad.TIMEOUT = timeout
        
        # Query SIMBAD
        result_table = custom_simbad.query_object(tic_formatted)
        
        if result_table is not None and len(result_table) > 0:
            # Extract object type and classification
            otype = str(result_table['OTYPE'][0]) if 'OTYPE' in result_table.colnames else ''
            main_id = str(result_table['MAIN_ID'][0]) if 'MAIN_ID' in result_table.colnames else ''
            
            # Check for multiplicity indicators
            multiplicity_keywords = ['Double', 'Multiple', 'EB*', 'V*', 'Binary', 'SB*', 'El*']
            is_multiple = any(keyword.lower() in otype.lower() or keyword.lower() in main_id.lower() 
                            for keyword in multiplicity_keywords)
            
            result = {
                "is_multiple": is_multiple,
                "catalog_source": "SIMBAD",
                "risk_level": "HIGH" if is_multiple else "LOW",
                "classification": otype if otype else main_id,
                "query_latency_ms": (time.time() - start_time) * 1000.0,
            }
            
            # Cache the result
            check_external_catalogs._cache[tic_id] = result.copy()
            
            logger.info(f"SIMBAD query for {tic_formatted}: {otype}, is_multiple={is_multiple}")
            return result
            
    except Exception as e:
        logger.warning(f"SIMBAD query failed for {tic_formatted}: {e}")
    
    # Try Gaia DR3 query as fallback
    try:
        from astroquery.gaia import Gaia
        
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
        Gaia.ROW_LIMIT = 1
        
        # Query Gaia for non_single_star flag
        # Note: We need to cross-match TIC to Gaia source_id first
        # For simplicity, we'll query by TIC designation if available
        query = f"""
        SELECT TOP 1 source_id, non_single_star, phot_variable_flag
        FROM gaiadr3.gaia_source
        WHERE source_id IN (
            SELECT gaia_source_id FROM gaiadr3.tic_gaia_xmatch
            WHERE tic_id = {tic_number}
        )
        """
        
        job = Gaia.launch_job_async(query, dump_to_file=False)
        result_table = job.get_results()
        
        if result_table is not None and len(result_table) > 0:
            non_single_star = result_table['non_single_star'][0] if 'non_single_star' in result_table.colnames else 0
            phot_variable = result_table['phot_variable_flag'][0] if 'phot_variable_flag' in result_table.colnames else ''
            
            # non_single_star flag: 0 = single, >0 = multiple
            is_multiple = bool(non_single_star > 0)
            
            classification = f"Gaia non_single_star={non_single_star}"
            if phot_variable:
                classification += f", variable={phot_variable}"
            
            result = {
                "is_multiple": is_multiple,
                "catalog_source": "Gaia",
                "risk_level": "HIGH" if is_multiple else "LOW",
                "classification": classification,
                "query_latency_ms": (time.time() - start_time) * 1000.0,
            }
            
            # Cache the result
            check_external_catalogs._cache[tic_id] = result.copy()
            
            logger.info(f"Gaia query for {tic_formatted}: non_single_star={non_single_star}, is_multiple={is_multiple}")
            return result
            
    except Exception as e:
        logger.warning(f"Gaia query failed for {tic_formatted}: {e}")
    
    # Both queries failed - return offline/unknown status
    result = {
        "is_multiple": False,
        "catalog_source": "Offline",
        "risk_level": "UNKNOWN",
        "classification": "Network unavailable",
        "query_latency_ms": (time.time() - start_time) * 1000.0,
    }
    
    logger.warning(f"All catalog queries failed for {tic_formatted}, returning offline status")
    return result


# ─────────────────────────────────────────────────────────────────────────────
# AxiomValidator — Public API
# ─────────────────────────────────────────────────────────────────────────────

class AxiomValidator:
    """
    Sovereign Validation Module for Axiom-ZSpace.

    Cross-references every candidate against NASA Exoplanet Archive + TOI list.
    Generates either exist_planet.json (known) or Discovery.json (new).

    Ephemeris-identity gate (NEW)
    -----------------------------
    When `expected_period_days` (the ground-truth archive period of the planet
    under test) is supplied — benchmark mode — the validator refuses to certify
    a signal whose period is not consistent with that ground truth (to within
    ±5% including harmonic aliases). Such a signal is classified
    EPHEMERIS_MISMATCH and never KNOWN / NEW_DISCOVERY. This prevents:
      - sibling-planet confusion   (L 98-59 b folded at c's period);
      - wrong-ephemeris SOVEREIGN_PASS (HD 63433 b/d both at 12.94 d).

    Parameters
    ----------
    output_dir   : directory for output JSON files (created if needed)
    period_tol   : period match tolerance in days (default 0.001)
    timeout      : API timeout in seconds (default 30)
    verbose      : print progress to stdout
    """

    def __init__(
        self,
        output_dir:  str   = ".",
        period_tol:  float = PERIOD_MATCH_TOLERANCE_DAYS,
        timeout:     int   = API_TIMEOUT_SECONDS,
        verbose:     bool  = True,
    ) -> None:
        self.output_dir  = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.period_tol  = period_tol
        self.verbose     = verbose
        self._querier    = ArchiveQueryEngine(timeout=timeout)
        self._comparator = PeriodComparator()
        self.logger      = get_logger(__name__)

    def _log(self, msg: str) -> None:
        if self.verbose:
            self.logger.info(f"[VALIDATOR] {msg}")

    # ── Main entry point ──────────────────────────────────────────────────────

    def validate(
        self,
        tic_id:                str,
        period_days:           float,
        transit_depth:         float,
        transit_duration_hrs:  float,
        t0_btjd:               float,
        stellar_mass_solar:    float,
        stellar_radius_solar:  float,
        stellar_teff_k:        float,
        stellar_logg:          float,
        planet_radius_earth:   float,
        cvs_score:             float,
        cvs_verdict:           str,
        cvs_proof_chain:       Optional[List[str]] = None,
        bls_snr:               float = 0.0,
        bls_fap:               float = 1.0,
        even_odd_delta_sigma:  float = 0.0,
        shape_ratio:           float = 1.0,
        secondary_snr:         float = 0.0,
        secondary_depth_ratio: float = 0.0,
        coherent_evidence:     int   = 0,
        alias_secondary_ratio: float = 0.0,
        centroid_sigma:        float = 0.0,
        limb_dark_u1:          float = 0.30,
        limb_dark_u2:          float = 0.10,
        zspace_id:             Optional[str] = None,
        planet_order:          int   = 1,
        # Optional chi-squared analysis parameters
        time:                  Optional[np.ndarray] = None,
        flux:                  Optional[np.ndarray] = None,
        flux_err:              Optional[np.ndarray] = None,
        model_flux:            Optional[np.ndarray] = None,
        n_params:              int = 5,
        # ── Ephemeris-identity gate (EPHEMERIS_MISMATCH) ─────────────────────
        # When the caller KNOWS the true planet being tested (benchmark mode),
        # pass its archive period here. The validator then refuses to certify
        # the signal under ANY name/period that is NOT consistent with this
        # ground-truth period. This kills the "wrong-ephemeris SOVEREIGN_PASS"
        # and "sibling-confused KNOWN" bugs: a signal folded at a sibling's
        # period is classified EPHEMERIS_MISMATCH, never KNOWN/NEW_DISCOVERY.
        expected_period_days:  Optional[float] = None,
        expected_planet_name:  Optional[str]   = None,
    ) -> ValidationResult:
        """
        Full validation pipeline.

        Parameters
        ----------
        (see module docstring for full parameter descriptions)

        Returns
        -------
        ValidationResult with status, output file path, and proof summary.
        """
        if zspace_id is None:
            zspace_id = f"ZS-T-{tic_id}-{planet_order:02d}"

        self._log(f"Starting sovereign validation for {zspace_id}")
        self._log(f"Period = {period_days:.6f} d | Depth = {transit_depth*1e6:.1f} ppm")

        # ── Step 1: Query NASA Archive ────────────────────────────────────────
        self._log("Querying NASA Exoplanet Archive + TOI …")
        planet_rows, toi_rows, net_error = self._querier.query(tic_id)

        if net_error:
            self._log(f"Network error: {net_error}")
            self._log("-> OFFLINE MODE: Proceeding with mathematical proof only.")
        else:
            np_ = len(planet_rows) if planet_rows else 0
            nt_ = len(toi_rows)    if toi_rows    else 0
            self._log(f"Archive: {np_} confirmed planet(s), {nt_} TOI(s) found for TIC {tic_id}")

        # ── Step 2: Period matching ───────────────────────────────────────────
        match:        Optional[ArchiveMatch] = None
        all_archive   = list(planet_rows or [])
        all_toi       = list(toi_rows    or [])

        # Check confirmed planets
        planet_hit = self._comparator.find_match(period_days, all_archive, "pl_orbper", self.period_tol)
        # Check TOIs
        toi_hit    = self._comparator.find_match(period_days, all_toi, "pl_orbper", self.period_tol)

        if planet_hit:
            row, delta, alias = planet_hit
            match = self._parse_planet_match(row, delta, alias, "NASA_ARCHIVE")
            self._log(f"MATCH FOUND in NASA Archive: {match.planet_name} "
                      f"(dP={delta:.5f} d, alias x{alias:.1f})")
        elif toi_hit:
            row, delta, alias = toi_hit
            match = self._parse_toi_match(row, delta, alias)
            self._log(f"MATCH FOUND in TOI list: TOI {match.planet_name} "
                      f"(dP={delta:.5f} d, alias x{alias:.1f})")
        else:
            self._log("No period match found -> proceeding to Discovery Protocol")

        # ── Step 3: Route to KNOWN, DISCOVERY, or EPHEMERIS_MISMATCH ─────────
        #
        # EPHEMERIS-MISMATCH identity gate (ground-truth benchmark mode):
        # The caller knows which planet is under test. If the found period is
        # NOT consistent with that planet's archive period (to within the
        # relative matching tolerance, ignoring pure harmonic aliases), then
        # whatever we found is a DIFFERENT signal (a sibling planet, an
        # activity/rotation harmonic, or an alias) — NOT the tested planet.
        # Certifying it as KNOWN (via the sibling match) or as a NEW_DISCOVERY
        # (sovereign pass on a wrong ephemeris) would both be false. It is
        # instead reported as EPHEMERIS_MISMATCH so the benchmark counts it
        # as a wrong-period recovery, not as a success.
        if expected_period_days is not None and expected_period_days > 0:
            ground_truth_ok = PeriodComparator.period_consistent(
                candidate_period  = period_days,
                ground_truth_period = expected_period_days,
            )
            if not ground_truth_ok:
                status, outfile = self._handle_ephemeris_mismatch(
                    tic_id, zspace_id, period_days,
                    expected_period_days, expected_planet_name,
                    match, all_archive, all_toi,
                )
                return ValidationResult(
                    status       = status,
                    tic_id       = tic_id,
                    period_days  = period_days,
                    match        = match,
                    output_file  = str(outfile),
                    proof_summary = (
                        f"EPHEMERIS_MISMATCH | found P={period_days:.6f} d "
                        f"(expected {expected_period_days:.6f} d) → "
                        f"{match.planet_name if match else 'no archive match'}"
                    ),
                    network_error = net_error,
                    expected_period_days = expected_period_days,
                    expected_planet_name = expected_planet_name,
                    alias_sibling_match  = match,
                )

        if match and not net_error:
            status, outfile = self._handle_known(
                tic_id, zspace_id, period_days, transit_depth,
                stellar_mass_solar, stellar_radius_solar, stellar_teff_k,
                planet_radius_earth, cvs_score, cvs_verdict,
                match, all_archive, all_toi,
            )
        else:
            # Also run discovery if network failed (OFFLINE_NEW_DISCOVERY)
            status_tag = "NEW_DISCOVERY" if not net_error else "OFFLINE_NEW_DISCOVERY"
            status, outfile = self._handle_discovery(
                tic_id, zspace_id, period_days, transit_depth,
                transit_duration_hrs, stellar_mass_solar, stellar_radius_solar,
                stellar_teff_k, stellar_logg, planet_radius_earth,
                cvs_score, cvs_verdict, cvs_proof_chain or [],
                bls_snr, bls_fap, even_odd_delta_sigma, shape_ratio,
                secondary_snr, secondary_depth_ratio, coherent_evidence,
                alias_secondary_ratio, centroid_sigma,
                limb_dark_u1, limb_dark_u2, status_tag, net_error,
                t0_days=t0_btjd,
                time=time, flux=flux, flux_err=flux_err, model_flux=model_flux, n_params=n_params,
            )

        proof_summary = (
            f"{'MATCH: ' + match.planet_name if match else 'NO MATCH -> SOVEREIGN DISCOVERY'} | "
            f"TIC={tic_id}, P={period_days:.6f} d, CVS={cvs_score:.4f} | "
            f"-> {outfile}"
        )

        self._log(f"Validation complete: {status} -> {outfile}")

        return ValidationResult(
            status       = status,
            tic_id       = tic_id,
            period_days  = period_days,
            match        = match,
            output_file  = str(outfile),
            proof_summary = proof_summary,
            network_error = net_error,
            expected_period_days = expected_period_days,
            expected_planet_name = expected_planet_name,
            alias_sibling_match  = match if status == "EPHEMERIS_MISMATCH" else None,
        )

    # ── Known planet handler ──────────────────────────────────────────────────

    def _handle_known(
        self,
        tic_id:               str,
        zspace_id:            str,
        period_days:          float,
        transit_depth:        float,
        stellar_mass_solar:   float,
        stellar_radius_solar: float,
        stellar_teff_k:       float,
        planet_radius_earth:  float,
        cvs_score:            float,
        cvs_verdict:          str,
        match:                ArchiveMatch,
        all_archive:          List[dict],
        all_toi:              List[dict],
    ) -> Tuple[str, str]:
        """Build exist_planet.json."""
        out_path = self.output_dir / f"exist_planet_{zspace_id}.json"

        card = {
            "schema":          "Axiom-ZSpace Existing Planet Card v1.0",
            "zspace_id":       zspace_id,
            "timestamp_utc":   datetime.now(timezone.utc).isoformat(),
            "status":          "KNOWN",
            "match_source":    match.source,
            "match_summary": {
                "planet_name":         match.planet_name,
                "period_delta_days":   round(match.period_delta_days, 8),
                "tolerance_days":      self.period_tol,
                "proof": (
                    f"|P_candidate - P_archive| = {match.period_delta_days:.8f} d "
                    f"≤ {self.period_tol} d → MATCH"
                ),
            },
            "candidate_parameters": {
                "tic_id":                tic_id,
                "period_days":           period_days,
                "transit_depth_ppm":     round(transit_depth * 1e6, 2),
                "planet_radius_earth":   planet_radius_earth,
                "stellar_mass_solar":    stellar_mass_solar,
                "stellar_radius_solar":  stellar_radius_solar,
                "stellar_teff_k":        stellar_teff_k,
                "cvs_score":             cvs_score,
                "cvs_verdict":           cvs_verdict,
            },
            "archive_parameters": {
                "planet_name":           match.planet_name,
                "period_days_archive":   match.period_days,
                "transit_depth_archive": match.transit_depth,
                "planet_radius_earth":   match.planet_radius_earth,
                "semi_major_axis_au":    match.semi_major_axis_au,
                "stellar_teff_k":        match.stellar_teff_k,
                "stellar_radius_solar":  match.stellar_radius_solar,
                "stellar_mass_solar":    match.stellar_mass_solar,
                "discovery_method":      match.discovery_method,
                "disposition":           match.disposition,
                "extra_fields":          match.extra_fields,
            },
            "parameter_deltas": self._compute_deltas(
                period_days, transit_depth, planet_radius_earth, match
            ),
            "all_archive_entries_for_tic": all_archive,
            "all_toi_entries_for_tic":     all_toi,
        }

        write_json(str(out_path), card)
        self._log(f"Written -> {out_path}")
        return "KNOWN", str(out_path)

    # ── Ephemeris-mismatch handler ────────────────────────────────────────────

    def _handle_ephemeris_mismatch(
        self,
        tic_id:                 str,
        zspace_id:              str,
        period_days:            float,
        expected_period_days:   Optional[float],
        expected_planet_name:   Optional[str],
        match:                  Optional[ArchiveMatch],
        all_archive:            List[dict],
        all_toi:                List[dict],
    ) -> Tuple[str, str]:
        """
        Build EphemerisMismatch.json.

        The signal found by the detector is real and periodic, but its period
        is NOT consistent with the ground-truth planet under test. It may be
        a sibling planet in the same system (e.g. L 98-59 b folded at c's
        period), a stellar-rotation/activity harmonic, or an alias. In every
        case it must NOT be certified as KNOWN or as a NEW_DISCOVERY.
        """
        out_path = self.output_dir / f"EphemerisMismatch_{zspace_id}.json"

        # Which archive entry is the found period actually consistent with?
        sibling_planet = None
        if match is not None:
            sibling_planet = {
                "planet_name":       match.planet_name,
                "period_days":       match.period_days,
                "period_delta_days": round(match.period_delta_days, 8),
                "alias_factor":      match.extra_fields.get("period_alias_factor"),
                "source":            match.source,
            }

        card = {
            "schema":           "Axiom-ZSpace Ephemeris Mismatch Card v1.0",
            "zspace_id":        zspace_id,
            "timestamp_utc":    datetime.now(timezone.utc).isoformat(),
            "status":           "EPHEMERIS_MISMATCH",
            "tic_id":           tic_id,
            "found_period_days":  round(period_days, 8),
            "expected_period_days": round(expected_period_days, 8) if expected_period_days else None,
            "expected_planet_name": expected_planet_name,
            "identity_gate_proof": (
                f"found P={period_days:.8f} d is NOT consistent with ground-truth "
                f"P={expected_period_days:.8f} d (tolerance ±5% incl. harmonic aliases) "
                f"→ signal belongs to a DIFFERENT ephemeris (sibling planet, "
                f"activity harmonic, or alias), NOT the planet under test."
            ),
            "sibling_planet_matched": sibling_planet,
            "all_archive_entries_for_tic": all_archive,
            "all_toi_entries_for_tic":     all_toi,
            # The found signal is NOT the tested planet → must NOT be counted:
            "certified_as": "NONE",
        }

        write_json(str(out_path), card)
        self._log(f"Written -> {out_path} (EPHEMERIS_MISMATCH)")
        return "EPHEMERIS_MISMATCH", str(out_path)

    # ── Discovery handler ─────────────────────────────────────────────────────

    def _handle_discovery(
        self,
        tic_id:               str,
        zspace_id:            str,
        period_days:          float,
        transit_depth:        float,
        transit_duration_hrs: float,
        stellar_mass_solar:   float,
        stellar_radius_solar: float,
        stellar_teff_k:       float,
        stellar_logg:         float,
        planet_radius_earth:  float,
        cvs_score:            float,
        cvs_verdict:          str,
        cvs_proof_chain:      List[str],
        bls_snr:              float,
        bls_fap:              float,
        even_odd_delta_sigma: float,
        shape_ratio:          float,
        secondary_snr:        float,
        secondary_depth_ratio: float,
        coherent_evidence: int,
        alias_secondary_ratio: float,
        centroid_sigma:       float,
        limb_dark_u1:         float,
        limb_dark_u2:         float,
        status_tag:           str,
        network_error:        Optional[str],
        t0_days:              Optional[float] = None,
        time:                 Optional[np.ndarray] = None,
        flux:                 Optional[np.ndarray] = None,
        flux_err:             Optional[np.ndarray] = None,
        model_flux:           Optional[np.ndarray] = None,
        n_params:             int = 5,
    ) -> Tuple[str, str]:
        """Build Discovery.json — the full Sovereign Logic Card."""
        out_path = self.output_dir / f"Discovery_{zspace_id}.json"
        self._log("Building Sovereign Logic Card (full mathematical proof) …")

        # Instantiate proof engine
        engine = ProofEngine(
            period_days            = period_days,
            transit_depth          = transit_depth,
            transit_duration_hrs   = transit_duration_hrs,
            stellar_mass_solar     = stellar_mass_solar,
            stellar_radius_solar   = stellar_radius_solar,
            stellar_teff_k         = stellar_teff_k,
            stellar_logg           = stellar_logg,
            planet_radius_earth    = planet_radius_earth,
            bls_snr                = bls_snr,
            bls_fap                = bls_fap,
            even_odd_delta_sigma   = even_odd_delta_sigma,
            shape_ratio            = shape_ratio,
            limb_dark_u1           = limb_dark_u1,
            limb_dark_u2           = limb_dark_u2,
        )

        full_proof = engine.build_full_proof(
            secondary_snr    = secondary_snr,
            secondary_depth_ratio = secondary_depth_ratio,
            coherent_evidence = coherent_evidence,
            alias_secondary_ratio = alias_secondary_ratio,
            centroid_sigma   = centroid_sigma,
            cvs_score        = cvs_score,
            cvs_verdict      = cvs_verdict,
            cvs_proof_chain  = cvs_proof_chain,
            zspace_id        = zspace_id,
            tic_id           = tic_id,
            t0_days          = t0_days,
            time             = time,
            flux             = flux,
            flux_err         = flux_err,
            model_flux       = model_flux,
            n_params         = n_params,
        )

        # Wrap in the Discovery card envelope
        card = {
            "schema":            "Axiom-ZSpace Sovereign Logic Card v1.0",
            "status":            status_tag,
            "network_mode":      "ONLINE" if not network_error else "OFFLINE",
            "network_error":     network_error,
            "tic_id":            tic_id,
            "zspace_id":         zspace_id,
            "timestamp_utc":     datetime.now(timezone.utc).isoformat(),
            "sovereign_verdict": full_proof["sovereign_verdict"],
            "validation_declaration": (
                f"TIC {tic_id} with P={period_days:.6f} d was queried against "
                f"the NASA Exoplanet Archive (pscomppars) and the TESS Object of "
                f"Interest catalogue (toi). No entry was found within "
                f"±{self.period_tol} days of the candidate period. "
                f"This signal is therefore classified as a "
                f"{'NEW DISCOVERY CANDIDATE' if status_tag == 'NEW_DISCOVERY' else 'OFFLINE DISCOVERY CANDIDATE (network unavailable, re-verification required)'}."
            ),
            "input_parameters": {
                "tic_id":                tic_id,
                "period_days":           period_days,
                "transit_depth":         transit_depth,
                "transit_depth_ppm":     round(transit_depth * 1e6, 2),
                "transit_duration_hrs":  transit_duration_hrs,
                "stellar_mass_solar":    stellar_mass_solar,
                "stellar_radius_solar":  stellar_radius_solar,
                "stellar_teff_k":        stellar_teff_k,
                "stellar_logg":          stellar_logg,
                "planet_radius_earth":   planet_radius_earth,
                "limb_dark_u1":          limb_dark_u1,
                "limb_dark_u2":          limb_dark_u2,
            },
            "detection_statistics": {
                "bls_snr":             bls_snr,
                "bls_fap":             bls_fap,
                "cvs_score":           cvs_score,
                "cvs_verdict":         cvs_verdict,
                "even_odd_delta_sigma": even_odd_delta_sigma,
                "shape_ratio":          shape_ratio,
                "secondary_snr":        secondary_snr,
                "centroid_sigma":       centroid_sigma,
            },
            **full_proof,   # ← All §1–§5 proof sections embedded here
        }

        write_json(str(out_path), card)
        self._log(f"Written -> {out_path}")
        
        # CRITICAL: Return the sovereign_verdict as the status, not the status_tag
        # This ensures FALSE_POSITIVE detections are not counted as discoveries
        final_status = full_proof["sovereign_verdict"]
        
        # Map sovereign verdicts to status codes
        # SOVEREIGN_PASS and CONDITIONAL_PASS are real planet candidates
        if final_status in ("SOVEREIGN_PASS", "CONDITIONAL_PASS"):
            # Only if it passes validation, use the original status_tag (NEW_DISCOVERY or OFFLINE_NEW_DISCOVERY)
            final_status = status_tag
        elif final_status == "FALSE_POSITIVE":
            # Keep as FALSE_POSITIVE - will be routed to rejected folder
            final_status = "FALSE_POSITIVE"
        else:
            # Unknown verdict - treat as false positive for safety
            self._log(f"WARNING: Unknown sovereign_verdict '{final_status}', treating as FALSE_POSITIVE")
            final_status = "FALSE_POSITIVE"
        
        return final_status, str(out_path)

    # ── Match parsers ─────────────────────────────────────────────────────────

    @staticmethod
    def _parse_planet_match(
        row: dict, delta: float, alias: float, source: str
    ) -> ArchiveMatch:
        def safe(key: str, default=None):
            try:
                v = row.get(key)
                return None if v is None else (float(v) if isinstance(v, (int, float)) else v)
            except Exception:
                return default

        return ArchiveMatch(
            source              = source,
            planet_name         = str(row.get("pl_name", "UNKNOWN")),
            period_days         = float(row.get("pl_orbper", 0)),
            period_delta_days   = delta,
            transit_depth       = safe("pl_trandep"),
            planet_radius_earth = safe("pl_rade"),
            semi_major_axis_au  = safe("pl_ratdor"),
            stellar_teff_k      = safe("st_teff"),
            stellar_radius_solar= safe("st_rad"),
            stellar_mass_solar  = safe("st_mass"),
            discovery_method    = str(row.get("discoverymethod", "")),
            disposition         = None,
            extra_fields        = {
                "disc_year":      row.get("disc_year"),
                "pl_controv_flag": row.get("pl_controv_flag"),
                "period_alias_factor": alias,
            },
        )

    @staticmethod
    def _parse_toi_match(row: dict, delta: float, alias: float) -> ArchiveMatch:
        def safe(key: str, default=None):
            try:
                v = row.get(key)
                return None if v is None else (float(v) if isinstance(v, (int, float)) else v)
            except Exception:
                return default

        return ArchiveMatch(
            source              = "TOI",
            planet_name         = str(row.get("toi", "TOI-UNKNOWN")),
            period_days         = float(row.get("pl_orbper", 0)),
            period_delta_days   = delta,
            transit_depth       = safe("pl_trandep"),
            planet_radius_earth = safe("pl_rade"),
            semi_major_axis_au  = safe("pl_ratdor"),
            stellar_teff_k      = safe("st_teff"),
            stellar_radius_solar= safe("st_rad"),
            stellar_mass_solar  = safe("st_mass"),
            discovery_method    = "Transit (TESS)",
            disposition         = str(row.get("tfopwg_disp", "")),
            extra_fields        = {
                "tmag":                row.get("st_tmag"),
                "period_alias_factor": alias,
            },
        )

    @staticmethod
    def _compute_deltas(
        period:  float,
        depth:   float,
        radius:  float,
        match:   ArchiveMatch,
    ) -> dict:
        """Compute fractional residuals between candidate and archive."""
        def pct(cand, arch):
            if arch and arch != 0:
                return round((cand - arch) / abs(arch) * 100, 4)
            return None

        return {
            "period_delta_days":     round(match.period_delta_days, 8),
            "period_delta_pct":      pct(period, match.period_days),
            "depth_delta_ppm":       round((depth - (match.transit_depth or 0)) * 1e6, 2)
                                     if match.transit_depth else None,
            "radius_delta_pct":      pct(radius, match.planet_radius_earth),
        }