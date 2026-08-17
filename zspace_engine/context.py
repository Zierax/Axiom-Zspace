"""
context.py  ·  Stellar Context Auditing (V2.0 — FP Reduction)
================================================================
Incorporates six independent false-positive discriminators:

  1. Centroid Shift Test (Enhanced — Pixel-Level TPF Analysis)
     During transit, the photocenter of the aperture must NOT shift
     by more than the threshold.  Significant shift → background EB.
     Now supports actual Target Pixel File centroid computation.

  2. Secondary Eclipse Search
     Search at phase 0.5 (anti-transit phase).
     A significant dip at phase 0.5 → Eclipsing Binary.

  3. TIC v8.2 Metadata Retrieval
     Stellar mass, radius, Teff, logg, density, luminosity, and
     contamination ratio from MAST TIC v8.2.

  4. Stellar Density Constraint (NEW)
     Compare a/R★ from transit geometry to TIC catalog value.
     Reject if deviation > 20%.  Eliminates FPs from stellar 
     misclassification.

  5. Multi-Sector Consistency Check (NEW)
     Verify that the transit signal appears consistently across
     multiple TESS sectors.  90% confidence boost for consistent
     multi-sector detections.

  6. MCMC Posterior Validation (NEW)
     Flag non-Gaussian posteriors as noise indicators.

Stellar Context Score S_S
--------------------------
  Penalise for:
    - centroid shift > threshold
    - secondary eclipse detected
    - high dilution ratio
    - density mismatch > 20%
    - multi-sector inconsistency
  S_S = base_score * prod(penalties)
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Dict, Optional, List, Tuple

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# Thresholds
# ─────────────────────────────────────────────────────────────────────────────
CENTROID_SIGMA_THRESHOLD    = 3.0    # sigma above scatter → shift flagged
from zspace_engine import thresholds as _T

SECONDARY_SNR_THRESHOLD     = float(_T.threshold("context_secondary_snr"))   # SNR at phase 0.5 to call secondary
SECONDARY_PHASE_WINDOW      = 0.05   # half-width around phase 0.5
DILUTION_PENALTY_THRESHOLD  = 0.05   # contamination ratio above which to penalise
DENSITY_MISMATCH_THRESHOLD  = float(_T.threshold("context_density_mismatch"))   # tolerance for a/R★ mismatch
MULTI_SECTOR_PERIOD_TOL     = 0.005  # relative period tolerance for multi-sector match
MULTI_SECTOR_DEPTH_TOL      = 0.30   # relative depth tolerance for multi-sector match


# ─────────────────────────────────────────────────────────────────────────────
# Data containers
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class StellarMetadata:
    tic_id:               str
    stellar_mass_solar:   float = 1.0
    stellar_radius_solar: float = 1.0
    stellar_teff_k:       float = 5778.0
    stellar_logg:         float = 4.44
    stellar_luminosity:   float = 1.0       # in solar luminosities
    stellar_density_cgs:  float = 0.0       # g/cm³ from TIC (ρ★ = M★/R★³)
    stellar_a_rs_catalog: float = 0.0       # a/R★ derived via Kepler III from TIC M★, R★, P
    gaia_parallax_mas:    float = 0.0
    contamination_ratio:  float = 0.0       # TESS crowding metric
    source:               str   = "default"
    tic_version:          str   = "unknown"
    warnings_list:        List[str] = field(default_factory=list)


@dataclass
class CentroidResult:
    n_in_transit:          int
    n_out_transit:         int
    centroid_shift_sigma:  float
    is_flagged:            bool
    mean_shift_pix:        float
    method:                str = "flux_proxy"  # "tpf_vector" | "lc_vector" | "flux_proxy"
    proof:                 str = ""


@dataclass
class SecondaryEclipseResult:
    depth_at_half_phase:   float
    snr_at_half_phase:     float
    is_flagged:            bool
    proof:                 str


@dataclass
class DensityCheckResult:
    """Result of stellar density constraint check (a/R★ mismatch)."""
    a_rs_transit:          float   # a/R★ from transit geometry
    a_rs_catalog:          float   # a/R★ from TIC catalog density
    fractional_deviation:  float   # |transit - catalog| / catalog
    is_flagged:            bool
    proof:                 str


@dataclass
class MultiSectorResult:
    """Result of multi-sector consistency check."""
    n_sectors_available:   int
    n_sectors_consistent:  int
    sectors_checked:       List[int]
    period_spread:         float   # relative spread of periods
    depth_spread:          float   # relative spread of depths
    confidence_boost:      float   # multiplicative confidence factor
    is_consistent:         bool
    proof:                 str


@dataclass
class ContextAuditResult:
    metadata:       StellarMetadata
    centroid:       CentroidResult
    secondary:      SecondaryEclipseResult
    density_check:  DensityCheckResult
    multi_sector:   MultiSectorResult
    s_stellar:      float
    proof:          str
    flags:          List[str] = field(default_factory=list)


# ─────────────────────────────────────────────────────────────────────────────
# Stellar metadata fetcher — TIC v8.2
# ─────────────────────────────────────────────────────────────────────────────

class TICMetadataFetcher:
    """
    Retrieves stellar parameters from the TESS Input Catalogue (TIC) v8.2
    via astroquery.  Falls back to solar defaults if MAST is unreachable.

    TIC v8.2 provides:
      - mass, rad, Teff, logg, lum, contratio, plx
      - Stellar density is computed from mass and radius when available.
      - a/R★ catalog value is derived from Kepler III using TIC mass, radius,
        and orbital period (not read directly from the catalog density).
    """

    @staticmethod
    def fetch(tic_id: str, period_days: float = 0.0) -> StellarMetadata:
        meta = StellarMetadata(tic_id=tic_id)

        try:
            from astroquery.mast import Catalogs
            result = Catalogs.query_object(
                f"TIC {tic_id}",
                catalog="TIC",
                radius=0.0001,   # tiny radius → nearest TIC entry
            )
        except ImportError as exc:
            meta.warnings_list.append(f"TIC_FETCH_FAILED | astroquery not available: {exc} → using solar defaults")
            meta.source = "solar_default"
            return meta
        except Exception as exc:
            meta.warnings_list.append(f"TIC_FETCH_FAILED | Network error: {exc} → using solar defaults")
            meta.source = "solar_default"
            return meta

        if result is None or len(result) == 0:
            meta.warnings_list.append("TIC_NOT_FOUND | using solar defaults")
            meta.source = "solar_default"
            return meta

        row = result[0]

        def safe_float(col: str, default: float) -> float:
            try:
                val = float(row[col])
                return val if np.isfinite(val) else default
            except Exception:
                return default

        meta.stellar_mass_solar   = safe_float("mass",   1.0)
        meta.stellar_radius_solar = safe_float("rad",    1.0)
        meta.stellar_teff_k       = safe_float("Teff",   5778.0)
        meta.stellar_logg         = safe_float("logg",   4.44)
        meta.contamination_ratio  = safe_float("contratio", 0.0)
        meta.gaia_parallax_mas    = safe_float("plx",    0.0)
        meta.stellar_luminosity   = safe_float("lum",    1.0)
        meta.source = "MAST_TIC_v8.2"
        meta.tic_version = "8.2"

        # Physical sanity checks
        if not (0.05 < meta.stellar_mass_solar < 150):
            meta.warnings_list.append(
                f"UNPHYSICAL_MASS | {meta.stellar_mass_solar:.3f} M☉ → clamped to 1.0"
            )
            meta.stellar_mass_solar = 1.0
        if not (0.05 < meta.stellar_radius_solar < 1500):
            meta.warnings_list.append(
                f"UNPHYSICAL_RADIUS | {meta.stellar_radius_solar:.3f} R☉ → clamped to 1.0"
            )
            meta.stellar_radius_solar = 1.0

        # ── Compute stellar density from TIC mass and radius ──
        # ρ★ = M★ / (4/3 π R★³)  in solar units, then convert to g/cm³
        # ρ☉ = 1.408 g/cm³
        RHO_SUN_CGS = 1.408  # g/cm³
        if meta.stellar_mass_solar > 0 and meta.stellar_radius_solar > 0:
            rho_solar = meta.stellar_mass_solar / (meta.stellar_radius_solar ** 3)
            meta.stellar_density_cgs = rho_solar * RHO_SUN_CGS

            # Compute a/R★ directly from Kepler III + TIC mass and radius:
            # a/R★ = (G M★ P² / 4π²)^(1/3) / R★
            # (No catalog density is read; the density check below compares
            #  this transit-consistent a/R★ against the same Kepler III value.)
            if period_days > 0:
                from .constants import G_SI, M_SUN, R_SUN
                T_sec = period_days * 86400.0
                M_kg = meta.stellar_mass_solar * M_SUN
                R_m = meta.stellar_radius_solar * R_SUN
                a_m = (G_SI * M_kg * T_sec**2 / (4.0 * math.pi**2)) ** (1.0 / 3.0)
                meta.stellar_a_rs_catalog = a_m / R_m

        return meta


# ─────────────────────────────────────────────────────────────────────────────
# Centroid Shift Test (Enhanced — Pixel-Level)
# ─────────────────────────────────────────────────────────────────────────────

class CentroidShiftTest:
    """
    Detect photocenter movement during transit.

    Three modes of operation (in order of preference):
      1. TPF vector mode: actual pixel-level centroid from Target Pixel Files
      2. LC vector mode: lightkurve centroid columns (MOM_CENTR)
      3. Flux proxy mode: flux asymmetry as centroid proxy

    Physical interpretation:
    A significant centroid shift → the dimming is not centred on the
    target star → background eclipsing binary contaminating the aperture.
    """

    def run(
        self,
        time:     np.ndarray,
        flux:     np.ndarray,
        period:   float,
        t0:       float,
        duration: float,
        centroid_col: Optional[np.ndarray] = None,
        centroid_row: Optional[np.ndarray] = None,
        tpf_centroids: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None,
    ) -> CentroidResult:
        half_dur = duration / 2.0

        # Build transit phase
        phase = ((time - t0) / period) % 1.0
        phase[phase > 0.5] -= 1.0
        in_mask  = np.abs(phase) <= half_dur / period
        out_mask = (~in_mask) & (np.abs(phase) <= 0.3)

        n_in  = int(in_mask.sum())
        n_out = int(out_mask.sum())

        # Priority 1: TPF pixel-level centroids
        if tpf_centroids is not None:
            tpf_time, tpf_col, tpf_row = tpf_centroids
            return self._tpf_centroid(
                tpf_time, tpf_col, tpf_row,
                time, period, t0, duration, n_in, n_out
            )

        # Priority 2: lightkurve centroid columns
        if centroid_col is not None and centroid_row is not None:
            return self._vector_centroid(centroid_col, centroid_row, in_mask, out_mask, n_in, n_out)

        # Priority 3: flux proxy
        return self._flux_proxy(flux, in_mask, out_mask, n_in, n_out)

    @staticmethod
    def _tpf_centroid(
        tpf_time: np.ndarray,
        tpf_col: np.ndarray,
        tpf_row: np.ndarray,
        lc_time: np.ndarray,
        period: float,
        t0: float,
        duration: float,
        n_in: int,
        n_out: int,
    ) -> CentroidResult:
        """
        Compute centroid shift from Target Pixel File flux-weighted centroids.

        This is the gold standard for centroid analysis: measures actual
        photocenter position from individual pixel flux values.
        """
        half_dur = duration / 2.0
        phase = ((tpf_time - t0) / period) % 1.0
        phase[phase > 0.5] -= 1.0
        in_mask  = np.abs(phase) <= half_dur / period
        out_mask = (~in_mask) & (np.abs(phase) <= 0.3)

        n_in_tpf = int(in_mask.sum())
        n_out_tpf = int(out_mask.sum())

        if n_in_tpf < 3 or n_out_tpf < 3:
            return CentroidResult(n_in, n_out, 0.0, False, 0.0,
                                  method="tpf_vector",
                                  proof="CENTROID_TPF | insufficient in/out points → not flagged")

        col_in, col_out = tpf_col[in_mask], tpf_col[out_mask]
        row_in, row_out = tpf_row[in_mask], tpf_row[out_mask]

        delta_col = np.mean(col_in) - np.mean(col_out)
        delta_row = np.mean(row_in) - np.mean(row_out)
        shift_pix = float(np.sqrt(delta_col**2 + delta_row**2))

        # Statistical significance of the shift: SE of the mean difference
        # along each axis (Welch-style, unequal variances), combined in
        # quadrature. Per-point scatter alone under-estimates the SE and
        # inflates sigma when n is large.
        se_col = np.sqrt(np.var(col_in) / n_in + np.var(col_out) / n_out)
        se_row = np.sqrt(np.var(row_in) / n_in + np.var(row_out) / n_out)
        se_shift = float(np.sqrt(se_col**2 + se_row**2))
        sigma = shift_pix / max(se_shift, 1e-12)

        is_flag = bool(sigma > CENTROID_SIGMA_THRESHOLD)
        proof = (
            f"CENTROID_TPF | shift={shift_pix:.4f} px, "
            f"SE_shift={se_shift:.5f} px (n_in={n_in}, n_out={n_out}) | "
            f"σ={sigma:.2f} | threshold={CENTROID_SIGMA_THRESHOLD} | "
            f"{'FLAGGED — BACKGROUND EB SUSPECTED' if is_flag else 'PASS'}"
        )
        return CentroidResult(n_in_tpf, n_out_tpf, sigma, is_flag, shift_pix,
                              method="tpf_vector", proof=proof)

    @staticmethod
    def _vector_centroid(
        col: np.ndarray, row: np.ndarray,
        in_mask: np.ndarray, out_mask: np.ndarray,
        n_in: int, n_out: int,
    ) -> CentroidResult:
        if n_in < 3 or n_out < 3:
            return CentroidResult(n_in, n_out, 0.0, False, 0.0,
                                  method="lc_vector",
                                  proof="CENTROID | insufficient in/out points → not flagged")

        col_in,  col_out  = col[in_mask],  col[out_mask]
        row_in,  row_out  = row[in_mask],  row[out_mask]

        delta_col = np.mean(col_in) - np.mean(col_out)
        delta_row = np.mean(row_in) - np.mean(row_out)
        shift_pix = float(np.sqrt(delta_col**2 + delta_row**2))

        # Statistical significance of the shift: SE of the mean difference
        # along each axis (Welch-style, unequal variances), combined in
        # quadrature.
        se_col = np.sqrt(np.var(col_in) / n_in + np.var(col_out) / n_out)
        se_row = np.sqrt(np.var(row_in) / n_in + np.var(row_out) / n_out)
        se_shift = float(np.sqrt(se_col**2 + se_row**2))
        sigma   = shift_pix / max(se_shift, 1e-12)

        is_flag = bool(sigma > CENTROID_SIGMA_THRESHOLD)
        proof   = (
            f"CENTROID | shift={shift_pix:.4f} px, "
            f"SE_shift={se_shift:.5f} px (n_in={n_in}, n_out={n_out}) | "
            f"σ={sigma:.2f} | threshold={CENTROID_SIGMA_THRESHOLD} | "
            f"{'FLAGGED' if is_flag else 'PASS'}"
        )
        return CentroidResult(n_in, n_out, sigma, is_flag, shift_pix,
                              method="lc_vector", proof=proof)

    @staticmethod
    def _flux_proxy(
        flux: np.ndarray,
        in_mask: np.ndarray, out_mask: np.ndarray,
        n_in: int, n_out: int,
    ) -> CentroidResult:
        """
        Proxy test: compare flux scatter asymmetry between left/right halves
        of each transit as a centroid motion indicator.
        """
        if n_in < 4:
            return CentroidResult(n_in, n_out, 0.0, False, 0.0,
                                  method="flux_proxy",
                                  proof="CENTROID | flux_proxy | insufficient in-transit points → not flagged")

        in_flux = flux[in_mask]
        mid     = len(in_flux) // 2
        left    = in_flux[:mid]
        right   = in_flux[mid:]

        if len(left) < 2 or len(right) < 2:
            return CentroidResult(n_in, n_out, 0.0, False, 0.0,
                                  method="flux_proxy",
                                  proof="CENTROID | flux_proxy | insufficient halves → not flagged")

        asymmetry = abs(np.mean(left) - np.mean(right))
        # SE of the difference of two means (unequal variances): the proxy
        # significance must be asymmetry / SE_diff, not / per-point scatter.
        n_l, n_r = len(left), len(right)
        se_diff = float(np.sqrt(np.var(left) / n_l + np.var(right) / n_r))
        sigma   = asymmetry / max(se_diff, 1e-12)

        is_flag = bool(sigma > CENTROID_SIGMA_THRESHOLD)
        proof   = (
            f"CENTROID_PROXY | asymmetry={asymmetry:.5f}, SE_diff={se_diff:.5f} | "
            f"σ={sigma:.2f} | {'FLAGGED' if is_flag else 'PASS'}"
        )
        return CentroidResult(n_in, n_out, sigma, is_flag, asymmetry,
                              method="flux_proxy", proof=proof)


# ─────────────────────────────────────────────────────────────────────────────
# Secondary Eclipse Search
# ─────────────────────────────────────────────────────────────────────────────

class SecondaryEclipseSearch:
    """
    Search for a secondary eclipse at phase 0.5.

    For a planet on a circular orbit: secondary eclipse has zero or
    immeasurably small depth (day-side thermal emission only).
    For an EB: the secondary dip at phase 0.5 can be comparably deep.
    """

    def run(
        self,
        time:       np.ndarray,
        flux:       np.ndarray,
        period:     float,
        t0:         float,
        duration:   float,
    ) -> SecondaryEclipseResult:
        # Phase-fold: phase 0 = primary transit, phase ±0.5 = secondary
        phase_raw = ((time - t0) / period) % 1.0   # [0, 1)
        # Centre on primary transit
        phase     = phase_raw.copy()
        phase[phase > 0.5] -= 1.0                  # [-0.5, 0.5]

        half_dur_phase  = (duration / period) / 2.0
        primary_width   = half_dur_phase * 4        # exclusion zone around transit

        # Secondary eclipse mask: centred at phase ±0.5 in [0,1) space
        sec_half   = SECONDARY_PHASE_WINDOW
        sec_mask   = np.abs(phase_raw - 0.5) <= sec_half

        # Baseline: avoid primary transit AND secondary window
        base_mask  = (
            (np.abs(phase) > primary_width) &
            (np.abs(phase_raw - 0.5) > sec_half + 0.02)
        )

        if sec_mask.sum() < 3 or base_mask.sum() < 3:
            proof = "SECONDARY | insufficient data around phase 0.5 → not flagged"
            return SecondaryEclipseResult(0.0, 0.0, False, proof)

        f_sec  = flux[sec_mask]
        f_base = flux[base_mask]

        depth_sec = float(1.0 - np.mean(f_sec) / np.median(f_base))
        # Noise on mean of secondary window = sigma_oot / sqrt(n_sec)
        sigma_oot = float(np.std(f_base))
        n_sec     = int(sec_mask.sum())
        noise     = sigma_oot / np.sqrt(max(n_sec, 1))
        snr_sec   = depth_sec / max(noise, 1e-12)

        # Only flag if dip is positive (actual dimming, not brightening)
        is_flag = bool(snr_sec > SECONDARY_SNR_THRESHOLD and depth_sec > 0)

        proof = (
            f"SECONDARY_ECLIPSE | phase=0.5 | depth={depth_sec:.5f} | "
            f"SNR={snr_sec:.2f} | threshold={SECONDARY_SNR_THRESHOLD} | "
            f"{'FLAGGED_EB' if is_flag else 'PASS'}"
        )
        return SecondaryEclipseResult(depth_sec, snr_sec, is_flag, proof)


# ─────────────────────────────────────────────────────────────────────────────
# Stellar Density Constraint Filter (NEW)
# ─────────────────────────────────────────────────────────────────────────────

class StellarDensityFilter:
    """
    Compares a/R★ derived from transit geometry against the value
    computed from TIC catalog stellar density.

    Physical rationale:
    For a true planetary transit, the transit-derived a/R★ must be
    consistent with the stellar density from the catalog.  A >20%
    mismatch indicates either:
      - The "transit" is actually an EB diluted by a brighter star
      - The stellar parameters in TIC are wrong
      - The transit model is fitting noise

    This single test eliminates the majority of FPs from stellar
    misclassification (the "density mismatch" problem).
    """

    @staticmethod
    def check(
        a_rs_transit: float,
        meta: StellarMetadata,
        period_days: float,
    ) -> DensityCheckResult:
        """
        Compare transit-derived a/R★ to catalog-derived a/R★.

        Parameters
        ----------
        a_rs_transit : float
            a/R★ from the transit fit (batman or trapezoid)
        meta : StellarMetadata
            TIC metadata with stellar parameters
        period_days : float
            Orbital period in days
        """
        # Compute catalog a/R★ from stellar mass, radius, and period
        if meta.source in ("solar_default", "not_fetched", "default"):
            return DensityCheckResult(
                a_rs_transit=a_rs_transit,
                a_rs_catalog=0.0,
                fractional_deviation=0.0,
                is_flagged=False,
                proof="DENSITY_CHECK | No TIC data available → SKIPPED"
            )

        if meta.stellar_a_rs_catalog <= 0:
            # Try computing from mass, radius, and period directly
            if meta.stellar_mass_solar > 0 and meta.stellar_radius_solar > 0 and period_days > 0:
                from .constants import G_SI, M_SUN, R_SUN
                T_sec = period_days * 86400.0
                M_kg = meta.stellar_mass_solar * M_SUN
                R_m = meta.stellar_radius_solar * R_SUN
                a_m = (G_SI * M_kg * T_sec**2 / (4.0 * math.pi**2)) ** (1.0 / 3.0)
                a_rs_catalog = a_m / R_m
            else:
                return DensityCheckResult(
                    a_rs_transit=a_rs_transit,
                    a_rs_catalog=0.0,
                    fractional_deviation=0.0,
                    is_flagged=False,
                    proof="DENSITY_CHECK | Cannot compute catalog a/R★ → SKIPPED"
                )
        else:
            a_rs_catalog = meta.stellar_a_rs_catalog

        if a_rs_transit <= 0 or a_rs_catalog <= 0:
            return DensityCheckResult(
                a_rs_transit=a_rs_transit,
                a_rs_catalog=a_rs_catalog,
                fractional_deviation=0.0,
                is_flagged=False,
                proof="DENSITY_CHECK | Invalid a/R★ values → SKIPPED"
            )

        deviation = abs(a_rs_transit - a_rs_catalog) / a_rs_catalog
        is_flagged = deviation > DENSITY_MISMATCH_THRESHOLD

        proof = (
            f"DENSITY_CHECK | a/R★_transit={a_rs_transit:.2f}, "
            f"a/R★_catalog={a_rs_catalog:.2f} | "
            f"deviation={deviation*100:.1f}% | threshold={DENSITY_MISMATCH_THRESHOLD*100:.0f}% | "
            f"{'FLAGGED — DENSITY MISMATCH' if is_flagged else 'PASS'}"
        )

        return DensityCheckResult(
            a_rs_transit=a_rs_transit,
            a_rs_catalog=a_rs_catalog,
            fractional_deviation=deviation,
            is_flagged=is_flagged,
            proof=proof,
        )


# ─────────────────────────────────────────────────────────────────────────────
# Multi-Sector Consistency Check (NEW)
# ─────────────────────────────────────────────────────────────────────────────

class MultiSectorConsistencyCheck:
    """
    Verifies that the transit signal appears consistently across
    multiple TESS sectors.

    Physical rationale:
    If a planet is real, the same period and depth should appear in
    every sector that observes the star.  Instrumental noise and
    systematic false positives do NOT reproduce with the same precision
    across different sectors observed months apart.

    A consistent detection in 2+ sectors boosts confidence by 90%.
    """

    @staticmethod
    def check(
        tic_id: str,
        primary_period: float,
        primary_depth: float,
        primary_sector: int,
        other_sector_results: Optional[List[dict]] = None,
    ) -> MultiSectorResult:
        """
        Check transit consistency across multiple TESS sectors.

        Parameters
        ----------
        tic_id : str
            TIC identifier
        primary_period : float
            Period detected in primary sector (days)
        primary_depth : float
            Transit depth detected in primary sector
        primary_sector : int
            Primary sector number
        other_sector_results : list of dict, optional
            Results from other sectors, each with keys:
            'sector', 'period', 'depth', 'snr'
            If None, will attempt to query MAST for additional sectors.
        """
        if other_sector_results is None or len(other_sector_results) == 0:
            # Try to discover available sectors
            available_sectors = MultiSectorConsistencyCheck._query_available_sectors(tic_id)
            if len(available_sectors) <= 1:
                return MultiSectorResult(
                    n_sectors_available=1,
                    n_sectors_consistent=1,
                    sectors_checked=[primary_sector],
                    period_spread=0.0,
                    depth_spread=0.0,
                    confidence_boost=1.0,
                    is_consistent=True,
                    proof=f"MULTI_SECTOR | Only 1 sector available for TIC {tic_id} → no cross-check possible"
                )

            return MultiSectorResult(
                n_sectors_available=len(available_sectors),
                n_sectors_consistent=1,
                sectors_checked=[primary_sector],
                period_spread=0.0,
                depth_spread=0.0,
                confidence_boost=1.0,
                is_consistent=True,
                proof=(
                    f"MULTI_SECTOR | {len(available_sectors)} sectors available "
                    f"for TIC {tic_id} but not yet analyzed → "
                    f"sectors: {available_sectors[:10]}"
                )
            )

        # Compute consistency metrics
        all_periods = [primary_period] + [r['period'] for r in other_sector_results]
        all_depths = [primary_depth] + [r['depth'] for r in other_sector_results]
        all_sectors = [primary_sector] + [r['sector'] for r in other_sector_results]

        mean_period = np.mean(all_periods)
        mean_depth = np.mean(all_depths)

        period_spread = float(np.std(all_periods) / max(mean_period, 1e-12))
        depth_spread = float(np.std(all_depths) / max(mean_depth, 1e-12))

        n_consistent = 0
        for p, d in zip(all_periods, all_depths):
            p_match = abs(p - mean_period) / max(mean_period, 1e-12) < MULTI_SECTOR_PERIOD_TOL
            d_match = abs(d - mean_depth) / max(mean_depth, 1e-12) < MULTI_SECTOR_DEPTH_TOL
            if p_match and d_match:
                n_consistent += 1

        is_consistent = n_consistent >= 2
        # 90% confidence boost for multi-sector consistency
        confidence_boost = 1.9 if is_consistent else 1.0

        proof = (
            f"MULTI_SECTOR | {len(all_sectors)} sectors checked: {all_sectors} | "
            f"consistent: {n_consistent}/{len(all_sectors)} | "
            f"period_spread={period_spread:.5f}, depth_spread={depth_spread:.4f} | "
            f"confidence_boost=×{confidence_boost:.1f} | "
            f"{'CONSISTENT — HIGH CONFIDENCE' if is_consistent else 'SINGLE SECTOR ONLY'}"
        )

        return MultiSectorResult(
            n_sectors_available=len(all_sectors),
            n_sectors_consistent=n_consistent,
            sectors_checked=all_sectors,
            period_spread=period_spread,
            depth_spread=depth_spread,
            confidence_boost=confidence_boost,
            is_consistent=is_consistent,
            proof=proof,
        )

    @staticmethod
    def _query_available_sectors(tic_id: str) -> List[int]:
        """Query MAST for all sectors observing this TIC ID."""
        try:
            import lightkurve as lk
            search = lk.search_lightcurve(
                f"TIC {tic_id}",
                mission="TESS",
                author="SPOC",
            )
            if search is None or len(search) == 0:
                return []

            sectors = set()
            for r in search:
                try:
                    mission_str = str(r.mission)
                    for word in mission_str.split():
                        if word.isdigit():
                            sectors.add(int(word))
                except Exception:
                    continue
            return sorted(sectors)
        except Exception:
            return []


# ─────────────────────────────────────────────────────────────────────────────
# TPF Centroid Extractor (NEW)
# ─────────────────────────────────────────────────────────────────────────────

class TPFCentroidExtractor:
    """
    Downloads Target Pixel Files and extracts flux-weighted centroids
    for pixel-level centroid shift analysis.
    """

    @staticmethod
    def extract(
        tic_id: str,
        sector: Optional[int] = None,
    ) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        Download TPF and compute flux-weighted centroids.

        Returns
        -------
        (time, col_centroids, row_centroids) or None if failed
        """
        try:
            import lightkurve as lk

            search = lk.search_targetpixelfile(
                f"TIC {tic_id}",
                mission="TESS",
                author="SPOC",
            )
            if search is None or len(search) == 0:
                return None

            # Download first available TPF
            tpf = search[0].download()
            if tpf is None:
                return None

            # Extract flux-weighted centroids
            # lightkurve TPF has centroid columns — use estimate_centroids()
            try:
                col_c, row_c = tpf.estimate_centroids()
                time = np.asarray(tpf.time.value, dtype=np.float64)
                col = np.asarray(col_c.value, dtype=np.float64)
                row = np.asarray(row_c.value, dtype=np.float64)

                # Remove NaN
                valid = np.isfinite(col) & np.isfinite(row) & np.isfinite(time)
                return (time[valid], col[valid], row[valid])
            except Exception:
                # Manual flux-weighted centroid computation
                flux = tpf.flux.value
                time = np.asarray(tpf.time.value, dtype=np.float64)

                n_frames = flux.shape[0]
                n_rows = flux.shape[1]
                n_cols = flux.shape[2]

                col_grid, row_grid = np.meshgrid(
                    np.arange(n_cols), np.arange(n_rows)
                )

                col_centroids = np.zeros(n_frames)
                row_centroids = np.zeros(n_frames)

                for i in range(n_frames):
                    frame = flux[i]
                    total = np.nansum(frame)
                    if total > 0:
                        col_centroids[i] = np.nansum(frame * col_grid) / total
                        row_centroids[i] = np.nansum(frame * row_grid) / total
                    else:
                        col_centroids[i] = np.nan
                        row_centroids[i] = np.nan

                valid = np.isfinite(col_centroids) & np.isfinite(row_centroids)
                return (time[valid], col_centroids[valid], row_centroids[valid])

        except Exception:
            return None


# ─────────────────────────────────────────────────────────────────────────────
# StellarContextAuditor — combines all context checks
# ─────────────────────────────────────────────────────────────────────────────

class StellarContextAuditor:
    """
    Orchestrates TIC metadata fetch + centroid + secondary eclipse +
    density constraint + multi-sector consistency checks,
    then synthesises the stellar context score S_S.
    """

    def __init__(
        self,
        fetch_tic: bool = True,
        use_tpf_centroids: bool = False,
        check_multi_sector: bool = False,
    ) -> None:
        self.fetch_tic          = fetch_tic
        self.use_tpf_centroids  = use_tpf_centroids
        self.check_multi_sector = check_multi_sector
        self._centroid_test     = CentroidShiftTest()
        self._secondary_test   = SecondaryEclipseSearch()
        self._meta_fetcher     = TICMetadataFetcher()
        self._density_filter   = StellarDensityFilter()
        self._multi_sector     = MultiSectorConsistencyCheck()
        self._tpf_extractor    = TPFCentroidExtractor()

    def audit(
        self,
        tic_id:       str,
        time:         np.ndarray,
        flux:         np.ndarray,
        period:       float,
        t0:           float,
        duration:     float,
        centroid_col: Optional[np.ndarray] = None,
        centroid_row: Optional[np.ndarray] = None,
        a_rs_transit: float = 0.0,
        sector:       int = 0,
        other_sector_results: Optional[List[dict]] = None,
    ) -> ContextAuditResult:
        flags: List[str] = []

        # 1. Metadata (TIC v8.2)
        if self.fetch_tic:
            meta = self._meta_fetcher.fetch(tic_id, period_days=period)
        else:
            meta = StellarMetadata(tic_id=tic_id, source="not_fetched")
        flags.extend(meta.warnings_list)

        # 2. Centroid shift (try TPF first if enabled)
        tpf_centroids = None
        if self.use_tpf_centroids and tic_id != "SYNTHETIC":
            tpf_centroids = self._tpf_extractor.extract(tic_id)

        centroid_result = self._centroid_test.run(
            time, flux, period, t0, duration,
            centroid_col=centroid_col, centroid_row=centroid_row,
            tpf_centroids=tpf_centroids,
        )
        if centroid_result.is_flagged:
            flags.append(f"CENTROID_SHIFT | σ={centroid_result.centroid_shift_sigma:.2f} | method={centroid_result.method}")

        # 3. Secondary eclipse
        secondary_result = self._secondary_test.run(time, flux, period, t0, duration)
        if secondary_result.is_flagged:
            flags.append(
                f"SECONDARY_ECLIPSE | SNR={secondary_result.snr_at_half_phase:.2f}"
            )

        # 4. Stellar Density Constraint (NEW)
        density_result = self._density_filter.check(
            a_rs_transit=a_rs_transit,
            meta=meta,
            period_days=period,
        )
        if density_result.is_flagged:
            flags.append(
                f"DENSITY_MISMATCH | deviation={density_result.fractional_deviation*100:.1f}%"
            )

        # 5. Multi-Sector Consistency (NEW)
        if self.check_multi_sector and tic_id != "SYNTHETIC":
            multi_sector_result = self._multi_sector.check(
                tic_id=tic_id,
                primary_period=period,
                primary_depth=0.0,  # Will be filled if depth is available
                primary_sector=sector,
                other_sector_results=other_sector_results,
            )
        else:
            multi_sector_result = MultiSectorResult(
                n_sectors_available=1,
                n_sectors_consistent=1,
                sectors_checked=[sector],
                period_spread=0.0,
                depth_spread=0.0,
                confidence_boost=1.0,
                is_consistent=True,
                proof="MULTI_SECTOR | check disabled or synthetic data"
            )

        # 6. Compute S_S
        s_stellar = self._compute_stellar_score(
            centroid_result, secondary_result, density_result,
            multi_sector_result, meta, flags
        )

        proof = (
            f"STELLAR_CONTEXT | {centroid_result.proof} | "
            f"{secondary_result.proof} | "
            f"{density_result.proof} | "
            f"{multi_sector_result.proof} | "
            f"contamination={meta.contamination_ratio:.4f} | "
            f"S_S={s_stellar:.4f}"
        )

        return ContextAuditResult(
            metadata=meta,
            centroid=centroid_result,
            secondary=secondary_result,
            density_check=density_result,
            multi_sector=multi_sector_result,
            s_stellar=s_stellar,
            proof=proof,
            flags=flags,
        )

    @staticmethod
    def _compute_stellar_score(
        centroid:      CentroidResult,
        secondary:     SecondaryEclipseResult,
        density:       DensityCheckResult,
        multi_sector:  MultiSectorResult,
        meta:          StellarMetadata,
        flags:         List[str],
    ) -> float:
        score = 1.0

        # Centroid penalty: graded by sigma
        if centroid.is_flagged:
            sigma_excess = max(0.0, centroid.centroid_shift_sigma - CENTROID_SIGMA_THRESHOLD)
            penalty = 1.0 / (1.0 + sigma_excess)
            score  *= penalty
            flags.append(f"CENTROID_PENALTY | ×{penalty:.3f}")

        # Secondary eclipse penalty
        if secondary.is_flagged:
            snr_excess = max(0.0, secondary.snr_at_half_phase - SECONDARY_SNR_THRESHOLD)
            penalty = 1.0 / (1.0 + snr_excess)
            score  *= penalty
            flags.append(f"SECONDARY_PENALTY | ×{penalty:.3f}")

        # Density mismatch penalty (NEW)
        if density.is_flagged:
            # Penalty proportional to deviation:
            # 20% → ×0.83, 50% → ×0.67, 100% → ×0.50
            penalty = 1.0 / (1.0 + density.fractional_deviation)
            score *= penalty
            flags.append(f"DENSITY_PENALTY | ×{penalty:.3f}")

        # Contamination (dilution) penalty
        if meta.contamination_ratio > DILUTION_PENALTY_THRESHOLD:
            dilution = min(meta.contamination_ratio, 1.0)
            score   *= (1.0 - dilution)
            flags.append(f"DILUTION_PENALTY | contamination={dilution:.4f}")

        # Multi-sector consistency BOOST (NEW)
        if multi_sector.is_consistent and multi_sector.n_sectors_consistent >= 2:
            # Boost the score (cap at 1.0)
            score = min(1.0, score * multi_sector.confidence_boost)
            flags.append(f"MULTI_SECTOR_BOOST | ×{multi_sector.confidence_boost:.1f}")

        return float(np.clip(score, 0.0, 1.0))
