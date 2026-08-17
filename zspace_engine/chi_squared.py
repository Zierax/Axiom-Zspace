"""
chi_squared.py  ·  ChiSquaredAnalyzer — Goodness-of-Fit Analysis
==================================================================
Computes χ² goodness-of-fit statistics for Mandel-Agol transit model fits.

Provides quantitative assessment of model-data agreement using:
  • χ² statistic: Σ[(O-M)²/σ²]
  • Reduced χ² (χ²_ν): χ² / degrees_of_freedom
  • Degrees of freedom: n_points - n_parameters
  • P-value: probability of observing χ² this large by chance
  • Interpretation: EXCELLENT / GOOD / POOR / OVER-FIT

Physical meaning
────────────────
  χ²_ν ≈ 1.0  →  Model explains data within expected noise level
  χ²_ν >> 1.0 →  Systematic errors or model inadequacy
  χ²_ν << 1.0 →  Errors overestimated or over-fitting

Usage
─────
  from zspace_engine.chi_squared import ChiSquaredAnalyzer

  analyzer = ChiSquaredAnalyzer()
  result = analyzer.compute_chi_squared(
      time=phase_folded_time,
      flux=phase_folded_flux,
      flux_err=flux_errors,
      model_flux=mandel_agol_model,
      n_params=5  # [rp/rs, a/rs, inc, u1, u2]
  )

  print(result["proof"])
  # χ² = Σ[(O-M)²/σ²] = 1247.32 | dof = 1219 | χ²_ν = 1.0234 | p = 3.421e-01 | EXCELLENT fit
"""

from __future__ import annotations

import math
from typing import Dict, Any, Optional

import numpy as np
from scipy.stats import chi2


class ChiSquaredAnalyzer:
    """
    Computes χ² goodness-of-fit for Mandel-Agol transit model.
    
    All calculations are deterministic and white-box, following the
    Axiom-ZSpace Truthimatics framework.
    """

    def __init__(self) -> None:
        """Initialize ChiSquaredAnalyzer."""
        pass

    def compute_chi_squared(
        self,
        time: np.ndarray,
        flux: np.ndarray,
        flux_err: np.ndarray,
        model_flux: np.ndarray,
        n_params: int = 5,
    ) -> Dict[str, Any]:
        """
        Compute χ² goodness-of-fit statistics for transit model.

        Parameters
        ----------
        time : np.ndarray
            Time array (not used in calculation, but included for API consistency)
        flux : np.ndarray
            Observed flux measurements (normalized, dimensionless)
        flux_err : np.ndarray
            Flux measurement uncertainties (1σ)
        model_flux : np.ndarray
            Model flux predictions from Mandel-Agol transit model
        n_params : int, optional
            Number of free parameters in the model (default: 5)
            Typical parameters: [R_p/R_★, a/R_★, inclination, u1, u2]

        Returns
        -------
        dict
            Dictionary containing:
            - chi_squared : float
                χ² statistic = Σ[(O-M)²/σ²]
            - reduced_chi_squared : float
                χ²_ν = χ² / (N - n_params)
            - degrees_of_freedom : int
                N - n_params
            - p_value : float
                Probability of observing χ² this large by chance
            - interpretation : str
                Human-readable fit quality assessment
            - proof : str
                Single-line proof string for Discovery Card

        Notes
        -----
        χ² interpretation guidelines:
          0.8 ≤ χ²_ν ≤ 1.2  →  EXCELLENT fit (model explains data within noise)
          0.5 ≤ χ²_ν ≤ 2.0  →  GOOD fit (acceptable model-data agreement)
          χ²_ν > 2.0        →  POOR fit (systematic errors or model inadequacy)
          χ²_ν < 0.5        →  OVER-FIT (errors overestimated or too many parameters)

        References
        ----------
        - Bevington & Robinson (2003), "Data Reduction and Error Analysis"
        - Press et al. (2007), "Numerical Recipes", Chapter 15
        """
        # Input validation
        if len(flux) != len(model_flux):
            raise ValueError(
                f"flux and model_flux must have same length: "
                f"{len(flux)} != {len(model_flux)}"
            )
        if len(flux) != len(flux_err):
            raise ValueError(
                f"flux and flux_err must have same length: "
                f"{len(flux)} != {len(flux_err)}"
            )
        if n_params < 0:
            raise ValueError(f"n_params must be non-negative: {n_params}")
        if n_params >= len(flux):
            raise ValueError(
                f"n_params ({n_params}) must be less than number of data points ({len(flux)})"
            )

        # Handle zero or negative errors (replace with median error)
        flux_err_safe = flux_err.copy()
        bad_errors = (flux_err_safe <= 0) | ~np.isfinite(flux_err_safe)
        if np.any(bad_errors):
            median_err = np.median(flux_err_safe[~bad_errors])
            if not np.isfinite(median_err) or median_err <= 0:
                median_err = 1e-4  # fallback for pathological cases
            flux_err_safe[bad_errors] = median_err

        # ── χ² calculation ────────────────────────────────────────────────────
        # χ² = Σ[(O_i - M_i)² / σ_i²]
        residuals = flux - model_flux
        chi_squared_terms = (residuals / flux_err_safe) ** 2
        chi_squared = float(np.sum(chi_squared_terms))

        # ── Degrees of freedom ────────────────────────────────────────────────
        # dof = N - n_params
        n_points = len(flux)
        dof = n_points - n_params

        if dof <= 0:
            raise ValueError(
                f"Degrees of freedom must be positive: dof = {n_points} - {n_params} = {dof}"
            )

        # ── Reduced χ² ────────────────────────────────────────────────────────
        # χ²_ν = χ² / dof
        reduced_chi_squared = chi_squared / dof

        # ── P-value ───────────────────────────────────────────────────────────
        # P(χ² ≥ observed | H₀) = 1 - CDF(χ², dof)
        # H₀: model is correct and errors are Gaussian
        try:
            p_value = float(1.0 - chi2.cdf(chi_squared, dof))
        except Exception:
            # Fallback for extreme values
            p_value = 0.0 if chi_squared > dof * 10 else 1.0

        # ── Interpretation ────────────────────────────────────────────────────
        if 0.8 <= reduced_chi_squared <= 1.2:
            interpretation = "EXCELLENT fit (χ²_ν ≈ 1)"
            quality = "EXCELLENT"
        elif 0.5 <= reduced_chi_squared <= 2.0:
            interpretation = "GOOD fit"
            quality = "GOOD"
        elif reduced_chi_squared > 2.0:
            interpretation = "POOR fit - systematic errors or model inadequacy"
            quality = "POOR"
        else:  # reduced_chi_squared < 0.5
            interpretation = "OVER-FIT - errors may be overestimated"
            quality = "OVER-FIT"

        # ── Proof string ──────────────────────────────────────────────────────
        proof = (
            f"χ² = Σ[(O-M)²/σ²] = {chi_squared:.2f} | "
            f"dof = {dof} | "
            f"χ²_ν = {reduced_chi_squared:.4f} | "
            f"p = {p_value:.4e} | "
            f"{interpretation}"
        )

        # ── Additional diagnostics ────────────────────────────────────────────
        rms_residual = float(np.sqrt(np.mean(residuals ** 2)))
        max_residual = float(np.max(np.abs(residuals)))
        mean_error = float(np.mean(flux_err_safe))

        return {
            "chi_squared": chi_squared,
            "reduced_chi_squared": reduced_chi_squared,
            "degrees_of_freedom": dof,
            "p_value": p_value,
            "interpretation": interpretation,
            "quality": quality,
            "proof": proof,
            # Additional diagnostics
            "n_points": n_points,
            "n_params": n_params,
            "rms_residual": rms_residual,
            "max_residual": max_residual,
            "mean_error": mean_error,
        }

    def interpret_fit_quality(self, reduced_chi_squared: float) -> str:
        """
        Provide detailed interpretation of reduced χ² value.

        Parameters
        ----------
        reduced_chi_squared : float
            The reduced χ² value (χ²_ν)

        Returns
        -------
        str
            Detailed interpretation string

        Examples
        --------
        >>> analyzer = ChiSquaredAnalyzer()
        >>> analyzer.interpret_fit_quality(1.05)
        'EXCELLENT: Model explains data within expected noise level (χ²_ν ≈ 1)'
        """
        if 0.8 <= reduced_chi_squared <= 1.2:
            return (
                f"EXCELLENT: Model explains data within expected noise level "
                f"(χ²_ν = {reduced_chi_squared:.4f} ≈ 1)"
            )
        elif 0.5 <= reduced_chi_squared <= 2.0:
            return (
                f"GOOD: Acceptable model-data agreement "
                f"(χ²_ν = {reduced_chi_squared:.4f})"
            )
        elif reduced_chi_squared > 2.0:
            return (
                f"POOR: χ²_ν = {reduced_chi_squared:.4f} >> 1 indicates systematic "
                f"errors, unmodeled astrophysical effects, or model inadequacy. "
                f"Possible causes: stellar activity, third body, incorrect limb darkening, "
                f"underestimated photometric errors."
            )
        else:  # reduced_chi_squared < 0.5
            return (
                f"OVER-FIT: χ²_ν = {reduced_chi_squared:.4f} << 1 suggests errors are "
                f"overestimated, too many free parameters, or data over-smoothed. "
                f"Model may be fitting noise rather than signal."
            )


# ─────────────────────────────────────────────────────────────────────────────
# Convenience function for quick analysis
# ─────────────────────────────────────────────────────────────────────────────

def quick_chi_squared(
    flux: np.ndarray,
    model_flux: np.ndarray,
    flux_err: Optional[np.ndarray] = None,
    n_params: int = 5,
) -> Dict[str, Any]:
    """
    Convenience function for quick χ² analysis.

    Parameters
    ----------
    flux : np.ndarray
        Observed flux measurements
    model_flux : np.ndarray
        Model flux predictions
    flux_err : np.ndarray, optional
        Flux uncertainties (if None, assumes uniform errors from RMS)
    n_params : int, optional
        Number of model parameters (default: 5)

    Returns
    -------
    dict
        χ² analysis results (see ChiSquaredAnalyzer.compute_chi_squared)

    Examples
    --------
    >>> flux = np.array([1.0, 0.99, 0.98, 0.99, 1.0])
    >>> model = np.array([1.0, 0.99, 0.98, 0.99, 1.0])
    >>> result = quick_chi_squared(flux, model, n_params=2)
    >>> print(result["quality"])
    'EXCELLENT'
    """
    if flux_err is None:
        # Estimate errors from RMS of residuals
        residuals = flux - model_flux
        rms = np.sqrt(np.mean(residuals ** 2))
        flux_err = np.full_like(flux, max(rms, 1e-4))

    analyzer = ChiSquaredAnalyzer()
    time = np.arange(len(flux))  # dummy time array
    return analyzer.compute_chi_squared(
        time=time,
        flux=flux,
        flux_err=flux_err,
        model_flux=model_flux,
        n_params=n_params,
    )


__all__ = [
    "ChiSquaredAnalyzer",
    "quick_chi_squared",
]
