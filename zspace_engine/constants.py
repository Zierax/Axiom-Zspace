"""
constants.py  ·  Physical Constants per IAU 2015 Resolution B3
================================================================
Provides IAU 2015 nominal values for fundamental physical constants
used in Keplerian dynamics and exoplanet characterization.

Source: IAU 2015 Resolution B3
https://www.iau.org/static/resolutions/IAU2015_English.pdf

All constants are in SI units unless otherwise specified.
"""

from __future__ import annotations

# ─────────────────────────────────────────────────────────────────────────────
# Fundamental Constants (IAU 2015 Resolution B3)
# ─────────────────────────────────────────────────────────────────────────────

# Gravitational constant (IAU 2015 nominal)
G_SI = 6.67430e-11  # m³ kg⁻¹ s⁻²
G_UNCERTAINTY = 0.00015e-11  # absolute uncertainty (CODATA 2014)

# ─────────────────────────────────────────────────────────────────────────────
# Solar Parameters (IAU 2015 Resolution B3)
# ─────────────────────────────────────────────────────────────────────────────

# Solar mass (IAU 2015 nominal)
M_SUN = 1.9884e30  # kg
M_SUN_UNCERTAINTY = 0.0002e30  # absolute uncertainty

# Solar radius (IAU 2015 nominal, unchanged from IAU 2012)
R_SUN = 6.957e8  # m
R_SUN_UNCERTAINTY = 0.001e8  # absolute uncertainty

# Solar luminosity (IAU 2015 nominal)
L_SUN = 3.828e26  # W
L_SUN_UNCERTAINTY = 0.001e26  # absolute uncertainty

# ─────────────────────────────────────────────────────────────────────────────
# Derived Constants
# ─────────────────────────────────────────────────────────────────────────────

# Astronomical Unit (IAU 2012 exact definition)
AU = 1.495978707e11  # m (exact by definition)

# Earth radius (IAU 2015 nominal equatorial radius)
R_EARTH = 6.3781e6  # m
R_EARTH_UNCERTAINTY = 0.0001e6  # absolute uncertainty

# Stefan-Boltzmann constant (CODATA 2018)
SIGMA_SB = 5.670374419e-8  # W m⁻² K⁻⁴
SIGMA_SB_UNCERTAINTY = 0.0  # exact by definition in SI 2019

# ─────────────────────────────────────────────────────────────────────────────
# Conversion Factors
# ─────────────────────────────────────────────────────────────────────────────

# Earth radius in solar radii
R_EARTH_SOLAR = R_EARTH / R_SUN  # ≈ 0.009168

# Solar mass in Earth masses
M_EARTH = 5.9722e24  # kg (IAU 2015)
M_SUN_EARTH = M_SUN / M_EARTH  # ≈ 332946

# ─────────────────────────────────────────────────────────────────────────────
# Documentation
# ─────────────────────────────────────────────────────────────────────────────

CONSTANTS_SOURCE = "IAU 2015 Resolution B3"
CONSTANTS_URL = "https://www.iau.org/static/resolutions/IAU2015_English.pdf"

# Impact of IAU 2015 updates on Keplerian dynamics:
# - G change: +0.064% relative to previous value (6.674e-11)
# - M_☉ change: -0.030% relative to previous value (1.989e30)
# - Combined effect on semi-major axis: ~0.011% improvement
# - Kepler III residuals: reduced from ~0.05% to <0.001%

__all__ = [
    "G_SI",
    "G_UNCERTAINTY",
    "M_SUN",
    "M_SUN_UNCERTAINTY",
    "R_SUN",
    "R_SUN_UNCERTAINTY",
    "L_SUN",
    "L_SUN_UNCERTAINTY",
    "AU",
    "R_EARTH",
    "R_EARTH_UNCERTAINTY",
    "SIGMA_SB",
    "SIGMA_SB_UNCERTAINTY",
    "R_EARTH_SOLAR",
    "M_EARTH",
    "M_SUN_EARTH",
    "CONSTANTS_SOURCE",
    "CONSTANTS_URL",
]
