#!/usr/bin/env python3
"""
synthetic.py  ·  Realistic Synthetic Dataset Generator
=======================================================
Produces TESS-like light curves with physically consistent injected
signals for honest recall / FPR evaluation of the Axiom-ZSpace pipeline.

Two target families:
  * "true"  : a real planet (Mandel-Agol transit) embedded in TESS-like noise
  * "false" : a contamination source that must NOT be certified as a planet
              (eclipsing binary, stellar rotation/activity, single-event dip,
               grazing EB, pure noise, harmonic alias traps)

Noise model mirrors TESS 2-min cadence SPOC pipeline:
  * per-sector continuous cadence 120 s
  * white gaussian noise (per-cadence, tuned to a target CDPP-equivalent)
  * red noise (1/f) + slow systematics trend per sector
  * sector-to-sector gaps and occasional flag-like bad cadences
  * quadratic limb-darkened Mandel-Agol transit via batman
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

CADENCE_S = 120.0            # TESS 2-min cadence
CADENCE_D = CADENCE_S / 86400.0
SECTOR_DAYS = 27.0           # one TESS sector ~27 d continuous
DEFAULT_RNG_SEED = 20260814


# ─────────────────────────────────────────────────────────────────────────────
# Noise: TESS-like correlated + white photometry
# ─────────────────────────────────────────────────────────────────────────────

def red_noise(n: int, alpha: float = 1.0, seed: Optional[int] = None,
              cadence_d: float = CADENCE_D, sector_days: float = SECTOR_DAYS) -> np.ndarray:
    """
    Correlated (pink-ish) noise via a stationary AR(1) process, per sector.
    FFT colouring is NOT used: on the gappy multi-sector axis an FFT-built
    spectrum develops a sharp quasi-periodic tone (measured: a spurious BLS
    peak at ~2.18 d with SNR>13 on pure noise) that swamps real transits.
    AR(1) gives broad 1/f-like power without periodic ringing.
    Correlation time is tuned to ~1 h (CDPP-like) — the TESS red-noise
    carpet operates on intra-hour timescales. Longer correlation or larger
    red fraction seeds spurious BLS peaks (measured SNR>8 on pure noise for
    tau>2 h) that would make the noise targets detectable. Returns
    zero-mean, unit-RMS series.
    """
    rng = np.random.default_rng(seed)
    n_sectors = max(1, int(round(n * cadence_d / sector_days)) + 1)
    pts = n // n_sectors
    rho = 1.0 - cadence_d / 0.042          # ~1 h correlation time
    rho = max(rho, 0.0)
    chunks = []
    for s in range(n_sectors):
        k = pts if s < n_sectors - 1 else n - pts * (n_sectors - 1)
        sig = 1.0
        x = np.empty(k)
        x[0] = rng.standard_normal() * sig / math.sqrt(max(1 - rho * rho, 1e-9))
        for i in range(1, k):
            x[i] = rho * x[i - 1] + rng.standard_normal() * sig
        chunks.append(x)
    colored = np.concatenate(chunks)
    std = np.std(colored) + 1e-12
    return (colored - np.mean(colored)) / std


def systematics_trend(t: np.ndarray, seed: Optional[int] = None,
                      sector_days: float = SECTOR_DAYS) -> np.ndarray:
    """
    Slow per-sector instrumental systematic. Deliberately aperiodic: a
    gentle quadratic drift with random sign/curvature, NO sinusoid — a
    periodic systematic at a period inside the BLS search window would be
    detected as a spurious transit (measured SNR>80 on pure systematics).
    """
    rng = np.random.default_rng(seed)
    n = t.size
    out = np.zeros(n)
    sector_boundaries = np.arange(t[0], t[-1] + sector_days, sector_days)
    for s0, s1 in zip(sector_boundaries[:-1], sector_boundaries[1:]):
        m = (t >= s0) & (t < s1)
        if m.sum() == 0:
            continue
        tt = (t[m] - s0) / (s1 - s0)
        amp = (0.5e-4 + 1.2e-4 * rng.random()) * rng.choice([-1.0, 1.0])
        curve = amp * (tt - 0.5) ** 2 - amp / 12.0      # zero-mean quadratic
        out[m] = curve
    return out


def make_noise(t: np.ndarray, white_ppm: float, seed: Optional[int] = None,
               red_alpha: float = 1.0, red_frac: float = 0.05) -> np.ndarray:
    """
    Combine white + red noise into flux units around 1.0.
    white_ppm : per-cadence white noise RMS in ppm.
    red_frac  : fraction of total variance carried by the 1/f component.
    """
    n = t.size
    white = np.random.default_rng((seed or 0) + 1).standard_normal(n) * white_ppm * 1e-6
    rng = np.random.default_rng((seed or 0) + 2)
    red = red_noise(n, alpha=red_alpha, seed=(seed or 0) + 3)
    sys = systematics_trend(t, seed=(seed or 0) + 4)
    # scale so total variance matches white_ppm budget: var = var_w + var_r
    red_target = red_frac * (white_ppm * 1e-6) ** 2
    if red_frac > 0:
        red = red * math.sqrt(red_target) / (np.std(red) + 1e-12)
    sys = sys - np.mean(sys)
    flux = 1.0 + white + red + sys
    return flux


# ─────────────────────────────────────────────────────────────────────────────
# Time axis: multi-sector TESS-like with gaps + flagged cadences
# ─────────────────────────────────────────────────────────────────────────────

def make_time_axis(n_sectors: int = 3, seed: Optional[int] = None) -> tuple[np.ndarray, np.ndarray]:
    """
    Build a TESS-like time axis: `n_sectors` continuous sectors of SECTOR_DAYS
    with gaps (data downlink) between them. Returns (time, quality_mask)
    with ~0.5% flagged bad cadences scattered through.
    """
    rng = np.random.default_rng(seed)
    times = []
    t0 = 1200.0
    gap = 1.0  # day gap between sectors
    for s in range(n_sectors):
        # Slightly randomise sector length (23-29 d) so the gappy axis has
        # no exact periodicity that could seed a spurious BLS tone.
        sector_len = SECTOR_DAYS * rng.uniform(0.85, 1.08)
        n_cad = int(round(sector_len / CADENCE_D))
        t = t0 + np.arange(n_cad) * CADENCE_D
        t += rng.normal(0, 0.05 * CADENCE_D, n_cad)  # tiny jitter
        times.append(t)
        t0 = times[-1][-1] + gap * rng.uniform(0.7, 1.4)
    t = np.concatenate(times)
    # ~0.5% flagged cadences (quality != 0) scattered, plus ~0.3% upward spikes
    n = t.size
    quality = np.zeros(n, dtype=int)
    n_flag = int(round(0.005 * n))
    idx = rng.choice(n, size=n_flag, replace=False)
    quality[idx] = 256
    return t, quality


# ─────────────────────────────────────────────────────────────────────────────
# Transit injection (Mandel-Agol via batman)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class InjectedTransit:
    period_days:   float
    t0_days:       float
    depth:         float           # dimensionless (Rp/Rs)^2 approx
    duration_hrs:  float
    rp_rs:         float
    a_rs:          float
    inclination:   float
    u1:            float
    u2:            float


def mandel_agol_transit(t: np.ndarray, inj: InjectedTransit) -> np.ndarray:
    """
    Full quadratic limb-darkened transit model via batman.
    Returns relative flux (1 outside transit, <1 inside).
    """
    import batman

    params = batman.TransitParams()
    params.t0 = inj.t0_days
    params.per = inj.period_days
    params.rp = inj.rp_rs
    params.a = inj.a_rs
    params.inc = inj.inclination
    params.ecc = 0.0
    params.w = 90.0
    params.u = [inj.u1, inj.u2]
    params.limb_dark = "quadratic"
    model = batman.TransitModel(params, t)
    return model.light_curve(params)


# ─────────────────────────────────────────────────────────────────────────────
# Target parameter generation
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class SyntheticTarget:
    kind: str                      # "true" | "false"
    subkind: str                   # e.g. "planet", "eb", "rotation", ...
    tic_id: str
    time: np.ndarray
    flux: np.ndarray               # normalised, ~1.0
    quality: np.ndarray
    injected: Optional[InjectedTransit] = None
    label_period: Optional[float] = None   # ground truth period if any
    injected_depth: Optional[float] = None
    meta: dict = field(default_factory=dict)
    stellar: dict = field(default_factory=dict)


# Physical consistency for generated planets
STELLAR_BANK = [
    # M0      M=0.60 R=0.60 Teff=3900
    dict(st_mass=0.60, st_rad=0.60, st_teff=3900.0, st_logg=4.66),
    # M3      M=0.36 R=0.36 Teff=3400
    dict(st_mass=0.36, st_rad=0.36, st_teff=3400.0, st_logg=4.87),
    # M5      M=0.21 R=0.21 Teff=3100
    dict(st_mass=0.21, st_rad=0.21, st_teff=3100.0, st_logg=5.10),
    # K5      M=0.70 R=0.72 Teff=4400
    dict(st_mass=0.70, st_rad=0.72, st_teff=4400.0, st_logg=4.57),
    # K2      M=0.80 R=0.80 Teff=5000
    dict(st_mass=0.80, st_rad=0.80, st_teff=5000.0, st_logg=4.55),
    # G2      M=1.00 R=1.00 Teff=5770
    dict(st_mass=1.00, st_rad=1.00, st_teff=5770.0, st_logg=4.44),
]


def random_planet_params(rng: np.random.Generator, period_days: float,
                         stellar: dict) -> InjectedTransit:
    """
    Derive a self-consistent transit geometry for a planet at `period_days`
    around `stellar` (Kepler + geometry). Depth from a drawn radius ratio.
    """
    G = 6.674e-11
    M_SUN = 1.989e30
    R_SUN = 6.957e8
    DAY = 86400.0
    AU = 1.496e11

    st_mass = stellar["st_mass"]
    st_rad = stellar["st_rad"]

    # semi-major axis (m) via Kepler III
    a_m = (G * st_mass * M_SUN * (period_days * DAY) ** 2 / (4 * math.pi ** 2)) ** (1 / 3.0)
    a_rs = a_m / (st_rad * R_SUN)

    # planet radius: 0.8 - 4.5 Earth radii
    R_EARTH = 6.371e6
    rp_earth = rng.uniform(0.8, 4.5)
    rp_rs = (rp_earth * R_EARTH) / (st_rad * R_SUN)

    # depth approx (Rp/Rs)^2 (quadratic LD slightly reduces)
    depth = rp_rs ** 2
    if depth > 0.02:
        depth = 0.02  # cap; very deep "planets" are unphysical for stellar hosts

    # inclination: random near-transit; duration from geometry
    cos_i = rng.uniform(0.0, 0.75)      # b = (a/Rs)cos i < 0.75
    inclination = math.degrees(math.acos(cos_i / a_rs)) if a_rs > 0 else 90.0
    inclination = min(90.0, inclination)

    # duration approx: T ≈ (P/π) * (R*/a) * sqrt(1-b^2)
    b = a_rs * math.cos(math.radians(inclination))
    b = min(b, 1.0 - 1e-3)
    T_frac = (1.0 / a_rs) * math.sqrt(1.0 - b ** 2) / math.pi
    duration_hrs = max(0.5, period_days * 24.0 * T_frac)

    u1 = rng.uniform(0.2, 0.45)
    u2 = rng.uniform(0.05, 0.20)
    t0 = rng.uniform(period_days * 0.0, period_days)

    return InjectedTransit(
        period_days=period_days, t0_days=t0, depth=depth, duration_hrs=duration_hrs,
        rp_rs=rp_rs, a_rs=a_rs, inclination=inclination, u1=u1, u2=u2,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Noise calibration: target SNR
# ─────────────────────────────────────────────────────────────────────────────

def calibrate_white_ppm(inj: InjectedTransit, target_snr: float,
                        t: np.ndarray, rng: np.random.Generator) -> float:
    """
    Choose the white-noise level so the injected transit achieves ~target_snr
    in a phase-folded detection sense.

    BLS-ish SNR estimate:
        SNR ≈ depth / sigma_white_pt * sqrt(N_in_transit) * sqrt(N_transits)

    Returns white_ppm such that the estimate matches target_snr (clamped).
    """
    n_in = max(1.0, inj.duration_hrs * 3600.0 / CADENCE_S)
    n_tr = max(1.0, (t[-1] - t[0]) / inj.period_days)
    sigma_needed = inj.depth * math.sqrt(n_in) * math.sqrt(n_tr) / target_snr
    white_ppm = sigma_needed * 1e6
    white_ppm = float(np.clip(white_ppm, 20.0, 50000.0))
    return white_ppm


# ─────────────────────────────────────────────────────────────────────────────
# Generators
# ─────────────────────────────────────────────────────────────────────────────

def _base_lc(n_sectors: int = 5, seed: Optional[int] = None):
    t, quality = make_time_axis(n_sectors=n_sectors, seed=seed)
    return t, quality


def generate_true_planet(idx: int, period_days: float, target_snr: float,
                         seed: Optional[int] = None,
                         n_sectors: int = 5) -> SyntheticTarget:
    rng = np.random.default_rng((seed or 0) + idx * 101)
    t, quality = _base_lc(n_sectors=n_sectors, seed=(seed or 0) + idx * 7)

    stellar = STELLAR_BANK[idx % len(STELLAR_BANK)]
    inj = random_planet_params(rng, period_days, stellar)

    flux = np.ones_like(t)
    transit = mandel_agol_transit(t, inj)
    flux = flux * transit

    white_ppm = calibrate_white_ppm(inj, target_snr, t, rng)
    noise = make_noise(t, white_ppm, seed=(seed or 0) + idx * 13)
    flux = flux + (noise - 1.0)

    tic_id = f"SYN{100000 + idx}"
    return SyntheticTarget(
        kind="true", subkind="planet", tic_id=tic_id, time=t, flux=flux, quality=quality,
        injected=inj, label_period=period_days, injected_depth=inj.depth,
        meta=dict(target_snr=target_snr, white_ppm=white_ppm, n_sectors=n_sectors),
        stellar=stellar,
    )


def generate_false_eb(idx: int, period_days: float, depth: float,
                      seed: Optional[int] = None,
                      n_sectors: int = 5) -> SyntheticTarget:
    """Eclipsing binary: two transits of different depths (even/odd mismatch)."""
    rng = np.random.default_rng((seed or 0) + idx * 11)
    t, quality = _base_lc(n_sectors=n_sectors, seed=(seed or 0) + idx * 7)

    flux = np.ones_like(t)
    phase = ((t / period_days) % 1.0)
    dur = 0.02
    even = (np.floor((t / period_days)) % 2) == 0
    dip1 = (phase < dur) | (phase > 1 - dur)          # primary
    dip2 = ((phase > 0.5 - dur) & (phase < 0.5 + dur))  # secondary
    flux[dip1 & even] -= depth
    flux[dip1 & ~even] -= depth * 0.5    # even/odd differ -> EB fingerprint
    flux[dip2] -= depth * 0.45

    noise = make_noise(t, 300.0, seed=(seed or 0) + idx * 13)
    flux = flux + (noise - 1.0)
    return SyntheticTarget(
        kind="false", subkind="eb", tic_id=f"SYN{200000 + idx}",
        time=t, flux=flux, quality=quality,
        label_period=period_days, injected_depth=depth,
        meta=dict(n_sectors=n_sectors), stellar=STELLAR_BANK[idx % len(STELLAR_BANK)],
    )


def generate_false_rotation(idx: int, period_days: float, amplitude: float,
                            seed: Optional[int] = None,
                            n_sectors: int = 5) -> SyntheticTarget:
    """Stellar rotation / spot modulation: broad quasi-sinusoidal signal."""
    rng = np.random.default_rng((seed or 0) + idx * 17)
    t, quality = _base_lc(n_sectors=n_sectors, seed=(seed or 0) + idx * 7)

    # rotation + a harmonic to mimic star spots/active regions
    rot = amplitude * np.sin(2 * np.pi * t / period_days)
    rot += amplitude * 0.5 * np.sin(2 * np.pi * t / (period_days / 2.0) + 1.0)
    flux = 1.0 + rot
    noise = make_noise(t, 150.0, seed=(seed or 0) + idx * 13)
    flux = flux + (noise - 1.0)
    return SyntheticTarget(
        kind="false", subkind="rotation", tic_id=f"SYN{300000 + idx}",
        time=t, flux=flux, quality=quality,
        label_period=period_days, injected_depth=amplitude,
        meta=dict(n_sectors=n_sectors), stellar=STELLAR_BANK[idx % len(STELLAR_BANK)],
    )


def generate_false_single_event(idx: int, depth: float,
                                seed: Optional[int] = None,
                                n_sectors: int = 5) -> SyntheticTarget:
    """One deep dip, non-repeating (instrumental / eclipsing single event)."""
    rng = np.random.default_rng((seed or 0) + idx * 19)
    t, quality = _base_lc(n_sectors=n_sectors, seed=(seed or 0) + idx * 7)
    flux = np.ones_like(t)
    center = t[0] + 0.4 * (t[-1] - t[0])
    width = 3.0 * CADENCE_D
    dip = np.abs(t - center) < width
    flux[dip] -= depth
    noise = make_noise(t, 200.0, seed=(seed or 0) + idx * 13)
    flux = flux + (noise - 1.0)
    return SyntheticTarget(
        kind="false", subkind="single_event", tic_id=f"SYN{400000 + idx}",
        time=t, flux=flux, quality=quality, label_period=None, injected_depth=depth,
        meta=dict(n_sectors=n_sectors), stellar=STELLAR_BANK[idx % len(STELLAR_BANK)],
    )


def generate_false_noise(idx: int, white_ppm: float,
                         seed: Optional[int] = None,
                         n_sectors: int = 5) -> SyntheticTarget:
    """Pure TESS-like noise with no injected signal."""
    t, quality = _base_lc(n_sectors=n_sectors, seed=(seed or 0) + idx * 7)
    flux = make_noise(t, white_ppm, seed=(seed or 0) + idx * 13)
    return SyntheticTarget(
        kind="false", subkind="noise", tic_id=f"SYN{500000 + idx}",
        time=t, flux=flux, quality=quality, label_period=None, injected_depth=None,
        meta=dict(white_ppm=white_ppm, n_sectors=n_sectors),
        stellar=STELLAR_BANK[idx % len(STELLAR_BANK)],
    )


def generate_false_grazing_eb(idx: int, period_days: float, depth: float,
                              seed: Optional[int] = None,
                              n_sectors: int = 5) -> SyntheticTarget:
    """Grazing EB: V-shaped (no flat bottom), depths differ even/odd."""
    rng = np.random.default_rng((seed or 0) + idx * 23)
    t, quality = _base_lc(n_sectors=n_sectors, seed=(seed or 0) + idx * 7)
    flux = np.ones_like(t)
    phase = ((t / period_days) % 1.0)
    # V-shape: triangular dip at both phases
    dur = 0.03
    def vshape(ph, ph0):
        d = np.abs(((ph - ph0 + 0.5) % 1.0) - 0.5)
        return np.clip(1 - d / dur, 0, 1)
    even = (np.floor(t / period_days) % 2) == 0
    v1 = vshape(phase, 0.0)
    v2 = vshape(phase, 0.5)
    flux -= v1 * depth * np.where(even, 1.0, 0.55)
    flux -= v2 * depth * 0.5
    noise = make_noise(t, 300.0, seed=(seed or 0) + idx * 13)
    flux = flux + (noise - 1.0)
    return SyntheticTarget(
        kind="false", subkind="grazing_eb", tic_id=f"SYN{600000 + idx}",
        time=t, flux=flux, quality=quality, label_period=period_days,
        injected_depth=depth, meta=dict(n_sectors=n_sectors),
        stellar=STELLAR_BANK[idx % len(STELLAR_BANK)],
    )
