"""
auditors.py  ·  Transit Vitality Auditing (V2.0 — FP Reduction)
=================================================================
Implements the physical discriminators that build the depth and
limb-shape component scores, plus new FP reduction filters.

Audit 1 — Even/Odd Transit Test
  Split transits into even/odd numbered events.
  If depth discrepancy Δσ > 3.0 → flag as Eclipsing Binary (EB).
  Score penalty: EB flags reduce S_δ.

Audit 2 — Depth Consistency (Coefficient of Variation)
  S_δ = max(0, 1 - CV / 0.10)
  where CV = std(individual depths) / mean(individual depths)

Audit 3 — Limb Shape Discriminator (Mandel-Agol via batman)
  Fit a Mandel-Agol transit model to the phase-folded light curve.
  Compute residuals at transit centre vs wings.
  U-shaped residuals → planet (score near 1)
  V-shaped residuals → EB / systematics (score near 0)
  Score S_τ encodes U-vs-V morphology.

Audit 4 — Ingress/Egress Duration Ratio (NEW)
  Computes the ratio of ingress/egress duration to total transit duration.
  Planets produce U-shaped transits (flat bottom, short ingress/egress).
  EBs produce V-shaped transits (no flat bottom, long ingress/egress).
  If ingress_fraction > 0.45 → high FP risk (V-shape).

Audit 5 — MCMC Posterior Validation (NEW)
  Uses emcee to sample the posterior distribution of transit parameters.
  Flags non-Gaussian posteriors as noise indicators.
  Provides proper uncertainties on all derived parameters.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from scipy.optimize import minimize, curve_fit
from scipy.stats import ttest_ind


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
from zspace_engine import thresholds as _T

EVEN_ODD_SIGMA_THRESHOLD = float(_T.threshold("ev_odd_sigma_eb"))
EVEN_ODD_P_VALUE_THRESHOLD = float(_T.threshold("ev_odd_pvalue_eb"))   # Welch t-test p-value ceiling for EB flag
CV_NORMALISATION         = 0.10   # CV at which S_δ → 0 (legacy reference)
CHI2_NORMALISATION       = 4.0    # chi²_red above which depths are "inconsistent"
MIN_TRANSITS_FOR_EOTEST  = 4      # minimum transits to run even/odd test
INGRESS_FRACTION_THRESHOLD = float(_T.threshold("ingress_fraction_vshape")) # V-shape threshold for ingress/egress ratio
MCMC_NWALKERS            = 32     # number of MCMC walkers
MCMC_NSTEPS              = 500    # number of MCMC steps
MCMC_BURNIN              = 200    # burn-in steps to discard


# ─────────────────────────────────────────────────────────────────────────────
# Audit result containers
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class EvenOddResult:
    n_even:          int
    n_odd:           int
    depth_even:      float
    depth_odd:       float
    depth_even_err:  float
    depth_odd_err:   float
    delta_sigma:     float
    is_eb_flag:      bool
    proof:           str
    t_stat:          float = 0.0
    p_value:         float = 1.0


@dataclass
class SecondaryEclipseResult:
    """
    Folded secondary-eclipse measurement (EB discriminator).
    secondary_ratio : secondary_depth / primary_depth at phase 0.5 vs 0.
    A real planet → ≈0; an EB → sizable (≥ SECONDARY_RATIO_VETO).
    """
    secondary_depth:  float
    primary_depth:    float
    secondary_ratio:  float
    secondary_snr:    float
    n_secondary:      int
    proof:            str


# EB veto: a phase-0.5 eclipse at ≥(threshold) of the primary transit depth
# means the companion emits light → eclipsing binary, not a dark planet.
SECONDARY_RATIO_VETO = float(_T.threshold("secondary_ratio_veto"))


@dataclass
class DepthConsistencyResult:
    depths:          np.ndarray
    mean_depth:      float
    std_depth:       float
    cv:              float
    s_depth:         float
    proof:           str
    flags:           List[str] = field(default_factory=list)


@dataclass
class LimbShapeResult:
    """
    Mandel-Agol fit result with V/U shape discrimination.
    model_params: fitted [rp/rs, a/rs, inc, u1, u2]
    residual_centre / residual_wings: RMS of residuals in each zone
    shape_ratio: residual_wings / residual_centre  (> 1 → U-shape → planet)
    s_limb: score [0, 1]
    """
    rp_rs:              float       # radius ratio R_p / R_★
    a_rs:               float       # scaled semi-major axis a / R_★
    inclination_deg:    float
    u1:                 float       # limb darkening coefficient 1
    u2:                 float       # limb darkening coefficient 2
    residual_rms:       float       # overall model residual RMS
    residual_centre:    float       # RMS in central 20% of transit duration
    residual_wings:     float       # RMS in outer 40% of transit duration
    shape_ratio:        float       # wings_rms / centre_rms
    s_limb:             float
    proof:              str
    flags:              List[str] = field(default_factory=list)


@dataclass
class IngressEgressResult:
    """
    Result of the ingress/egress duration ratio test.
    
    Planets (U-shape): ingress_fraction ≈ 0.1-0.3 (short ingress, flat bottom)
    EBs (V-shape):     ingress_fraction ≈ 0.4-0.5 (long ingress, no flat bottom)
    """
    ingress_duration:    float    # hours
    egress_duration:     float    # hours
    flat_duration:       float    # hours
    total_duration:      float    # hours
    ingress_fraction:    float    # ingress_dur / total_dur
    flat_fraction:       float    # flat_dur / total_dur
    is_v_shape:          bool     # True if V-shaped (EB candidate)
    fp_risk:             str      # "LOW", "MEDIUM", "HIGH"
    proof:               str
    flags:               List[str] = field(default_factory=list)


@dataclass
class MCMCResult:
    """
    Result of MCMC posterior sampling.
    
    If posteriors are non-Gaussian (high skewness/kurtosis),
    this indicates the transit model is fitting noise rather
    than a coherent signal.
    """
    n_walkers:           int
    n_steps:             int
    acceptance_fraction: float
    params_median:       dict      # median of each parameter
    params_stddev:       dict      # std dev of each parameter  
    is_gaussian:         bool      # True if posteriors are roughly Gaussian
    skewness_max:        float     # max |skewness| across parameters
    kurtosis_max:        float     # max |excess kurtosis| across parameters
    gelman_rubin_max:    float     # max Gelman-Rubin statistic
    proof:               str
    flags:               List[str] = field(default_factory=list)


# ─────────────────────────────────────────────────────────────────────────────
# Individual transit depth extractor
# ─────────────────────────────────────────────────────────────────────────────

def extract_individual_transit_depths(
    time:     np.ndarray,
    flux:     np.ndarray,
    period:   float,
    t0:       float,
    duration: float,
    cadence_est: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Measure the depth of each individual transit event using template
    matched-filtering.

    A high-SNR phase-folded template of the transit shape is built from the
    full light curve, then each epoch is aligned to it by cross-correlation.
    This is far more robust than window-mean or argmin anchoring: at low
    per-transit SNR a single least-squares fitted amplitude against the
    template enjoys the noise-averaging power of every in-transit point.

    Returns
    -------
    depths     : array of per-transit depths
    transit_ns : array of transit numbers (0=first, 1=second, ...)
    depth_errs : array of per-transit depth uncertainties
    """
    # Transit epoch numbers
    n_min = int(np.floor((time[0] - t0) / period))
    n_max = int(np.ceil((time[-1] - t0) / period)) + 1
    epochs = np.arange(n_min, n_max + 1)

    ts_sorted = np.sort(time)
    if cadence_est is None or cadence_est <= 0:
        dts = np.diff(ts_sorted)
        cadence_est = float(np.median(dts[dts > 0])) if dts.size else duration
    cadence_est = max(cadence_est, 1e-9)

    half_dur      = duration / 2.0
    oot_half      = duration * 3.0      # out-of-transit window half-width
    temp_half     = half_dur * 2.0      # template half-width (phases beyond are 1.0)
    phase_half    = half_dur / period
    temp_phase_h  = temp_half / period

    # ── Build folded template on the joint (phase) grid ──────────────────────
    phases = ((time - t0) / period) % 1.0
    phases[phases > 0.5] -= 1.0
    in_temp = np.abs(phases) <= temp_phase_h
    templ_x_b = 100                      # bins across the transit window
    bin_edges = np.linspace(-temp_phase_h, temp_phase_h, templ_x_b + 1)
    bin_idx   = np.digitize(phases[in_temp], bin_edges) - 1
    bin_idx   = np.clip(bin_idx, 0, templ_x_b - 1)
    t_dip = np.zeros(templ_x_b)
    t_cnt = np.zeros(templ_x_b)
    for b in range(templ_x_b):
        m = bin_idx == b
        if m.sum() > 0:
            t_dip[b] = float(np.mean(flux[in_temp][m]))
            t_cnt[b] = float(m.sum())
    good = t_cnt > 0
    if good.sum() < templ_x_b * 0.4:
        # Not enough in-transit phase coverage → fall back to a top-hat.
        t_dip[:] = 1.0 - 0.0
        return np.array([]), np.array([]), np.array([])
    t_dip[~good] = np.interp(np.flatnonzero(~good), np.flatnonzero(good), t_dip[good])
    depth_ref = float(np.mean(1.0 - t_dip))
    if depth_ref <= 0:
        depth_ref = 1e-8
    # Normalised dip template: 1.0 at full immersion, 0.0 far from transit.
    t_dip_norm = (1.0 - t_dip) / depth_ref
    t_xs = bin_edges[:-1] + 0.5 * (bin_edges[1] - bin_edges[0])

    def template_at(phase: np.ndarray) -> np.ndarray:
        """Return the normalised dip template evaluated at given phases."""
        return np.clip(np.interp(phase, t_xs, t_dip_norm, left=0.0, right=0.0), 0.0, None)

    depths     = []
    depth_errs = []
    ns_used    = []
    shift_step = cadence_est
    shifts     = np.arange(-temp_half, temp_half + shift_step, shift_step)

    for n in epochs:
        t_centre = t0 + n * period
        # In/out masks around predicted centre (expand the outer mask to the
        # template coverage; alignment shifts are absorbed by the template).
        in_mask = np.abs(time - t_centre) <= temp_half
        search  = np.abs(time - t_centre) <= oot_half
        if in_mask.sum() < 3 or search.sum() < 6:
            continue

        ft = flux[in_mask]
        tt = time[in_mask]

        # Cross-correlate the in-transit residual against the template over
        # candidate shifts — matches the epoch's true centre robustly.
        best = None
        for s in shifts:
            rel = (tt - (t_centre + s)) / period
            tpl = template_at(rel)
            num = float(np.sum(tpl * (1.0 - ft)))
            den = float(np.sum(tpl * tpl)) + 1e-12
            if den <= 0:
                continue
            if best is None or num / den > best[0]:
                best = (num / den, s)
        if best is None:
            continue
        _, shift = best
        t_align = t_centre + shift

        # Depth = amplitude of matched-filter fit over the aligned window.
        inm = np.abs(time - t_align) <= temp_half
        if inm.sum() < 3:
            continue
        rel = (time[inm] - t_align) / period
        tpl = template_at(rel)
        oot_mask = (
            (np.abs(time - t_align) > temp_half) &
            (np.abs(time - t_align) <= oot_half)
        )
        if oot_mask.sum() < 3:
            continue

        f_in  = flux[inm]
        f_oot = flux[oot_mask]
        baseline = float(np.median(f_oot))
        if baseline <= 0:
            continue

        yy = 1.0 - f_in / baseline          # in-transit dip (flux units)
        den = float(np.sum(tpl * tpl)) + 1e-12
        depth = float(np.sum(tpl * yy) / den)
        if depth < -0.01 or depth > 0.5:
            continue

        # Uncertainty: propagate per-point noise sigma_oot through the same
        # linear filter, plus the two-mean baseline term.
        sigma_oot = float(np.std(f_oot, ddof=1))
        depth_err = (sigma_oot / baseline) / den * float(np.sqrt(np.sum(tpl * tpl)))
        depth_err = float(np.sqrt(depth_err ** 2 + (sigma_oot / baseline) ** 2))
        if not math.isfinite(depth_err) or depth_err <= 0:
            continue

        depths.append(depth)
        depth_errs.append(depth_err)
        ns_used.append(n)

    return np.array(depths), np.array(ns_used), np.array(depth_errs)


# ─────────────────────────────────────────────────────────────────────────────
# Ingress/Egress Duration Ratio Test (NEW)
# ─────────────────────────────────────────────────────────────────────────────

class IngressEgressTest:
    """
    Measures the ratio of ingress/egress duration to total transit duration.
    
    Physical motivation:
    - Planets: short ingress (ratio ~0.1-0.3), long flat bottom (U-shape)
    - EBs: long ingress (ratio ~0.4-0.5), no flat bottom (V-shape)
    
    The transit is fitted with a trapezoid model to extract:
    - T_ingress: time from first contact to full immersion
    - T_flat: time of full immersion (flat bottom)
    - T_egress: time from full immersion end to last contact
    
    ingress_fraction = T_ingress / T_total
    """

    @staticmethod
    def test(
        bin_phase: np.ndarray,
        bin_flux: np.ndarray,
        period: float,
        duration: float,
        transit_depth: float,
    ) -> IngressEgressResult:
        """
        Fit a trapezoid to the phase-folded transit and measure
        ingress/egress fractions.
        """
        flags: List[str] = []
        valid = np.isfinite(bin_flux)
        
        if valid.sum() < 10:
            return IngressEgressResult(
                ingress_duration=0, egress_duration=0, flat_duration=0,
                total_duration=duration * 24, ingress_fraction=0,
                flat_fraction=0.5, is_v_shape=False, fp_risk="UNKNOWN",
                proof="INGRESS_EGRESS | insufficient data → SKIPPED",
                flags=["INSUFFICIENT_DATA"],
            )

        ph = bin_phase[valid]
        fl = bin_flux[valid]
        half_dur = (duration / period) / 2.0

        # Fit trapezoid model: parameters = [depth, ingress_frac, flat_frac]
        def trapezoid_model(phase, depth, ingress_f, flat_f):
            model = np.ones_like(phase)
            ingress_w = ingress_f * half_dur
            flat_w = flat_f * half_dur
            for i, p in enumerate(phase):
                ap = abs(p)
                if ap >= half_dur:
                    model[i] = 1.0
                elif ap >= flat_w:
                    frac = (ap - flat_w) / max(ingress_w, 1e-9)
                    model[i] = 1.0 - depth * (1.0 - np.clip(frac, 0, 1))
                else:
                    model[i] = 1.0 - depth
            return model

        try:
            popt, pcov = curve_fit(
                trapezoid_model, ph, fl,
                p0=[transit_depth, 0.2, 0.6],
                bounds=([0, 0.01, 0.01], [0.5, 0.99, 0.99]),
                maxfev=3000,
            )
            depth_fit, ingress_fit, flat_fit = popt

            # Ensure ingress + flat ≤ 1.0
            if ingress_fit + flat_fit > 1.0:
                scale = 0.99 / (ingress_fit + flat_fit)
                ingress_fit *= scale
                flat_fit *= scale

            # Convert to hours
            total_dur_hrs = duration * 24.0
            ingress_hrs = ingress_fit * total_dur_hrs
            flat_hrs = flat_fit * total_dur_hrs
            egress_hrs = ingress_hrs  # Symmetric transit assumption

            ingress_fraction = ingress_fit
            flat_fraction = flat_fit

            # V-shape classification
            is_v_shape = ingress_fraction > INGRESS_FRACTION_THRESHOLD

            if ingress_fraction > 0.45:
                fp_risk = "HIGH"
                flags.append(f"V_SHAPE_HIGH_RISK | ingress_frac={ingress_fraction:.3f}")
            elif ingress_fraction > 0.35:
                fp_risk = "MEDIUM"
                flags.append(f"V_SHAPE_MEDIUM_RISK | ingress_frac={ingress_fraction:.3f}")
            else:
                fp_risk = "LOW"

            proof = (
                f"INGRESS_EGRESS | depth={depth_fit:.5f} | "
                f"ingress_frac={ingress_fraction:.3f}, flat_frac={flat_fraction:.3f} | "
                f"T_ingress={ingress_hrs:.2f}h, T_flat={flat_hrs:.2f}h, T_total={total_dur_hrs:.2f}h | "
                f"shape={'V-SHAPE (EB)' if is_v_shape else 'U-SHAPE (PLANET)'} | "
                f"FP_risk={fp_risk}"
            )

            return IngressEgressResult(
                ingress_duration=ingress_hrs,
                egress_duration=egress_hrs,
                flat_duration=flat_hrs,
                total_duration=total_dur_hrs,
                ingress_fraction=ingress_fraction,
                flat_fraction=flat_fraction,
                is_v_shape=is_v_shape,
                fp_risk=fp_risk,
                proof=proof,
                flags=flags,
            )

        except Exception as exc:
            return IngressEgressResult(
                ingress_duration=0, egress_duration=0, flat_duration=0,
                total_duration=duration * 24, ingress_fraction=0,
                flat_fraction=0.5, is_v_shape=False, fp_risk="UNKNOWN",
                proof=f"INGRESS_EGRESS | fit failed: {exc} → SKIPPED",
                flags=["FIT_FAILED"],
            )


# ─────────────────────────────────────────────────────────────────────────────
# MCMC Posterior Validator (NEW)
# ─────────────────────────────────────────────────────────────────────────────

class MCMCValidator:
    """
    Uses Markov Chain Monte Carlo (emcee) to sample the posterior distribution
    of transit parameters.  Flags non-Gaussian posteriors as noise indicators.

    Physical rationale:
    A real transit signal produces well-constrained, approximately Gaussian
    posteriors.  Noise or systematic artifacts produce:
    - Multi-modal posteriors
    - Highly skewed distributions
    - Very broad (unconstrained) parameters
    
    We flag candidates where |skewness| > 1.0 or |kurtosis| > 6.0
    as potentially fitting noise rather than real signals.
    """

    @staticmethod
    def validate(
        bin_phase: np.ndarray,
        bin_flux: np.ndarray,
        period: float,
        duration: float,
        transit_depth: float,
        n_steps: int = 0,
    ) -> MCMCResult:
        """
        Run MCMC sampling on transit parameters and validate posteriors.

        Parameters
        ----------
        n_steps : int, optional
            Override the default MCMC_NSTEPS. If 0, uses the module default.
        """
        try:
            import emcee
        except ImportError:
            return MCMCResult(
                n_walkers=0, n_steps=0, acceptance_fraction=0.0,
                params_median={}, params_stddev={},
                is_gaussian=True, skewness_max=0.0, kurtosis_max=0.0,
                gelman_rubin_max=1.0,
                proof="MCMC | emcee not installed → SKIPPED",
                flags=["EMCEE_NOT_AVAILABLE"],
            )

        valid = np.isfinite(bin_flux)
        if valid.sum() < 20:
            return MCMCResult(
                n_walkers=0, n_steps=0, acceptance_fraction=0.0,
                params_median={}, params_stddev={},
                is_gaussian=True, skewness_max=0.0, kurtosis_max=0.0,
                gelman_rubin_max=1.0,
                proof="MCMC | insufficient data points → SKIPPED",
                flags=["INSUFFICIENT_DATA"],
            )

        ph = bin_phase[valid]
        fl = bin_flux[valid]
        half_dur = (duration / period) / 2.0
        flux_err = float(np.std(fl[np.abs(ph) > half_dur * 2])) if np.sum(np.abs(ph) > half_dur * 2) > 5 else 1e-4

        # Trapezoid model for MCMC (faster than batman)
        def trapezoid_model(phase, depth, ingress_f, flat_f):
            model = np.ones_like(phase)
            ingress_w = ingress_f * half_dur
            flat_w = flat_f * half_dur
            for i, p in enumerate(phase):
                ap = abs(p)
                if ap >= half_dur:
                    model[i] = 1.0
                elif ap >= flat_w:
                    frac = (ap - flat_w) / max(ingress_w, 1e-9)
                    model[i] = 1.0 - depth * (1.0 - np.clip(frac, 0, 1))
                else:
                    model[i] = 1.0 - depth
            return model

        def log_likelihood(theta):
            depth, ingress_f, flat_f = theta
            if depth <= 0 or depth > 0.5:
                return -np.inf
            if ingress_f <= 0.01 or ingress_f >= 0.99:
                return -np.inf
            if flat_f <= 0.01 or flat_f >= 0.99:
                return -np.inf
            if ingress_f + flat_f > 0.99:
                return -np.inf

            model = trapezoid_model(ph, depth, ingress_f, flat_f)
            residuals = fl - model
            return -0.5 * np.sum((residuals / flux_err) ** 2)

        def log_prior(theta):
            depth, ingress_f, flat_f = theta
            if 0 < depth < 0.5 and 0.01 < ingress_f < 0.99 and 0.01 < flat_f < 0.99:
                if ingress_f + flat_f < 0.99:
                    return 0.0
            return -np.inf

        def log_probability(theta):
            lp = log_prior(theta)
            if not np.isfinite(lp):
                return -np.inf
            ll = log_likelihood(theta)
            if not np.isfinite(ll):
                return -np.inf
            return lp + ll

        # Initialize walkers
        ndim = 3
        nwalkers = min(MCMC_NWALKERS, max(2 * ndim + 2, 10))
        nsteps = n_steps if n_steps > 0 else MCMC_NSTEPS

        p0_depth = max(transit_depth, 1e-4)
        p0 = np.array([p0_depth, 0.2, 0.6])

        # Add small perturbations
        pos = p0 + 1e-3 * np.random.randn(nwalkers, ndim)
        # Ensure all starting positions are valid
        for i in range(nwalkers):
            pos[i, 0] = np.clip(pos[i, 0], 1e-5, 0.49)
            pos[i, 1] = np.clip(pos[i, 1], 0.02, 0.98)
            pos[i, 2] = np.clip(pos[i, 2], 0.02, 0.98)
            if pos[i, 1] + pos[i, 2] > 0.98:
                pos[i, 2] = 0.98 - pos[i, 1]

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)
                sampler.run_mcmc(pos, nsteps, progress=False)

            # Discard burn-in (proportional to step count, at most 40% of steps)
            burnin = min(MCMC_BURNIN, int(nsteps * 0.4))
            samples = sampler.get_chain(discard=burnin, flat=True)
            acceptance = float(np.mean(sampler.acceptance_fraction))

            if len(samples) < 50:
                return MCMCResult(
                    n_walkers=nwalkers, n_steps=nsteps,
                    acceptance_fraction=acceptance,
                    params_median={}, params_stddev={},
                    is_gaussian=True, skewness_max=0.0, kurtosis_max=0.0,
                    gelman_rubin_max=1.0,
                    proof=f"MCMC | insufficient samples after burn-in (n={len(samples)}) → SKIPPED",
                    flags=["LOW_SAMPLES"],
                )

            # Compute statistics
            from scipy.stats import skew, kurtosis as scipy_kurtosis
            param_names = ['depth', 'ingress_frac', 'flat_frac']
            medians = {}
            stddevs = {}
            skewnesses = []
            kurtoses = []

            for i, name in enumerate(param_names):
                chain = samples[:, i]
                medians[name] = float(np.median(chain))
                stddevs[name] = float(np.std(chain))
                skewnesses.append(float(abs(skew(chain))))
                kurtoses.append(float(abs(scipy_kurtosis(chain))))

            skewness_max = max(skewnesses)
            kurtosis_max = max(kurtoses)

            # Gaussianity check:
            # |skewness| > 1.0 or |kurtosis| > 6.0 → non-Gaussian
            is_gaussian = skewness_max < 1.0 and kurtosis_max < 6.0

            flags = []
            if not is_gaussian:
                flags.append(
                    f"NON_GAUSSIAN_POSTERIOR | skew_max={skewness_max:.2f}, "
                    f"kurt_max={kurtosis_max:.2f}"
                )
            if acceptance < 0.15:
                flags.append(f"LOW_ACCEPTANCE | {acceptance:.3f}")
            if acceptance > 0.75:
                flags.append(f"HIGH_ACCEPTANCE | {acceptance:.3f}")

            proof = (
                f"MCMC | walkers={nwalkers}, steps={nsteps}, burn-in={burnin} | "
                f"acceptance={acceptance:.3f} | "
                f"depth={medians.get('depth', 0):.5f}±{stddevs.get('depth', 0):.5f} | "
                f"skew_max={skewness_max:.2f}, kurt_max={kurtosis_max:.2f} | "
                f"{'GAUSSIAN — CONSISTENT SIGNAL' if is_gaussian else 'NON-GAUSSIAN — NOISE INDICATOR'}"
            )

            return MCMCResult(
                n_walkers=nwalkers,
                n_steps=nsteps,
                acceptance_fraction=acceptance,
                params_median=medians,
                params_stddev=stddevs,
                is_gaussian=is_gaussian,
                skewness_max=skewness_max,
                kurtosis_max=kurtosis_max,
                gelman_rubin_max=1.0,  # Simplified for now
                proof=proof,
                flags=flags,
            )

        except Exception as exc:
            return MCMCResult(
                n_walkers=nwalkers, n_steps=nsteps,
                acceptance_fraction=0.0,
                params_median={}, params_stddev={},
                is_gaussian=True, skewness_max=0.0, kurtosis_max=0.0,
                gelman_rubin_max=1.0,
                proof=f"MCMC | sampling failed: {exc} → SKIPPED",
                flags=["MCMC_FAILED"],
            )


# ─────────────────────────────────────────────────────────────────────────────
# Transit Auditor
# ─────────────────────────────────────────────────────────────────────────────

class TransitAuditor:
    """
    Executes all vitality audits and returns:
      - EvenOddResult
      - DepthConsistencyResult  →  S_δ
      - LimbShapeResult         →  S_τ
      - IngressEgressResult     →  V/U shape
      - MCMCResult              →  Posterior validation
    """

    def __init__(self, verbose: bool = True, run_mcmc: bool = False) -> None:
        self.verbose = verbose
        self.run_mcmc = run_mcmc
        self._ie_test = IngressEgressTest()
        self._mcmc_validator = MCMCValidator()

    # ── Audit 1: Even/Odd ────────────────────────────────────────────────────

    def secondary_eclipse_test(
        self,
        time:     np.ndarray,
        flux:     np.ndarray,
        period:   float,
        t0:       float,
        duration: float,
    ) -> SecondaryEclipseResult:
        """
        Fold at `period` and compare the mean depth at phase 0.5 (secondary
        eclipse) against the primary transit at phase 0.

        A transiting PLANET has negligible flux at phase 0.5 (the planet emits
        no light), so secondary_depth ≈ 0. An ELLIPSING BINARY shows a real
        secondary eclipse at phase 0.5 because the fainter star still emits
        light — giving a large secondary/primary depth ratio. Robust EB
        discriminator that does not rely on per-transit depth extraction.
        """
        phase = ((time - t0) / period) % 1.0
        phase[phase > 0.5] -= 1.0
        half_h = (duration / period) / 2.0

        pri_mask  = np.abs(phase) <= half_h
        sec_mask  = np.abs(np.abs(phase) - 0.5) <= half_h
        oot_mask  = (
            (np.abs(phase) > half_h * 2.5) &
            (np.abs(np.abs(phase) - 0.5) > half_h * 2.5)
        )

        n_pri, n_sec, n_oot = int(pri_mask.sum()), int(sec_mask.sum()), int(oot_mask.sum())
        if n_pri < 3 or n_oot < 3:
            return SecondaryEclipseResult(
                secondary_depth=0.0, primary_depth=0.0, secondary_ratio=0.0,
                secondary_snr=0.0, n_secondary=0, proof="SEC_ECLIPSE | insufficient phase coverage",
            )

        baseline = float(np.mean(flux[oot_mask]))
        primary_depth   = float(np.mean(flux[oot_mask]) - np.mean(flux[pri_mask]))
        secondary_depth = float(np.mean(flux[oot_mask]) - np.mean(flux[sec_mask]))

        primary_depth   = max(primary_depth, 0.0)
        secondary_depth = max(secondary_depth, 0.0)

        if primary_depth <= 0.0:
            return SecondaryEclipseResult(
                secondary_depth=secondary_depth, primary_depth=0.0, secondary_ratio=0.0,
                secondary_snr=0.0, n_secondary=n_sec, proof="SEC_ECLIPSE | no primary transit depth",
            )

        # Secondary SNR: eclipse depth against the scatter of the phase-0.5 window.
        sec_std  = float(np.std(flux[sec_mask]))
        sec_snr  = secondary_depth / max(sec_std / np.sqrt(n_sec), 1e-12)
        sec_ratio = min(secondary_depth / primary_depth, 3.0)

        proof = (
            f"SEC_ECLIPSE | fold P={period:.5f} d: primary_depth={primary_depth:.6f}"
            f" (n={n_pri}), secondary_depth={secondary_depth:.6f} (n={n_sec})"
            f" → ratio={sec_ratio:.3f}, second_SNR={sec_snr:.2f}"
        )
        return SecondaryEclipseResult(
            secondary_depth=secondary_depth, primary_depth=primary_depth,
            secondary_ratio=sec_ratio, secondary_snr=sec_snr,
            n_secondary=n_sec, proof=proof,
        )

    def even_odd_test(
        self,
        time:     np.ndarray,
        flux:     np.ndarray,
        period:   float,
        t0:       float,
        duration: float,
    ) -> EvenOddResult:
        depths, ns, _ = extract_individual_transit_depths(time, flux, period, t0, duration)

        n_total = len(depths)
        if n_total < MIN_TRANSITS_FOR_EOTEST:
            # Not enough transits — pass with neutral result
            proof = (
                f"EVEN_ODD | n_transits={n_total} < {MIN_TRANSITS_FOR_EOTEST} → "
                "INSUFFICIENT_DATA → not flagged"
            )
            return EvenOddResult(
                n_even=0, n_odd=0,
                depth_even=np.nan, depth_odd=np.nan,
                depth_even_err=np.nan, depth_odd_err=np.nan,
                delta_sigma=0.0, is_eb_flag=False, proof=proof,
            )

        even_depths = depths[ns % 2 == 0]
        odd_depths  = depths[ns % 2 == 1]

        if len(even_depths) < 2 or len(odd_depths) < 2:
            proof = "EVEN_ODD | insufficient transits per parity → not flagged"
            return EvenOddResult(0, 0, np.nan, np.nan, np.nan, np.nan, 0.0, False, proof)

        mu_e = float(np.mean(even_depths))
        mu_o = float(np.mean(odd_depths))
        sig_e = float(np.std(even_depths, ddof=1))
        sig_o = float(np.std(odd_depths,  ddof=1))

        # Welch's t-test (unequal variances, small samples): the per-transit
        # depth scatter means n is small, so the difference must be compared
        # against the t-distribution, not the normal. p-value is reported
        # honestly and the EB flag requires BOTH |t| >= 3.0 and p < 0.01.
        t_stat, p_value = ttest_ind(even_depths, odd_depths, equal_var=False)

        se_e, se_o = sig_e / np.sqrt(len(even_depths)), sig_o / np.sqrt(len(odd_depths))
        combined_err = np.sqrt(se_e**2 + se_o**2)
        delta_sigma  = abs(mu_e - mu_o) / max(combined_err, 1e-12)

        is_eb = bool(delta_sigma >= EVEN_ODD_SIGMA_THRESHOLD and p_value < EVEN_ODD_P_VALUE_THRESHOLD)
        flag_str = f"EB_FLAG={'YES' if is_eb else 'NO'}"

        proof = (
            f"EVEN_ODD | "
            f"even: n={len(even_depths)}, μ={mu_e:.5f}±{se_e:.5f} | "
            f"odd:  n={len(odd_depths)},  μ={mu_o:.5f}±{se_o:.5f} | "
            f"Welch t={abs(t_stat):.2f}, p={p_value:.4f} | "
            f"|t|={abs(t_stat):.2f}>={EVEN_ODD_SIGMA_THRESHOLD} and "
            f"p={p_value:.4f}<{EVEN_ODD_P_VALUE_THRESHOLD} → {flag_str}"
        )
        return EvenOddResult(
            n_even=len(even_depths), n_odd=len(odd_depths),
            depth_even=mu_e, depth_odd=mu_o,
            depth_even_err=se_e, depth_odd_err=se_o,
            delta_sigma=delta_sigma, is_eb_flag=is_eb, proof=proof,
            t_stat=abs(float(t_stat)), p_value=float(p_value),
        )

    # ── Audit 2: Depth Consistency ───────────────────────────────────────────

    def depth_consistency_score(
        self,
        time:     np.ndarray,
        flux:     np.ndarray,
        period:   float,
        t0:       float,
        duration: float,
        eo_result: Optional[EvenOddResult] = None,
    ) -> DepthConsistencyResult:
        depths, _, depth_errs = extract_individual_transit_depths(time, flux, period, t0, duration)
        flags: List[str] = []

        if len(depths) < 2:
            # Single transit — use a default score with a flag
            flags.append("SINGLE_TRANSIT | CV undefined, S_δ set to 0.50")
            proof = "DEPTH_CONSISTENCY | single transit → S_δ=0.50"
            return DepthConsistencyResult(
                depths=depths, mean_depth=float(np.mean(depths)) if len(depths) else 0,
                std_depth=0.0, cv=0.0, s_depth=0.50, proof=proof, flags=flags,
            )

        mu  = float(np.mean(depths))
        std = float(np.std(depths, ddof=1))
        # Empirical CV (reported for completeness / unchanged semantics).
        cv  = std / max(mu, 1e-12)

        # Chi²-reduced consistency: (1/(n-1)) Σ (d_i - <d>)² / σ_i².
        # For a genuine planet whose per-transit depths agree to within the
        # expected noise, chi2_red ≈ 1. Large chi2_red → real depth variation
        # (EB, blended secondary) that the noise model cannot explain.
        sigma_med = float(np.median(depth_errs))
        if sigma_med > 0:
            chi2_red = float(np.mean(((depths - mu) ** 2) / depth_errs ** 2))
        else:
            chi2_red = float(np.mean(((depths - mu) ** 2))) / max(std ** 2, 1e-12)
        s_delta = max(0.0, 1.0 - (chi2_red - 1.0) / (CHI2_NORMALISATION - 1.0)) if chi2_red >= 1.0 else 1.0
        s_delta = float(np.clip(s_delta, 0.0, 1.0))

        # EB flag penalty: if even/odd test flagged EB, halve the score
        if eo_result is not None and eo_result.is_eb_flag:
            s_delta *= 0.5
            flags.append(f"EB_PENALTY | S_δ halved due to even/odd EB flag")

        proof = (
            f"DEPTH_CONSISTENCY | n={len(depths)}, μ={mu:.5f}, σ={std:.5f}, "
            f"median σ_i={sigma_med:.5f} | χ²_red={chi2_red:.3f} | "
            f"S_δ = clip(1-(χ²_red-1)/{CHI2_NORMALISATION - 1:.1f},0,1) = {s_delta:.4f}"
        )
        if chi2_red > CHI2_NORMALISATION:
            flags.append(f"HIGH_CHI2 | χ²_red={chi2_red:.2f} > {CHI2_NORMALISATION:.1f} → inconsistent depths")

        return DepthConsistencyResult(
            depths=depths, mean_depth=mu, std_depth=std,
            cv=cv, s_depth=s_delta, proof=proof, flags=flags,
        )

    # ── Audit 3: Limb Shape / Mandel-Agol ────────────────────────────────────

    @staticmethod
    def _bin_transit_zoom(
        time:     np.ndarray,
        flux:     np.ndarray,
        period:   float,
        t0:       float,
        duration: float,
        n_wing:   float = 1.5,
        n_bins:   int = 64,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Phase-fold the raw light curve and re-bin only the transit window
        ([-n_wing·T, +n_wing·T]) into `n_bins` bins.

        The transit spans a phase width of only T/P ≈ 0.004 for a 1 h transit on
        a 10 d orbit; a full-orbit 200-bin fold then has ~1 bin across the whole
        transit, leaving the centre/wing residual zones nearly empty.  Zooming
        on the window recovers real shape information while keeping the
        out-of-transit wings as the "flat" reference.

        Returns (bin_phase, bin_flux) over the zoomed window (transit at 0).
        """
        phase = ((time - t0) / period) % 1.0
        phase[phase > 0.5] -= 1.0
        T_ph = max(duration / period, 1e-5)
        win  = n_wing * T_ph
        sel  = (phase >= -win) & (phase <= win)
        phs  = phase[sel]
        fls  = flux[sel]
        if phs.size < 12:
            # Too few points in/near transit — degrade to a flat window.
            return np.linspace(-win, win, n_bins), np.ones(n_bins)

        edges   = np.linspace(-win, win, n_bins + 1)
        centres = 0.5 * (edges[:-1] + edges[1:])
        bf      = np.full(n_bins, np.nan)
        for i in range(n_bins):
            m = (phs >= edges[i]) & (phs < edges[i + 1])
            if m.sum() >= 3:
                bf[i] = np.mean(fls[m])
        good = np.isfinite(bf)
        if good.sum() < 10:
            return centres, bf
        return centres[good], bf[good]

    def limb_shape_score(
        self,
        period:       float,
        duration:     float,
        transit_depth: float,
        bin_phase:    Optional[np.ndarray] = None,
        bin_flux:     Optional[np.ndarray] = None,
        time:         Optional[np.ndarray] = None,
        flux:         Optional[np.ndarray] = None,
        t0:           Optional[float] = None,
    ) -> LimbShapeResult:
        """
        Fit a Mandel-Agol transit model to the phase-folded light curve.
        Uses the batman transit modelling package.

        Shape discrimination:
          ratio = RMS_wings / RMS_centre
          ratio > 1 → residuals larger in wings → flat-bottomed (U-shape) → planet
          ratio < 1 → residuals larger in centre → V-shape → EB or systematics

        Binning: a fixed 200-bin fold across the full [-0.5, 0.5] phase gives a
        bin width of 0.005 in phase — wider than a 1 h transit on a ~10 d period
        (phase span ≈ 0.004).  A zoomed re-binning around the transit window is
        therefore used when raw `time`/`flux`/`t0` are supplied, so the centre
        and wing zones contain enough points for a meaningful residual ratio.
        """
        flags: List[str] = []

        if time is not None and flux is not None and t0 is not None:
            zoomed = self._bin_transit_zoom(time, flux, period, t0, duration)
            bin_phase, bin_flux = zoomed

        if bin_phase is None or bin_flux is None:
            flags.append("LIMB_NO_BINS | no phase bins supplied")
            return self._trapezoid_fallback(
                np.linspace(-0.02, 0.02, 41), np.ones(41),
                period, duration, transit_depth,
            )

        try:
            import batman
        except ImportError:
            # Graceful degradation: run a simple trapezoidal fit instead
            return self._trapezoid_fallback(
                bin_phase, bin_flux, period, duration, transit_depth
            )

        valid = np.isfinite(bin_flux)
        if valid.sum() < 10:
            flags.append("INSUFFICIENT_BINS | <10 valid phase bins")
            proof = "LIMB_SHAPE | insufficient bins → S_τ=0.50"
            return LimbShapeResult(
                rp_rs=0, a_rs=0, inclination_deg=90, u1=0.3, u2=0.1,
                residual_rms=0, residual_centre=0, residual_wings=0,
                shape_ratio=1.0, s_limb=0.50, proof=proof, flags=flags,
            )

        ph = bin_phase[valid]
        fl = bin_flux[valid]

        # Initial parameter guesses from BLS result
        rp0  = np.sqrt(max(transit_depth, 1e-6))
        dur_phase = duration / period
        # Estimate a/R★ from transit duration geometry. For a full transit:
        #   T ≈ (P/π) · (R★/a) · sqrt(1-b²)   →   a/R★ ≈ P/(π·T)
        # in units where T is in days. This stays well inside the fit bounds
        # for both short and long-orbit cases (the old 2/(π·dur_phase) form
        # blew up to ~160 for a ~1 h transit on a ~10 d orbit, pinning the
        # optimiser at the a_max bound → degenerate wings/centre residuals).
        T_days   = duration / max(period, 1e-9) * period          # duration in days
        a_est    = period / (np.pi * max(T_days, 1e-6))
        a0       = float(np.clip(a_est, 2.0, 150.0))
        inc0 = 89.0

        # Build batman params
        params = batman.TransitParams()
        params.t0  = 0.0
        params.per = 1.0        # everything is in phase units
        params.rp  = rp0
        params.a   = a0
        params.inc = inc0
        params.ecc = 0.0
        params.w   = 90.0
        params.u   = [0.3, 0.1]
        params.limb_dark = "quadratic"

        # Phase grid for model
        t_model = batman.TransitModel(params, ph)

        def residual_rms(x: np.ndarray) -> float:
            rp, a, inc, u1, u2 = x
            if rp <= 0 or a < 1.0 or not (0 < inc <= 90) or u1 < 0 or u2 < 0:
                return 1e9
            params.rp  = rp
            params.a   = a
            params.inc = inc
            params.u   = [u1, u2]
            try:
                model_flux = t_model.light_curve(params)
            except Exception:
                return 1e9
            return float(np.sum((fl - model_flux) ** 2))

        x0     = np.array([rp0, a0, inc0, 0.3, 0.1])
        bounds = [(0.001, 0.5), (1.0, 200.0), (60.0, 90.0), (0.0, 1.0), (0.0, 1.0)]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = minimize(residual_rms, x0, method="L-BFGS-B", bounds=bounds,
                              options={"maxiter": 500, "ftol": 1e-12})

        rp_fit, a_fit, inc_fit, u1_fit, u2_fit = result.x

        # Guard: non-finite / out-of-range fitted limb-darkening coefficients
        # are unphysical — refuse them and fall back to the prior values so
        # the batman model stays well-posed (batman raises for invalid u).
        if not (np.isfinite([u1_fit, u2_fit]).all() and 0.0 <= u1_fit <= 1.0 and 0.0 <= u2_fit <= 1.0):
            flags.append(
                f"NONFINITE_LIMB_DARK | fitted u1={u1_fit:.4f}, u2={u2_fit:.4f} → "
                f"fallback to prior [0.3, 0.1]"
            )
            u1_fit, u2_fit = 0.3, 0.1

        params.rp  = rp_fit
        params.a   = a_fit
        params.inc = inc_fit
        params.u   = [u1_fit, u2_fit]

        model_flux = t_model.light_curve(params)
        resid      = fl - model_flux
        rms_total  = float(np.sqrt(np.mean(resid**2)))

        # Zonal residual analysis
        half_dur_phase = (duration / period) / 2.0
        centre_mask = np.abs(ph) <= 0.20 * half_dur_phase * 2
        wings_mask  = (np.abs(ph) > 0.30 * half_dur_phase * 2) & (np.abs(ph) <= half_dur_phase * 2)

        rms_c = float(np.sqrt(np.mean(resid[centre_mask]**2))) if centre_mask.sum() > 2 else rms_total
        rms_w = float(np.sqrt(np.mean(resid[wings_mask] **2))) if wings_mask.sum()  > 2 else rms_total
        ratio = rms_w / max(rms_c, 1e-12)

        # U-shape score: a planet-like (flat-bottomed, U-shaped) transit has
        # wing residuals ≈/slightly above centre residuals, giving ratio ≈ 1.0-1.3
        # on real folded data. A V-shaped EB leaves the fit unable to match the
        # sharp centre, pushing centre residuals up and ratio clearly below 1.
        # Empirically (controlled runs) genuine planets land in 1.0-1.4, so the
        # score must saturate near ~1.2 rather than the old "/1.5" ramp that
        # capped strong planets at ~0.6. Calibrated on synthetic TRUTH set:
        #   ratio = 1.20  →  S_τ = 1.00
        #   ratio = 0.85  →  S_τ ≈ 0.38   (V-ish, flagged below)
        #   ratio = 0.70  →  0.0
        s_limb = float(np.clip((ratio - 0.70) / 0.50, 0.0, 1.0))

        if ratio < 0.8:
            flags.append(f"V_SHAPE | ratio={ratio:.3f} < 0.8 → EB/systematic morphology")
        if rp_fit > 0.2:
            flags.append(f"LARGE_PLANET | R_p/R★={rp_fit:.4f} > 0.2 → potential EB")

        proof = (
            f"LIMB_SHAPE | batman fit: rp/rs={rp_fit:.5f}, a/rs={a_fit:.2f}, "
            f"inc={inc_fit:.2f}° | RMS_total={rms_total:.5f} | "
            f"RMS_centre={rms_c:.5f}, RMS_wings={rms_w:.5f} | "
            f"ratio=wings/centre={ratio:.3f} | "
            f"S_τ={s_limb:.4f}"
        )

        return LimbShapeResult(
            rp_rs=rp_fit, a_rs=a_fit, inclination_deg=inc_fit,
            u1=u1_fit, u2=u2_fit,
            residual_rms=rms_total,
            residual_centre=rms_c, residual_wings=rms_w,
            shape_ratio=ratio, s_limb=s_limb, proof=proof, flags=flags,
        )

    # ── Audit 4: Ingress/Egress ──────────────────────────────────────────────

    def ingress_egress_test(
        self,
        bin_phase:    np.ndarray,
        bin_flux:     np.ndarray,
        period:       float,
        duration:     float,
        transit_depth: float,
    ) -> IngressEgressResult:
        """Run the V-shape vs U-shape ingress/egress ratio test."""
        return self._ie_test.test(bin_phase, bin_flux, period, duration, transit_depth)

    # ── Audit 5: MCMC ────────────────────────────────────────────────────────

    def mcmc_validate(
        self,
        bin_phase:    np.ndarray,
        bin_flux:     np.ndarray,
        period:       float,
        duration:     float,
        transit_depth: float,
        n_steps:      int = 0,
    ) -> MCMCResult:
        """Run MCMC posterior validation.
        
        Parameters
        ----------
        n_steps : int, optional
            Override default MCMC step count. If 0, uses module default (500).
            Use lower values (50-100) for fast-reject on obvious false positives.
        """
        if not self.run_mcmc:
            return MCMCResult(
                n_walkers=0, n_steps=0, acceptance_fraction=0.0,
                params_median={}, params_stddev={},
                is_gaussian=True, skewness_max=0.0, kurtosis_max=0.0,
                gelman_rubin_max=1.0,
                proof="MCMC | disabled in config -> SKIPPED",
                flags=["MCMC_DISABLED"],
            )
        return self._mcmc_validator.validate(bin_phase, bin_flux, period, duration, transit_depth, n_steps=n_steps)

    # ── Trapezoidal fallback (no batman) ─────────────────────────────────────

    @staticmethod
    def _trapezoid_fallback(
        bin_phase:    np.ndarray,
        bin_flux:     np.ndarray,
        period:       float,
        duration:     float,
        transit_depth: float,
    ) -> LimbShapeResult:
        """
        Pure-numpy trapezoidal transit model.
        Used when batman is not installed.

        Trapezoid parameters: [depth, ingress_fraction, flat_fraction]
        """
        valid = np.isfinite(bin_flux)
        ph    = bin_phase[valid]
        fl    = bin_flux[valid]
        flags: List[str] = ["BATMAN_UNAVAILABLE | using trapezoidal fallback"]

        half_dur = (duration / period) / 2.0

        def trapezoid_model(phase: np.ndarray, depth: float, ingress_f: float, flat_f: float) -> np.ndarray:
            model = np.ones_like(phase)
            ingress_w = ingress_f * half_dur
            flat_w    = flat_f    * half_dur
            for i, p in enumerate(phase):
                ap = abs(p)
                if ap >= half_dur:
                    model[i] = 1.0
                elif ap >= flat_w:
                    frac = (ap - flat_w) / max(ingress_w, 1e-9)
                    model[i] = 1.0 - depth * (1.0 - np.clip(frac, 0, 1))
                else:
                    model[i] = 1.0 - depth
            return model

        try:
            popt, _ = curve_fit(
                trapezoid_model, ph, fl,
                p0=[transit_depth, 0.2, 0.6],
                bounds=([0, 0.01, 0.01], [0.5, 1.0, 1.0]),
                maxfev=2000,
            )
            model_fl = trapezoid_model(ph, *popt)
            resid    = fl - model_fl
            rms      = float(np.sqrt(np.mean(resid**2)))
            depth_fit, ingress_fit, flat_fit = popt

            # V vs U: flat fraction — high flat fraction → more U-shaped
            s_limb = float(np.clip(flat_fit, 0.0, 1.0))

            proof = (
                f"TRAPEZOID_FIT | depth={depth_fit:.5f}, ingress_frac={ingress_fit:.3f}, "
                f"flat_frac={flat_fit:.3f} | RMS={rms:.5f} | S_τ={s_limb:.4f}"
            )
        except Exception as exc:
            s_limb = 0.50
            proof  = f"TRAPEZOID_FIT | FAILED ({exc}) → S_τ=0.50"

        return LimbShapeResult(
            rp_rs=np.sqrt(transit_depth), a_rs=0, inclination_deg=90,
            u1=0.3, u2=0.1,
            residual_rms=0, residual_centre=0, residual_wings=0,
            shape_ratio=1.0, s_limb=s_limb, proof=proof, flags=flags,
        )
