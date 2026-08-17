"""
ephemeris.py  ·  Harmonics / Alias Ephemeris Resolution
========================================================
Discovers whether a BLS-reported period P_best is the FUNDAMENTAL orbital
period or an integer sub-harmonic alias (P_true/2, P_true/3) of the same
physical signal, then yields the resolved candidate for re-validation.

WHY THIS EXISTS
---------------
A periodic transit detected at a sub-harmonic (e.g. P_true/3) folds every
event into ONE dip in its own fold — so it looks, to every single-fold test,
exactly like a genuine single-transit planet and certifies at the WRONG
ephemeris. The distinguishing physics appears in the INTEGER-MULTIPLE folds:

    fold at  m*P_best  for m = 2, 3   (when inside the search window)

  * FUNDAMENTAL (P_best == P_true) :
        2P fold → 2 equal dips (0, 0.5)
        3P fold → 3 equal dips (0, 1/3, 2/3)
  * P_TRUE/2 alias                 :
        own fold → 1 dip
        2P fold → 1 dip (restores P_true)
        3P fold → 2 dips
  * P_TRUE/3 alias                 :
        own fold → 1 dip
        3P fold → 1 dip (restores P_true)
  * Detached / grazing EB folded at half the true period :
        own fold → 2 DIPs OF UNEQUAL DEPTH (primary + secondary) →
        NOT a single-dip candidate → this resolver ignores it and the
        FP-5/FP-5b/FP-5c gates handle it.

RESOLUTION POLICY (fail-safe, zero-regression)
----------------------------------------------
Only a candidate that (a) has a single significant dip in its OWN fold, and
(b) already passed the full sovereign validator at P_best, may be resolved to
an integer multiple. The resolved period is then RE-VALIDATED by the complete
validator from scratch; the resolved ephemeris is adopted only if that second
validation certifies. Any ambiguity, noise-dominated fold, or out-of-window
multiple aborts back to the reported period.

Classifiers (own fold counted first; dips counted in the 2P and 3P folds)
------------------------------------------------------------------------
Folding any signal of true period P_t at a multiple m*P (P = reported period)
gives, for the three interesting cases:

    case         2P-fold dips (N2)    3P-fold dips (N3)
    FUNDAMENTAL       2               3
    P_TRUE/2 alias    1               3
    P_TRUE/3 alias    2               1

so the decision rule (own == 1 dip) is:
    N3 == 1                          -> P_TRUE/3  (physical = 3*P)
    N3 == 3 and N2 == 1              -> P_TRUE/2  (physical = 2*P)
    N3 == 3 and N2 == 2              -> FUNDAMENTAL (P_best is physical)
    anything else / any own-fold>1   -> FUNDAMENTAL (no resolution)

Note: a detached/grazing EB detected at HALF its true period folds primary
and secondary onto the same phase 0 (own fold = 1 dip — the same trap FP-5c
defends) but its 2P fold restores the EB's unequal primary/secondary pair,
so N2=2 with unequal depths; the resolver returns FUNDAMENTAL and the
FP-5/FP-5b/FP-5c gates reject it.

Over-harmonics (P_best == m*P_true, found LONGER than truth, e.g. true*3)
are deliberately left unresolvable in this stage: they present m-equal dips
in the own fold and resolving downward requires a divide branch that risks
re-certifying a fundamental planet at a fractional period. Recorded as an
honest wrong-ephemeris detection (known limitation).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# Fold dip detection (phase-folded, binned)
# ─────────────────────────────────────────────────────────────────────────────

# A dip counts as a real event only when BOTH hold:
#   * depth >= DIP_MIN_SNR × per-bin noise (matches a >3σ transit image), and
#   * depth >= REL_MIN_FRAC × deepest group in the SAME fold (kills sub-percent
#     noise nicks that would inflate the multiplicity and break the N2/N3 test).
DIP_MIN_SNR   = 5.0
REL_MIN_FRAC  = 0.25

@dataclass
class FoldDip:
    phase_center:  float              # weighted phase of the dip centre in [-0.5,0.5]
    phase_lo:      float              # phase of the first dip bin
    phase_hi:      float              # phase of the last dip bin
    depth:         float              # baseline - mean(binned flux inside group) (>=0)
    depth_snr:     float              # matched sig = depth_integrated / per-bin noise
    width_phase:   float              # (phase_hi - phase_lo) + bin width


@dataclass
class FoldSignature:
    period:          float
    n_bins:          int
    dips:            List[FoldDip]
    baseline:        float
    per_bin_noise:   float
    covered_fraction: float            # fraction of phase bins with enough cadences


@dataclass
class ResolvedEphemeris:
    """
    Resolution decision for one BLS candidate.
    multiple == 1  → candidate was already fundamental (no override).
    multiple  n    → physical period is n*period_best (n = 2 or 3).
    """
    candidate_period:  float
    physical_period:   float
    multiple:          int                       # 1, 2 or 3
    classifier:        str
    pattern:           dict                      # {'own':n1,'at_3p':N3,'at_2p':N2}
    confidence:        float                     # 0..1 of the deciding fold
    evidence:          str
    flags:             List[str] = field(default_factory=list)


def _fold_dip_signature(
    time:    np.ndarray,
    flux:    np.ndarray,
    period:  float,
    t0:      float,
    n_bins:  int = 200,
) -> FoldSignature:
    """Phase-fold `time/flux` at `period`, bin, and return the recognised dips.

    Two-pass dip counting:
      1. candidate groups: bins > ABS_SNR_MAD (3 MAD) below the robust median;
      2. confirmation: a group counts as a real event only if its depth is both
         >= DIP_MIN_SNR per-bin significance AND >= REL_MIN_FRAC of the deepest
         group in the SAME fold. The relative floor kills sub-percent noise
         nicks (which otherwise inflate the multiplicity and break the N2/N3
         classifier) while preserving EB secondaries (~0.5x primary) and the
         equal-depth images of a periodic signal.
    """
    from zspace_engine.detectors import BLSDetector

    ph, bf, berr = BLSDetector.fold_and_bin(time, flux, period, t0, n_bins=n_bins)
    valid = ~np.isnan(bf)
    covered = float(np.mean(valid))

    baseline = float(np.nanmedian(bf))
    dev = bf - baseline
    per_bin = 1.4826 * float(np.nanmedian(np.abs(dev)))
    if not (np.isfinite(per_bin) and per_bin > 1e-15):
        per_bin = float(np.nanstd(dev)) if np.isfinite(np.nanstd(dev)) else 1.0
        per_bin = max(per_bin, 1e-15)

    # Pass 1: candidate groups at ABS_SNR_MAD MAD below robust median.
    abs_snr_mad = 3.0
    thr = baseline - abs_snr_mad * per_bin
    below = np.flatnonzero(valid & (bf < thr))
    bin_w = 1.0 / n_bins

    if below.size == 0:
        return FoldSignature(period, n_bins, [], baseline, per_bin, covered)

    # Cluster contiguous/near-contiguous dip bins (gap <= 2 bins).
    groups: List[List[int]] = [[int(below[0])]]
    for j in below[1:]:
        if j - groups[-1][-1] <= 2:
            groups[-1].append(int(j))
        else:
            groups.append([int(j)])

    # Expand each 3-MAD core out to the full transit wings: any adjacent bin
    # still >= INCLUDE_SNR below baseline belongs to the same physical event.
    # Without this, a wide (alias-diluted) transit contributes only its 1-2
    # deepest bins to the group and the matched SNR undercounts the signal.
    include_snr = 1.0
    include_thr = baseline - include_snr * per_bin
    expanded: List[List[int]] = []
    for grp in groups:
        g = set(grp)
        grown = True
        while grown:
            grown = False
            for b in list(g):
                for nb in (b - 1, b + 1):
                    nb_w = nb % n_bins
                    if nb_w in g:
                        continue
                    if not valid[nb_w]:
                        continue
                    if bf[nb_w] <= include_thr:
                        g.add(nb_w)
                        grown = True
        expanded.append(sorted(g))
    groups = expanded

    cands: List[FoldDip] = []
    for grp in groups:
        g = np.asarray(grp, dtype=int)
        gphase = ph[g]
        gflux  = bf[g]
        weights = np.maximum((baseline - gflux), 0.0) + 1e-12
        center = float(np.average(gphase, weights=weights))
        depth  = float(np.mean(baseline - gflux))
        depth  = max(depth, 0.0)
        lo, hi = float(gphase.min()) - 0.5 * bin_w, float(gphase.max()) + 0.5 * bin_w
        # Matched significance: a real transit has the SAME integrated deficit
        # whatever phase width it spreads over (aliases put it in more bins).
        # depth × sqrt(n_bins) / per-bin-noise recovers the fundamental SNR
        # even when per-bin depth alone halves at an N*P_fold alias.
        n_bins_group = float(len(g))
        depth_snr = depth * math.sqrt(n_bins_group) / max(per_bin, 1e-15)
        cands.append(FoldDip(
            phase_center=center, phase_lo=lo, phase_hi=hi,
            depth=depth, depth_snr=depth_snr, width_phase=hi - lo,
        ))

    # Pass 2: relative-depth + significance confirmation against the deepest
    # group in this fold.
    max_depth = max((c.depth for c in cands), default=0.0)
    dips: List[FoldDip] = []
    for c in cands:
        if c.depth_snr < DIP_MIN_SNR:
            continue
        if max_depth > 1e-15 and c.depth < REL_MIN_FRAC * max_depth:
            continue
        dips.append(c)

    dips.sort(key=lambda d: d.phase_center)
    return FoldSignature(period, n_bins, dips, baseline, per_bin, covered)


def _dips_equal_depth(dips: List[FoldDip], max_ratio: float = 1.7) -> bool:
    """True when all dip depths agree within `max_ratio` (planet-like transit)."""
    if not dips:
        return False
    depths = np.array([max(d.depth, 0.0) for d in dips])
    dmax, dmin = float(depths.max()), float(depths.min())
    return dmin > 1e-12 and (dmax / dmin) <= max_ratio


def _min_dip_snr(dips: List[FoldDip]) -> float:
    return float(min((d.depth_snr for d in dips), default=0.0))


def _spacing_presence(dips: List[FoldDip], expected_spacing: float) -> bool:
    """Check that the dip centres lie approximately at multiples of spacing."""
    if len(dips) < 2:
        return True
    centers = np.sort([d.phase_center for d in dips])
    tol = 0.5 * expected_spacing
    for i in range(len(centers) - 1):
        # allow wrapping 0.5 -> -0.5
        d = abs(centers[i + 1] - centers[i])
        d = min(d, 1.0 - d)
        if abs(d - expected_spacing) > tol:
            return False
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Resolver
# ─────────────────────────────────────────────────────────────────────────────

class EphemerisResolver:
    """
    Classifies a single already-validated BLS candidate against its 2P/3P
    folds and returns the physical (fundamental) ephemeris when the candidate
    is an integer sub-harmonic alias.
    """

    DECIDE_MIN_SNR    = 4.0     # deciding multiple-fold must be this strong
    EQUAL_DEPTH_RATIO = 1.7     # max depth spread for "equal dips"
    PERIOD_TOL        = 0.05    # resolved BLS must land within 5% of expected

    def resolve(
        self,
        time:   np.ndarray,
        flux:   np.ndarray,
        period_best: float,
        t0:     float,
        period_min:  float = 0.5,
        period_max:  float = 13.5,
        transit_duration: Optional[float] = None,
    ) -> ResolvedEphemeris:
        """
        Return the physical ephemeris for a candidate reported at
        `period_best`. `transit_duration` is unused but accepted for API
        symmetry with the detector (kept for future scale-aware folds).
        """
        sig_own = _fold_dip_signature(time, flux, period_best, t0)
        n1 = len(sig_own.dips)
        pattern = {"own": n1}

        # OVER-harmonic probe (independent of the raw own-fold count): BLS can
        # lock a MULTIPLE of the physical period (P_best = n*P_true); folding at
        # P_best then shows n real, equally-spaced, equally-deep images. The raw
        # BLS t0 + coarse bins can hide or split those images, so re-fold with a
        # centred t0 and fine bins before deciding. A real multiplicity must be
        # spaced ≈ 1/n (a single transit split by binning sits ~1 transit-width
        # apart, which fails the spacing check).
        t0_fine = self._center_t0(time, flux, period_best, t0)
        sig_fine = _fold_dip_signature(time, flux, period_best, t0_fine, n_bins=600)
        fine_dips = self._merge_near_dips(sig_fine.dips)
        if 1 < len(fine_dips) <= 3:
            down = self._try_down_resolve(
                time, flux, period_best, t0, len(fine_dips), fine_dips, pattern, period_min,
            )
            if down is not None:
                return down

        if n1 != 1:
            return ResolvedEphemeris(
                candidate_period=period_best,
                physical_period=period_best,
                multiple=1,
                classifier="FUNDAMENTAL",
                pattern=pattern,
                confidence=0.0,
                evidence=(
                    f"own fold has {n1} dips (need exactly 1 to test sub-harmonics) "
                    f"→ no resolution"
                ),
            )

        # Count dips in the 2P and 3P folds (physical multiples inside window).
        n2_sig = None
        pat2 = None
        if 2.0 * period_best <= period_max * 1.02:
            n2_sig = _fold_dip_signature(time, flux, 2.0 * period_best, t0)
            pat2 = len(n2_sig.dips)
        pattern["at_2p"] = pat2

        n3_sig = None
        pat3 = None
        if 3.0 * period_best <= period_max * 1.02:
            n3_sig = _fold_dip_signature(time, flux, 3.0 * period_best, t0)
            pat3 = len(n3_sig.dips)
        pattern["at_3p"] = pat3

        # Classification table (own == 1 dip):
        #   FUNDAMENTAL      : N2=2, N3=3
        #   P_TRUE/2 alias   : N2=1, N3=3
        #   P_TRUE/3 alias   : N2=2, N3=1
        if pat3 is None:
            return self._resolved(
                period_best, 1.0, "FUNDAMENTAL", pattern, 0.0,
                "3P outside search window → cannot verify sub-harmonic → fundamental",
            )

        # ── P_TRUE/3 alias: N3 == 1 ──────────────────────────────────────────
        if pat3 == 1:
            conf = _min_dip_snr(n3_sig.dips) if n3_sig else 0.0
            if conf < self.DECIDE_MIN_SNR:
                return self._weak_pattern(pattern, "N3=1", conf, period_best)
            return self._resolved(
                period_best, 3.0, "P_TRUE/3_ALIAS", pattern, conf,
                f"own fold 1 dip; 3P fold restores 1 dip (SNR={conf:.1f}) → P_true = 3*P",
            )

        if pat3 == 3:
            # Distinguish FUNDAMENTAL (N2=2) from P_TRUE/2 alias (N2=1).
            if pat2 == 1:
                conf = _min_dip_snr(n2_sig.dips) if n2_sig else 0.0
                if conf < self.DECIDE_MIN_SNR:
                    return self._weak_pattern(pattern, "N3=3,N2=1", conf, period_best)
                return self._resolved(
                    period_best, 2.0, "P_TRUE/2_ALIAS", pattern, conf,
                    f"own fold 1 dip; 3P fold shows 3 dips, 2P fold restores 1 dip "
                    f"(SNR={conf:.1f}) → P_true = 2*P",
                )
            equal2 = _dips_equal_depth(n2_sig.dips) if n2_sig else False
            equal3 = _dips_equal_depth(n3_sig.dips) and _spacing_presence(n3_sig.dips, 1.0 / 3.0)
            return ResolvedEphemeris(
                candidate_period=period_best,
                physical_period=period_best,
                multiple=1,
                classifier="FUNDAMENTAL",
                pattern=pattern,
                confidence=_min_dip_snr(n3_sig.dips) if n3_sig else 0.0,
                evidence=(
                    f"own fold 1 dip; 3P fold shows 3 dip(s) "
                    f"{'equal+spaced' if equal3 else 'irregular'}; "
                    f"2P fold shows {pat2} dip(s) {'equal' if equal2 else 'unequal'} "
                    f"→ P_best is fundamental"
                ),
            )

        # Ambiguous multiple count (0, 2, or >3 at 3P) — noise dominated.
        return self._weak_pattern(pattern, f"N3={pat3}", 0.0, period_best)

    # ── helpers ──────────────────────────────────────────────────────────────

    @staticmethod
    def _center_t0(
        time:    np.ndarray,
        flux:    np.ndarray,
        period:  float,
        t0:      float,
        n_bins:  int = 400,
    ) -> float:
        """Shift `t0` so the deepest bin of the fold sits at phase 0.

        Bin-edge smearing with an arbitrary reference epoch hides shallow
        transit images (alias/harmonic folds); centering the fold on the
        deepest event makes the dip multiplicity robust.
        """
        from zspace_engine.detectors import BLSDetector
        ph, bf, _ = BLSDetector.fold_and_bin(time, flux, period, t0, n_bins=n_bins)
        if not np.any(np.isfinite(bf)):
            return t0
        i = int(np.nanargmin(np.where(np.isfinite(bf), bf, np.inf)))
        return t0 + float(ph[i]) * period

    @staticmethod
    def _merge_near_dips(dips: List[FoldDip], tol: float = 0.02) -> List[FoldDip]:
        """Merge dips whose centres are within `tol` phase (a single inline
        transit split across a bin edge produces duplicate near-identical
        events). Returns a new sorted list."""
        centers = sorted([d.phase_center for d in dips])
        merged: List[FoldDip] = []
        for d in dips:
            if not merged:
                merged.append(d)
                continue
            last = merged[-1]
            sep = min(abs(d.phase_center - last.phase_center),
                      1.0 - abs(d.phase_center - last.phase_center))
            if sep < tol:
                # keep the deeper of the two near-identical events
                if d.depth > last.depth:
                    merged[-1] = d
            else:
                merged.append(d)
        merged.sort(key=lambda x: x.phase_center)
        return merged

    def _try_down_resolve(
        self,
        time:   np.ndarray,
        flux:   np.ndarray,
        period_best: float,
        t0:     float,
        n1:     int,
        own_dips: List[FoldDip],
        pattern: dict,
        period_min: float,
    ) -> Optional[ResolvedEphemeris]:
        """
        OVER-harmonic branch: BLS locked P_best = n1 * P_true and the own fold
        therefore shows n1 regular, equally-deep dips. Confirm by folding DOWN
        at P_best/n1 and requiring exactly one significant dip to survive.
        """
        # Only 2P and 3P over-harmonics (n1 == 2 or 3) are physically sane.
        if n1 not in (2, 3):
            return None
        if not _dips_equal_depth(own_dips):
            return None
        if not _spacing_presence(own_dips, 1.0 / n1):
            return None

        p_down = period_best / float(n1)
        if p_down < period_min * 0.98:
            return None

        t0_down = self._center_t0(time, flux, p_down, t0)
        sig_down = _fold_dip_signature(time, flux, p_down, t0_down, n_bins=400)
        nd = len(self._merge_near_dips(sig_down.dips))
        pattern[f"at_down_{n1}p"] = nd
        if nd != 1:
            return None

        conf = _min_dip_snr(sig_down.dips)
        if conf < self.DECIDE_MIN_SNR:
            return None

        return ResolvedEphemeris(
            candidate_period=period_best,
            physical_period=p_down,
            multiple=int(n1),
            classifier="OVER_HARMONIC",
            pattern=pattern,
            confidence=conf,
            evidence=(
                f"own fold {n1} equally-spaced dips; P/{n1} fold restores 1 dip "
                f"(SNR={conf:.1f}) → P_true = P/{n1}"
            ),
            flags=["EPHEMERIS_RESOLVED", f"RESOLVE_DOWN_{n1}"],
        )

    def _weak_pattern(self, pattern: dict, pat3, conf: float, period_best: float) -> ResolvedEphemeris:
        return ResolvedEphemeris(
            candidate_period=period_best,
            physical_period=period_best,
            multiple=1,
            classifier="FUNDAMENTAL",
            pattern=pattern,
            confidence=conf,
            evidence=(
                f"3P fold weak/ambiguous ({pat3}, min dip SNR={conf:.1f} < "
                f"{self.DECIDE_MIN_SNR}) → no resolution (keep P)"
            ),
            flags=["LOW_CONFIDENCE_MULTIPLE_FOLD"],
        )

    def _resolved(
        self,
        period_best: float,
        mult: float,
        classifier: str,
        pattern: dict,
        conf: float,
        evidence: str,
    ) -> ResolvedEphemeris:
        return ResolvedEphemeris(
            candidate_period=period_best,
            physical_period=mult * period_best,
            multiple=int(math.ceil(mult)),
            classifier=classifier,
            pattern=pattern,
            confidence=conf,
            evidence=evidence,
            flags=["EPHEMERIS_RESOLVED", f"RESOLVE_X{int(math.ceil(mult))}"],
        )