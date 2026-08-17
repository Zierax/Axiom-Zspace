"""
report.py  ·  The Truthimatics JSON Report Generator
======================================================
Assembles all pipeline results into a fully documented, machine-readable
Discovery Card conforming to the Truthimatics V3.0 specification.

Every numeric result is accompanied by its "Logic Proof" string.
No result is opaque; every number traces to a physical law.

Output schema (top-level keys)
------------------------------
  schema_version
  zspace_id
  tic_id
  planet_order
  timestamp_utc
  pipeline_version
  verdict
  composite_vitality_score    (CVS breakdown)
  orbital_mechanics
  bls_detection
  transit_audits
    even_odd
    depth_consistency
    limb_shape
  stellar_context
    metadata
    centroid_shift
    secondary_eclipse
  ingestion_audit
  all_flags
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

import numpy as np

from .core       import VitalityMatrix
from .ingestion  import LightCurveProduct
from .detectors  import BLSResult
from .auditors   import EvenOddResult, DepthConsistencyResult, LimbShapeResult
from .context    import ContextAuditResult
from .logging_config import get_logger

SCHEMA_VERSION   = "3.0"
PIPELINE_VERSION = "1.0.0"

# Module logger
logger = get_logger(__name__)


class TruthimaticsReport:
    """
    Assembles all pipeline artifacts into the Truthimatics Discovery Card.

    Usage
    -----
    report = TruthimaticsReport(vitality_matrix)
    report.attach_ingestion(lc_product)
    report.attach_bls(bls_result)
    report.attach_audits(eo, depth, limb)
    report.attach_context(context_result)
    card = report.emit()
    """

    def __init__(self, matrix: VitalityMatrix) -> None:
        self._matrix   = matrix
        self._lc:       Optional[LightCurveProduct]    = None
        self._bls:      Optional[BLSResult]            = None
        self._eo:       Optional[EvenOddResult]        = None
        self._depth:    Optional[DepthConsistencyResult] = None
        self._limb:     Optional[LimbShapeResult]      = None
        self._context:  Optional[ContextAuditResult]   = None

    def attach_ingestion(self, lc: LightCurveProduct) -> "TruthimaticsReport":
        self._lc = lc
        return self

    def attach_bls(self, bls: BLSResult) -> "TruthimaticsReport":
        self._bls = bls
        return self

    def attach_audits(
        self,
        eo:    EvenOddResult,
        depth: DepthConsistencyResult,
        limb:  LimbShapeResult,
    ) -> "TruthimaticsReport":
        self._eo    = eo
        self._depth = depth
        self._limb  = limb
        return self

    def attach_context(self, ctx: ContextAuditResult) -> "TruthimaticsReport":
        self._context = ctx
        return self

    def emit(self) -> Dict[str, Any]:
        """
        Finalise the VitalityMatrix and assemble the complete Discovery Card.
        Returns a dict that is JSON-serialisable.
        """
        core_dict = self._matrix.finalize()

        # Collect every flag from all subsystems
        all_flags: List[str] = []
        if self._bls:
            all_flags.extend(self._bls.flags)
        if self._eo:
            pass  # EB flag is in depth result
        if self._depth:
            all_flags.extend(self._depth.flags)
        if self._limb:
            all_flags.extend(self._limb.flags)
        if self._context:
            all_flags.extend(self._context.flags)

        card: Dict[str, Any] = {
            "schema_version":   SCHEMA_VERSION,
            "pipeline_version": PIPELINE_VERSION,
            "zspace_id":        core_dict["zspace_id"],
            "tic_id":           core_dict["tic_id"],
            "planet_order":     core_dict["planet_order"],
            "timestamp_utc":    datetime.now(timezone.utc).isoformat(),
            "verdict":          core_dict["cvs"]["verdict"],

            "composite_vitality_score": self._format_cvs(core_dict),
            "orbital_mechanics":        self._format_orbital(core_dict),
            "bls_detection":            self._format_bls(),
            "transit_audits": {
                "even_odd":          self._format_even_odd(),
                "depth_consistency": self._format_depth(),
                "limb_shape":        self._format_limb(),
            },
            "stellar_context":   self._format_context(),
            "ingestion_audit":   self._format_ingestion(),
            "all_flags":         all_flags,
            "truthimatics_seal": self._generate_seal(core_dict["cvs"]["cvs"]),
        }

        return card

    # ── Formatters ────────────────────────────────────────────────────────────

    def _format_cvs(self, core_dict: dict) -> dict:
        cvs_data = core_dict["cvs"]
        return {
            "value":     cvs_data["cvs"],
            "verdict":   cvs_data["verdict"],
            "threshold_map": {
                "PLANET_CANDIDATE":           "CVS ≥ 0.80",
                "LIKELY_PLANET_CANDIDATE":    "CVS ≥ 0.55",
                "AMBIGUOUS":                  "CVS ≥ 0.35",
                "FALSE_POSITIVE":             "CVS  < 0.35",
            },
            "proof": f"CVS = Σ(w·S)/Σ(w) = {cvs_data['cvs']:.6f}",
            "components": cvs_data["components"],
            "proof_chain": cvs_data["proof_chain"],
        }

    def _format_orbital(self, core_dict: dict) -> dict:
        om = core_dict["orbital_mechanics"]
        return {
            "period_days": {
                "value": om["period_days"],
                "proof": f"BLS best-fit period = {om['period_days']:.6f} d",
            },
            "semi_major_axis_au": {
                "value": om["semi_major_axis_au"],
                "proof": om["proof_semi_major"],
            },
            "equilibrium_temperature_k": {
                "value": om["equilibrium_temp_k"],
                "proof": om["proof_teq"],
            },
            "planet_radius_earth": {
                "value": om["planet_radius_earth"],
                "proof": om["proof_radius"],
            },
            "transit_depth_ppm": {
                "value": om["transit_depth_ppm"],
                "proof": f"transit_depth = {om['transit_depth_ppm']:.2f} ppm",
            },
            "albedo_assumed": om["albedo"],
        }

    def _format_bls(self) -> dict:
        if not self._bls:
            return {"status": "NOT_RUN"}
        b = self._bls
        return {
            "period_days":      b.period_best,
            "transit_depth":    b.transit_depth,
            "transit_duration_hours": round(b.transit_duration * 24, 4),
            "t0":               b.t0,
            "snr": {
                "value": round(b.snr, 4),
                "proof": f"{b.snr:.2f} > {b.snr_threshold} → {'PASS' if b.snr > b.snr_threshold else 'FAIL'}",
            },
            "fap": {
                "value": b.fap,
                "proof": f"FAP={b.fap:.3e} < {b.fap_threshold:.0e} → {'PASS' if b.fap < b.fap_threshold else 'FAIL'}",
            },
            "n_trial_periods":  b.n_trial_periods,
            "detection_gate":   "PASS" if b.passed_detection_gate() else "FAIL",
            "proof":            b.proof,
            "flags":            b.flags,
        }

    def _format_even_odd(self) -> dict:
        if not self._eo:
            return {"status": "NOT_RUN"}
        e = self._eo
        n_total = e.n_even + e.n_odd
        if n_total == 0:
            eb_proof = f"INSUFFICIENT_DATA — test skipped (total transits={n_total})"
        elif e.delta_sigma == 0.0 and n_total < 4:
            eb_proof = f"INSUFFICIENT_DATA — Δσ=0.00, only {n_total} transit(s) available"
        else:
            eb_proof = f"Δσ={e.delta_sigma:.2f} > 3.0 → {'EB_FLAG' if e.is_eb_flag else 'PASS'}"
        return {
            "n_even":          e.n_even,
            "n_odd":           e.n_odd,
            "depth_even":      round(e.depth_even, 7) if not np.isnan(e.depth_even) else None,
            "depth_odd":       round(e.depth_odd,  7) if not np.isnan(e.depth_odd)  else None,
            "delta_sigma":     round(e.delta_sigma, 4),
            "eb_flagged": {
                "value": e.is_eb_flag,
                "proof": eb_proof,
            },
            "proof": e.proof,
        }

    def _format_depth(self) -> dict:
        if not self._depth:
            return {"status": "NOT_RUN"}
        d = self._depth
        return {
            "n_transits":  len(d.depths),
            "mean_depth":  round(d.mean_depth, 7),
            "std_depth":   round(d.std_depth,  7),
            "cv": {
                "value": round(d.cv, 6),
                "proof": f"CV = σ/μ = {d.std_depth:.5f}/{d.mean_depth:.5f} = {d.cv:.6f}",
            },
            "s_depth": {
                "value": round(d.s_depth, 6),
                "proof": f"S_δ = max(0, 1-CV/0.10) = max(0, 1-{d.cv:.4f}/0.10) = {d.s_depth:.4f}",
            },
            "proof": d.proof,
            "flags": d.flags,
        }

    def _format_limb(self) -> dict:
        if not self._limb:
            return {"status": "NOT_RUN"}
        l = self._limb
        return {
            "model": "Mandel-Agol (batman)" if l.a_rs > 0 else "Trapezoidal",
            "rp_rs":              round(l.rp_rs, 6),
            "a_rs":               round(l.a_rs, 4),
            "inclination_deg":    round(l.inclination_deg, 4),
            "limb_darkening":     {"u1": round(l.u1, 4), "u2": round(l.u2, 4)},
            "residual_rms":       round(l.residual_rms, 7),
            "shape_analysis": {
                "residual_centre": round(l.residual_centre, 7),
                "residual_wings":  round(l.residual_wings, 7),
                "shape_ratio": {
                    "value": round(l.shape_ratio, 4),
                    "proof": (
                        f"ratio = wings_rms/centre_rms = "
                        f"{l.residual_wings:.5f}/{l.residual_centre:.5f} = {l.shape_ratio:.4f} | "
                        f"{'ratio>1 → U-shape → PLANET-LIKE' if l.shape_ratio > 1.05 else 'ratio≈1 → NEUTRAL SHAPE' if l.shape_ratio >= 0.95 else 'ratio<1 → V-shape → EB/ARTIFACT'}"
                    ),
                },
            },
            "s_limb": {
                "value": round(l.s_limb, 6),
                "proof": l.proof,
            },
            "flags": l.flags,
        }

    def _format_context(self) -> dict:
        if not self._context:
            return {"status": "NOT_RUN"}
        ctx = self._context
        m   = ctx.metadata
        c   = ctx.centroid
        s   = ctx.secondary

        result = {
            "stellar_metadata": {
                "source":                m.source,
                "tic_version":           getattr(m, 'tic_version', 'unknown'),
                "stellar_mass_solar":    m.stellar_mass_solar,
                "stellar_radius_solar":  m.stellar_radius_solar,
                "stellar_teff_k":        m.stellar_teff_k,
                "stellar_logg":          m.stellar_logg,
                "stellar_density_cgs":   getattr(m, 'stellar_density_cgs', 0.0),
                "contamination_ratio":   m.contamination_ratio,
                "gaia_parallax_mas":     m.gaia_parallax_mas,
            },
            "centroid_shift_test": {
                "shift_sigma":  round(c.centroid_shift_sigma, 4),
                "is_flagged":   c.is_flagged,
                "method":       getattr(c, 'method', 'flux_proxy'),
                "proof":        c.proof,
            },
            "secondary_eclipse_search": {
                "depth_at_phase_0.5":    round(s.depth_at_half_phase, 7),
                "snr_at_phase_0.5":      round(s.snr_at_half_phase, 4),
                "is_flagged":            s.is_flagged,
                "proof":                 s.proof,
            },
            "density_constraint": {
                "a_rs_transit":          round(ctx.density_check.a_rs_transit, 4),
                "a_rs_catalog":          round(ctx.density_check.a_rs_catalog, 4),
                "deviation_percent":     round(ctx.density_check.fractional_deviation * 100, 2),
                "is_flagged":            ctx.density_check.is_flagged,
                "proof":                 ctx.density_check.proof,
            },
            "multi_sector_consistency": {
                "n_sectors_available":   ctx.multi_sector.n_sectors_available,
                "n_sectors_consistent":  ctx.multi_sector.n_sectors_consistent,
                "sectors_checked":       ctx.multi_sector.sectors_checked,
                "confidence_boost":      ctx.multi_sector.confidence_boost,
                "is_consistent":         ctx.multi_sector.is_consistent,
                "proof":                 ctx.multi_sector.proof,
            },
            "s_stellar": {
                "value": round(ctx.s_stellar, 6),
                "proof": ctx.proof,
            },
            "flags": ctx.flags,
        }
        return result

    def _format_ingestion(self) -> dict:
        if not self._lc:
            return {"status": "NOT_RUN"}
        lc = self._lc
        return {
            "tic_id":              lc.tic_id,
            "sector":              lc.sector,
            "n_points_raw":        lc.n_points_raw,
            "n_points_cleaned":    lc.n_points_cleaned,
            "n_dropped_quality":   lc.n_dropped_quality,
            "n_dropped_sigma":     lc.n_dropped_sigma,
            "coverage_fraction": {
                "value": round(lc.coverage_fraction, 4),
                "proof": (
                    f"coverage = {lc.n_points_cleaned} / {lc.n_points_raw} = "
                    f"{lc.coverage_fraction:.4f}"
                ),
            },
            "cadence_days": round(lc.cadence_days, 6),
            "audit_log":          lc.audit_log,
        }

    @staticmethod
    def _generate_seal(cvs: float) -> str:
        """A short human-readable certification string."""
        if cvs is None:
            return "TRUTHIMATICS-SEAL | UNCOMPUTED"
        cvs_pct = round(cvs * 100, 2)
        return (
            f"TRUTHIMATICS-SEAL | CVS={cvs_pct}% | "
            f"WHITE-BOX | NO-ML | NO-IMPUTATION | "
            f"PHYSICS-DERIVED | Axiom-ZSpace v1.0"
        )

    # ── Serialisation helpers ─────────────────────────────────────────────────

    @staticmethod
    def to_json(card: dict, indent: int = 2) -> str:
        """JSON serialise, handling numpy types."""
        def default(obj: Any) -> Any:
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.bool_):
                return bool(obj)
            raise TypeError(f"Object of type {type(obj)} is not JSON serialisable")

        return json.dumps(card, indent=indent, default=default)

    def save(self, card: dict, filepath: str) -> None:
        """Write the Discovery Card to a JSON file with error handling."""
        try:
            with open(filepath, "w", encoding="utf-8") as fh:
                fh.write(self.to_json(card))
            logger.info(f"[ZSPACE] Discovery Card saved -> {filepath}")
        except IOError as e:
            logger.error(f"File I/O error writing Discovery Card to {filepath}: {e}")
            raise IOError(f"Failed to write Discovery Card to {filepath}: {e}") from e
        except Exception as e:
            logger.error(f"Unexpected error writing Discovery Card to {filepath}: {e}")
            raise RuntimeError(f"Failed to write Discovery Card to {filepath}: {e}") from e
