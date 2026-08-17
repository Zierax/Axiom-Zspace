"""
╔══════════════════════════════════════════════════════════════════════╗
║          PROJECT AXIOM-ZSPACE  ·  Deterministic Exoplanet Engine     ║
║          Version 1.0  ·  "Truthimatics" Framework                    ║
╚══════════════════════════════════════════════════════════════════════╝

White-box, physics-first exoplanet detection pipeline.
Every inference is traceable to a physical law or mathematical identity.
"""

from .core import VitalityMatrix, CompositeVitalityScore
from .ingestion import LightCurveIngester
from .detectors import BLSDetector, FAPValidator
from .auditors import TransitAuditor
from .context import StellarContextAuditor
from .report import TruthimaticsReport

__version__ = "1.0.0"
__framework__ = "Truthimatics"
__author__ = "Axiom-ZSpace Constructor"

__all__ = [
    "VitalityMatrix",
    "CompositeVitalityScore",
    "LightCurveIngester",
    "BLSDetector",
    "FAPValidator",
    "TransitAuditor",
    "StellarContextAuditor",
    "TruthimaticsReport",
]
