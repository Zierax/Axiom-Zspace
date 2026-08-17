"""
config.py · Central configuration loader
==========================================
Loads detection/validation thresholds from the production YAML so the
same gates are tunable at the configuration level (no code edits).

Priority order for locating the config file:
  1. $AXIOM_CONFIG (explicit env override)
  2. <repo root>/config/production.yaml
  3. <repo root>/Axiom-Zspace/config.yaml
  4. built-in defaults (fall back gracefully)
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Optional

import yaml

# ── Built-in defaults (used only when no config file is found) ─────────────
DEFAULTS: Dict[str, Any] = {
    "detection": {
        "bls_snr_threshold": 5.5,
        "fap_threshold": 5.0e-2,     # 5% FAP gate (config-tunable)
        "cvs_planet_threshold": 0.80,
        "snr_ref": 15.0,             # SNR at which S_P saturates to 1.0
    },
}

_CACHE: Optional[Dict[str, Any]] = None
_CONFIG_PATH: Optional[str] = None


def _candidate_paths() -> list:
    explicit = os.environ.get("AXIOM_CONFIG")
    here = Path(__file__).resolve().parent  # zspace_engine/
    repo = here.parent
    return [
        Path(explicit) if explicit else None,
        repo / "config" / "production.yaml",
        repo / "Axiom-Zspace" / "config.yaml",
    ]


def find_config() -> Optional[str]:
    """Locate the first existing config file, or None."""
    for p in _candidate_paths():
        if p is not None and p.is_file():
            return str(p)
    return None


def load_config(force_reload: bool = False) -> Dict[str, Any]:
    """Load + cache the merged config dict (file over defaults)."""
    global _CACHE, _CONFIG_PATH
    if _CACHE is not None and not force_reload:
        return _CACHE

    merged: Dict[str, Any] = {"detection": dict(DEFAULTS["detection"])}
    path = find_config()
    if path is not None:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                raw = yaml.safe_load(fh) or {}
            det = raw.get("detection") or {}
            for k, v in det.items():
                if v is not None:
                    merged["detection"][k] = v
            _CONFIG_PATH = path
        except Exception:
            _CONFIG_PATH = None
    else:
        _CONFIG_PATH = None

    _CACHE = merged
    return merged


def get_config() -> Dict[str, Any]:
    return load_config()


def config_path() -> Optional[str]:
    load_config()
    return _CONFIG_PATH


# ── Convenience getters (used by detectors / validator / benchmarks) ──────

def fap_threshold() -> float:
    return float(load_config()["detection"].get("fap_threshold", 5.0e-2))


def snr_threshold() -> float:
    return float(load_config()["detection"].get("bls_snr_threshold", 5.5))


def snr_ref() -> float:
    return float(load_config()["detection"].get("snr_ref", 15.0))


def cvs_planet_threshold() -> float:
    return float(load_config()["detection"].get("cvs_planet_threshold", 0.80))
