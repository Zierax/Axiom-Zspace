#!/usr/bin/env python3
"""
c99_bridge.py  ·  C99-Version Sovereign Card bridge for the Python pipeline.

Builds a candidate bundle + light-curve CSV for the C99 engine
(bin/zspace_card, built from C99-Version/) and returns the parsed
Sovereign Logic Card.  Run via WSL when the binary is ELF (default).

Usage (inside run_pipeline.py):
    from c99_bridge import run_c99_sovereign
    card = run_c99_sovereign(candidate, time, flux, ...)

Never imported by default; only when --engine c99 is requested.
"""
from __future__ import annotations

import csv
import json
import os
import shlex
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Sequence

_ROOT = Path(__file__).resolve().parent
_BIN = _ROOT / "C99-Version" / "build" / "zspace_card"
_BIN_C99 = _ROOT / "C99-Version" / "build" / "zspace_card.exe"
_LINUX_BIN = "/mnt/d/Axioms/Axiom-Zspace-CODE/C99-Version/build/zspace_card"
_WSL_CD = "cd /mnt/d/Axioms/Axiom-Zspace-CODE/C99-Version/build &&"

def _is_wsl() -> bool:
    return os.name == "posix" and Path("/proc/version").exists() and "microsoft" in Path("/proc/version").read_text(errors="ignore").lower()

def _wsl_prefix() -> str:
    return "" if _is_wsl() else f"{_WSL_CD} "


def _find_binary() -> Optional[str]:
    if _BIN.exists():
        return str(_BIN)
    if _BIN_C99.exists():
        return str(_BIN_C99)
    return None


def _to_wsl(p: str) -> str:
    p = p.replace("\\", "/")
    for drive, mount in (("D:", "/mnt/d"), ("C:", "/mnt/c"), ("E:", "/mnt/e")):
        if p.startswith(drive):
            return mount + p[len(drive):]
    return p


def _run(cand_path: str, lc_path: Optional[str]) -> dict:
    quote = shlex.quote
    bin_path = _find_binary()
    args = []
    if bin_path and bin_path.endswith(".exe"):
        args = [bin_path, cand_path]
        if lc_path:
            args.append(lc_path)
        proc = subprocess.run(args, capture_output=True, text=True, timeout=600)
    elif _is_wsl():
        args = ["/mnt/d/Axioms/Axiom-Zspace-CODE/C99-Version/build/zspace_card",
                cand_path]
        if lc_path:
            args.append(lc_path)
        proc = subprocess.run(args, capture_output=True, text=True, timeout=600)
    else:
        cmd = f"{_WSL_CD} {_LINUX_BIN} {quote(_to_wsl(cand_path))}"
        if lc_path:
            cmd += f" {quote(_to_wsl(lc_path))}"
        proc = subprocess.run(["wsl", "bash", "-lc", cmd],
                              capture_output=True, text=True, timeout=600)
    if proc.returncode != 0:
        raise RuntimeError(f"zspace_card failed (rc={proc.returncode}): "
                           f"{proc.stderr[:1500]}")
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as e:
        with open(os.path.join(tempfile.gettempdir(), "zspace_c99_stdout.txt"), "w") as fh:
            fh.write(proc.stdout)
        raise RuntimeError(f"zspace_card emitted invalid JSON ({e}); "
                           f"stdout dumped to %TEMP%\\zspace_c99_stdout.txt")


def run_c99_bls(
    time: Sequence[float],
    flux: Sequence[float],
    period_min: float = 0.5,
    period_max: float = 13.5,
    flux_err: Optional[Sequence[float]] = None,
) -> dict:
    """
    Run the C99 BLS periodogram engine (zspace_card bls) on a light curve.

    Returns a dict with the C99 keys:
      period_days power snr fap duration_hrs t0_days depth lc_points
    """
    tmp = tempfile.mkdtemp(prefix="zspace_c99_bls_")
    lc_path = os.path.join(tmp, "lightcurve.csv")
    with open(lc_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        if flux_err is not None and len(flux_err) == len(time):
            w.writerow(["time", "flux", "flux_err"])
            for t, f, e in zip(time, flux, flux_err):
                w.writerow([f"{t:.8f}", f"{f:.10f}", f"{e:.10f}"])
        else:
            w.writerow(["time", "flux"])
            for t, f in zip(time, flux):
                w.writerow([f"{t:.8f}", f"{f:.10f}"])

    bin_path = _find_binary()
    if bin_path and bin_path.endswith(".exe"):
        proc = subprocess.run(
            [bin_path, "bls", lc_path, str(period_min), str(period_max)],
            capture_output=True, text=True, timeout=600)
    elif _is_wsl():
        proc = subprocess.run(
            ["/mnt/d/Axioms/Axiom-Zspace-CODE/C99-Version/build/zspace_card",
             "bls", lc_path, str(period_min), str(period_max)],
            capture_output=True, text=True, timeout=600)
    else:
        quote = shlex.quote
        cmd = (f"{_WSL_CD} {_LINUX_BIN} bls "
               f"{quote(_to_wsl(lc_path))} {period_min} {period_max}")
        proc = subprocess.run(["wsl", "bash", "-lc", cmd],
                              capture_output=True, text=True, timeout=600)
    if proc.returncode != 0:
        raise RuntimeError(f"zspace_card bls failed (rc={proc.returncode}): "
                           f"{proc.stderr[:1500]}")
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as e:
        raise RuntimeError(f"zspace_card bls emitted invalid JSON ({e}): "
                           f"{proc.stdout[:1500]}")


def run_c99_audit(
    time: Sequence[float],
    flux: Sequence[float],
    period: float,
    t0: float,
    duration_hrs: float,
    transit_depth: float,
) -> dict:
    """
    Run the C99 transit audits (zspace_card audit) on a light curve.

    Returns dict with keys:
      even_odd.{n_even,n_odd,depth_even,depth_odd,depth_even_err,depth_odd_err,
                delta_sigma,t_stat,p_value,is_eb_flag}
      depth_consistency.{n_transits,mean_depth,std_depth,cv,sigma_med,chi2_red,s_depth}
      secondary_eclipse.{primary_depth,secondary_depth,secondary_ratio,
                         secondary_snr,n_primary,n_secondary,ok}
      ingress_egress.{depth_fit,ingress_fraction,flat_fraction,ingress_hrs,
                      flat_hrs,is_v_shape,fp_risk,fit_ok}
    """
    tmp = tempfile.mkdtemp(prefix="zspace_c99_audit_")
    lc_path = os.path.join(tmp, "lightcurve.csv")
    with open(lc_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time", "flux"])
        for t, f in zip(time, flux):
            w.writerow([f"{t:.8f}", f"{f:.10f}"])

    args = ["audit", lc_path,
            f"{period:.8f}", f"{t0:.8f}",
            f"{duration_hrs:.6f}", f"{transit_depth:.8f}"]
    bin_path = _find_binary()
    if bin_path and bin_path.endswith(".exe"):
        proc = subprocess.run([bin_path] + args,
                              capture_output=True, text=True, timeout=600)
    elif _is_wsl():
        proc = subprocess.run(
            ["/mnt/d/Axioms/Axiom-Zspace-CODE/C99-Version/build/zspace_card"] + args,
            capture_output=True, text=True, timeout=600)
    else:
        quote = shlex.quote
        cmd = (f"{_WSL_CD} {_LINUX_BIN} "
               + " ".join(quote(_to_wsl(a)) if a == lc_path else a
                          for a in args))
        proc = subprocess.run(["wsl", "bash", "-lc", cmd],
                              capture_output=True, text=True, timeout=600)
    if proc.returncode != 0:
        raise RuntimeError(f"zspace_card audit failed (rc={proc.returncode}): "
                           f"{proc.stderr[:1500]}")
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as e:
        raise RuntimeError(f"zspace_card audit emitted invalid JSON ({e}): "
                           f"{proc.stdout[:1500]}")


def run_c99_flatten(
    time: Sequence[float],
    flux: Sequence[float],
    period_days: float,
) -> Sequence[float]:
    """
    Run the C99 Savitzky-Golay flattening (zspace_card flatten) on a light curve.

    Matches LightCurveIngester._savgol_flatten (window = 0.75*P / cadence,
    min 51, odd, polyorder <= 3, scipy-compatible edges).  Returns the
    flattened flux in the original time order.
    """
    tmp = tempfile.mkdtemp(prefix="zspace_c99_flat_")
    in_path = os.path.join(tmp, "lightcurve.csv")
    out_path = os.path.join(tmp, "flattened.csv")
    with open(in_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time", "flux"])
        for t, f in zip(time, flux):
            w.writerow([f"{t:.8f}", f"{f:.10f}"])

    args = ["flatten", in_path, out_path, f"{period_days:.8f}"]
    bin_path = _find_binary()
    if bin_path and bin_path.endswith(".exe"):
        proc = subprocess.run([bin_path] + args,
                              capture_output=True, text=True, timeout=600)
    elif _is_wsl():
        proc = subprocess.run(
            ["/mnt/d/Axioms/Axiom-Zspace-CODE/C99-Version/build/zspace_card"] + args,
            capture_output=True, text=True, timeout=600)
    else:
        quote = shlex.quote
        cmd = (f"{_WSL_CD} {_LINUX_BIN} flatten "
               + " ".join(quote(_to_wsl(a)) if a in (in_path, out_path) else a
                          for a in args[1:]))
        proc = subprocess.run(["wsl", "bash", "-lc", cmd],
                              capture_output=True, text=True, timeout=600)
    if proc.returncode != 0:
        raise RuntimeError(f"zspace_card flatten failed (rc={proc.returncode}): "
                           f"{proc.stderr[:1500]}")
    try:
        with open(out_path, "r", encoding="utf-8") as fh:
            rows = [r for r in csv.reader(fh) if r and r[0] != "time"]
        return [float(r[1]) for r in rows]
    except (OSError, IndexError, ValueError) as e:
        raise RuntimeError(f"zspace_card flatten output unreadable ({e})")


def run_c99_sovereign(
    candidate: dict,
    time: Optional[Sequence[float]] = None,
    flux: Optional[Sequence[float]] = None,
) -> dict:
    """
    Execute the C99 Sovereign engine on a candidate bundle.

    candidate keys (subset of ZSCandidate, floats):
      period_days transit_depth transit_duration_hrs t0_days
      stellar_mass_solar stellar_radius_solar stellar_teff_k stellar_logg
      planet_radius_earth bls_snr bls_fap even_odd_delta_sigma shape_ratio
      secondary_snr secondary_depth_ratio alias_secondary_ratio
      coherent_evidence centroid_sigma limb_dark_u1 limb_dark_u2
      s_periodicity s_depth s_limb s_stellar

    Returns the full Sovereign Logic Card dict.
    """
    tmp = tempfile.mkdtemp(prefix="zspace_c99_")
    cand_path = os.path.join(tmp, "candidate.txt")
    with open(cand_path, "w", encoding="utf-8") as fh:
        for k, v in candidate.items():
            if v is None:
                continue
            if isinstance(v, bool):
                v = int(v)
            fh.write(f"{k}={v}\n")

    lc_path = None
    if time is not None and flux is not None and len(time) > 0:
        lc_path = os.path.join(tmp, "lightcurve.csv")
        with open(lc_path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["time", "flux"])
            for t, f in zip(time, flux):
                w.writerow([f"{t:.8f}", f"{f:.10f}"])

    return _run(cand_path, lc_path)