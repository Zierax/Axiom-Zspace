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
_BIN_MAKE = _ROOT / "C99-Version" / "bin" / "zspace_card"
_BIN_MAKE_EXE = _ROOT / "C99-Version" / "bin" / "zspace_card.exe"
_BIN = _ROOT / "C99-Version" / "build" / "zspace_card"
_BIN_C99 = _ROOT / "C99-Version" / "build" / "zspace_card.exe"
# Derived WSL paths (computed dynamically from _ROOT, not hardcoded)
def _path_to_wsl_posix(p: Path | str) -> str:
    """Convert a Windows or POSIX path to WSL POSIX; handles any drive letter A-Z."""
    s = str(p).replace("\\", "/")
    # Already WSL absolute
    if s.startswith("/mnt/"):
        return s
    if s.startswith("/"):
        return s
    # Windows drive letter e.g. D:/...
    if len(s) >= 2 and s[1] == ":" and s[0].isalpha():
        drive = s[0].lower()
        rest = s[2:]
        if not rest.startswith("/"):
            rest = "/" + rest
        return f"/mnt/{drive}{rest}"
    return s

_LINUX_BIN = _path_to_wsl_posix(_BIN)
_WSL_BUILD_DIR = _path_to_wsl_posix(_BIN.parent)
_WSL_CD = f"cd {_WSL_BUILD_DIR} &&"

def _is_wsl() -> bool:
    try:
        return os.name == "posix" and Path("/proc/version").exists() and "microsoft" in Path("/proc/version").read_text(errors="ignore").lower()
    except OSError:
        return False

def _wsl_prefix() -> str:
    return "" if _is_wsl() else f"{_WSL_CD} "


def _find_binary() -> Optional[str]:
    for cand in (_BIN_MAKE, _BIN_MAKE_EXE, _BIN, _BIN_C99):
        if cand.exists():
            return str(cand)
    return None


def _to_wsl(p: str) -> str:
    # Prefer wslpath when available (handles any mount, spaces, etc.)
    try:
        # Use wsl wslpath -a if not inside WSL and wsl is available
        if not _is_wsl() and os.name != "posix":
            proc = subprocess.run(["wsl", "wslpath", "-a", str(p)], capture_output=True, text=True, timeout=5)
            if proc.returncode == 0 and proc.stdout.strip():
                return proc.stdout.strip()
    except (OSError, subprocess.SubprocessError, subprocess.TimeoutExpired):
        pass
    return _path_to_wsl_posix(p)


def _wsl_bin_posix() -> str:
    """Return the POSIX path of the C99 binary as seen inside WSL."""
    # Prefer bin/ over build/ when both exist inside WSL
    for cand in (_BIN_MAKE, _BIN):
        wsl = _path_to_wsl_posix(cand)
        # Quick existence check inside WSL is not possible here; return first candidate
        # that exists on Windows side
        if cand.exists():
            return wsl
    return _LINUX_BIN

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
        wsl_bin = _wsl_bin_posix()
        args = [wsl_bin, cand_path]
        if lc_path:
            args.append(lc_path)
        proc = subprocess.run(args, capture_output=True, text=True, timeout=600)
    else:
        cmd = f"{_WSL_CD} {_wsl_bin_posix()} {quote(_to_wsl(cand_path))}"
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
            [_wsl_bin_posix(),
             "bls", lc_path, str(period_min), str(period_max)],
            capture_output=True, text=True, timeout=600)
    else:
        quote = shlex.quote
        cmd = (f"{_WSL_CD} {_wsl_bin_posix()} bls "
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
            [_wsl_bin_posix()] + args,
            capture_output=True, text=True, timeout=600)
    else:
        quote = shlex.quote
        cmd = (f"{_WSL_CD} {_wsl_bin_posix()} "
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
            [_wsl_bin_posix()] + args,
            capture_output=True, text=True, timeout=600)
    else:
        quote = shlex.quote
        cmd = (f"{_WSL_CD} {_wsl_bin_posix()} flatten "
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