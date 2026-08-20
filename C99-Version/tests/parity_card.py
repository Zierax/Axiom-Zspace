#!/usr/bin/env python3
"""
parity_card.py — C99 vs Python parity harness for the Sovereign Logic Card.

Generates N random candidate profiles, runs ./bin/zspace_card for each,
computes the same quantities with the actual Python engine
(zspace_engine.validator.ProofEngine) and compares them field by field.

Usage:
    python parity_card.py [--n 60] [--lc] [--seed 20260817]
"""
import argparse
import json
import math
import os
import random
import re
import shlex
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(ROOT, "bin", "zspace_card")

sys.path.insert(0, os.path.abspath(os.path.join(ROOT, "..")))
from zspace_engine.validator import ProofEngine          # noqa: E402
import zspace_engine.validator as _V                     # noqa: E402

# No network in parity harness: FP-9 behaves as CATALOG_OFFLINE (PASS),
# which is exactly what the C engine implements.
_V.check_external_catalogs = lambda tic_id: {
    "is_multiple": False, "catalog_source": "Offline",
    "risk_level": "UNKNOWN", "classification": "offline",
    "query_latency_ms": 0.0,
}

from zspace_engine import thresholds as _T               # noqa: E402
_T = _T.threshold

FIELDS_NUM = [
    "a_au", "k", "delta_ld_corrected", "rp_earth",
    "rho_transit_gcc", "rho_tic_gcc", "density_ratio",
    "P_tr", "impact_parameter_b", "ingress_hrs", "i_min_deg",
]
FIELDS_INT = ["n_tests", "n_pass", "n_fail", "n_critical", "n_critical_pass", "n_transits"]
FIELDS_STR = ["sovereign_verdict", "overall_verdict", "cvs_verdict"]


def make_candidate(rng, lc=False):
    M = round(rng.uniform(0.3, 2.0), 5)
    P = round(rng.uniform(0.6, 30.0), 5)
    depth = round(rng.uniform(1e-4, 0.03), 8)
    dur = round(rng.uniform(1.0, 16.0), 3)
    snr = round(rng.uniform(3.0, 22.0), 3)
    fap = round(10 ** rng.uniform(-8, -1), 10)
    eo = round(rng.uniform(0.0, 6.0), 3)
    shape = round(rng.uniform(0.2, 3.0), 3)
    sec_snr = round(rng.uniform(0.0, 8.0), 3)
    sec_ratio = round(rng.uniform(0.0, 0.5), 4)
    alias = round(rng.uniform(0.0, 1.0), 4)
    coherent = 1 if rng.random() < 0.3 else 0
    centroid = round(rng.uniform(0.0, 4.0), 3)
    u1 = round(rng.uniform(0.1, 0.7), 3)
    u2 = round(rng.uniform(0.0, 0.3), 3)
    rp = round(rng.uniform(0.5, 20.0), 3)
    lines = {
        "period_days": P, "transit_depth": depth, "transit_duration_hrs": dur,
        "t0_days": 1.0, "stellar_mass_solar": M,
        "stellar_radius_solar": round(rng.uniform(0.5, 1.8), 4),
        "stellar_teff_k": round(rng.uniform(3500, 7000), 0),
        "stellar_logg": round(rng.uniform(4.0, 4.6), 3),
        "planet_radius_earth": rp, "bls_snr": snr, "bls_fap": fap,
        "even_odd_delta_sigma": eo, "shape_ratio": shape,
        "secondary_snr": sec_snr, "secondary_depth_ratio": sec_ratio,
        "alias_secondary_ratio": alias, "coherent_evidence": coherent,
        "centroid_sigma": centroid, "limb_dark_u1": u1, "limb_dark_u2": u2,
        "s_periodicity": 0.85, "s_depth": 0.85, "s_limb": 0.85, "s_stellar": 0.85,
    }
    extra = {}
    if lc:
        n_base = 120
        lines["lc"] = {
            "t": [], "f": [],
        }
        # baseline + 2 transits at t0=1.0 and 1.0+P
        rng_flux = rng.uniform
        for i in range(n_base):
            t = rng.uniform(0.0, 2.0 * P)
            lines["lc"]["t"].append(round(t, 4))
            lines["lc"]["f"].append(round(1.0, 6))
        for k in (0, 1):
            e = 1.0 + k * P
            for off in (-0.04, -0.02, 0.0, 0.02, 0.04, -0.03, -0.01, 0.01, 0.03):
                lines["lc"]["t"].append(round(e + off, 4))
                lines["lc"]["f"].append(round(1.0 - depth, 6))
    return lines, extra


def run_c(lines, lc_path=None):
    tmp = os.path.join(HERE, "tmp")
    os.makedirs(tmp, exist_ok=True)
    cand = os.path.join(tmp, "cand.txt")
    with open(cand, "w") as fh:
        for k, v in lines.items():
            if k == "lc":
                continue
            fh.write(f"{k}={v}\n")
    # Run binary directly (we're in WSL/Linux); convert paths to /mnt/ form.
    def to_wsl(p):
        if not p: return ""
        p = p.replace("\\", "/")
        if p.startswith("D:/"): return p.replace("D:/", "/mnt/d/")
        if p.startswith("C:/"): return p.replace("C:/", "/mnt/c/")
        return p
    bin_path = os.path.join(ROOT, "bin", "zspace_card")
    p = subprocess.run(
        [bin_path, cand, to_wsl(lc_path) if lc_path else ""],
        capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        raise RuntimeError(p.stderr[:2000])
    out = json.loads(p.stdout)
    # flatten
    flat = {
        "a_au": out["section_1_kepler"]["a_au"],
        "residual_si_pct": out["section_1_kepler"]["residual_si_pct"],
        "residual_solar_pct": out["section_1_kepler"]["residual_solar_pct"],
        "k": out["section_2_geometry"]["k"],
        "delta_ld_corrected": out["section_2_geometry"]["delta_ld_corrected"],
        "consistency_residual_pct": out["section_2_geometry"]["consistency_residual_pct"],
        "rp_earth": out["section_2_geometry"]["rp_earth"],
        "rho_transit_gcc": out["section_3_density"]["rho_transit_gcc"],
        "rho_tic_gcc": out["section_3_density"]["rho_tic_gcc"],
        "density_ratio": out["section_3_density"]["density_ratio"],
        "logg_calc": out["section_3_density"]["logg_calc"],
        "logg_residual": out["section_3_density"]["logg_residual"],
        "P_tr": out["section_4_probability"]["P_tr"],
        "impact_parameter_b": out["section_4_probability"]["impact_parameter_b"],
        "ingress_hrs": out["section_4_probability"]["ingress_hrs"],
        "i_min_deg": out["section_4_probability"]["i_min_deg"],
        "n_tests": out["section_5_fp_ruling"]["n_tests"],
        "n_pass": out["section_5_fp_ruling"]["n_pass"],
        "n_fail": out["section_5_fp_ruling"]["n_fail"],
        "n_critical": out["section_5_fp_ruling"]["n_critical"],
        "n_critical_pass": out["section_5_fp_ruling"]["n_critical_pass"],
        "overall_verdict": out["section_5_fp_ruling"]["overall_verdict"],
        "sovereign_verdict": out["sovereign_verdict"],
        "cvs_verdict": out["cvs_verdict"],
        "n_transits": out["n_transits"],
        "fp_verdicts": out["section_5_fp_ruling"]["fp_verdicts"],
    }
    return flat


def python_ref(lines, extra, lc_path=None):
    eng = ProofEngine(
        period_days=lines["period_days"], transit_depth=lines["transit_depth"],
        transit_duration_hrs=lines["transit_duration_hrs"],
        stellar_mass_solar=lines["stellar_mass_solar"],
        stellar_radius_solar=lines["stellar_radius_solar"],
        stellar_teff_k=lines["stellar_teff_k"],
        stellar_logg=lines["stellar_logg"],
        planet_radius_earth=lines["planet_radius_earth"],
        bls_snr=lines["bls_snr"], bls_fap=lines["bls_fap"],
        even_odd_delta_sigma=lines["even_odd_delta_sigma"],
        shape_ratio=lines["shape_ratio"],
        limb_dark_u1=lines["limb_dark_u1"], limb_dark_u2=lines["limb_dark_u2"],
    )
    k1 = eng.kepler_third_law()
    a_m = k1["a_m"]
    g2 = eng.geometric_consistency()
    d3 = eng.density_constraint(a_m)
    p4 = eng.transit_probability(a_m)

    density_ratio = d3["density_ratio"]

    # FP-10: compute the real count with the reference implementation
    from zspace_engine.validator import count_observed_transits as _cot
    n_transits_ref = None
    if lines.get("lc") and lc_path:
        import csv as _csv
        t_lc, f_lc = [], []
        with open(lc_path) as fh:
            rd = _csv.reader(fh)
            next(rd, None)
            for row in rd:
                try:
                    t_lc.append(float(row[0])); f_lc.append(float(row[1]))
                except (ValueError, IndexError):
                    continue
        n_transits_ref = _cot(t_lc, f_lc, lines["period_days"],
                              lines["t0_days"], lines["transit_duration_hrs"] / 24.0)
    elif lines.get("lc"):
        n_transits_ref = 2  # degenerate fallback (no file written)

    fp = eng.false_positive_ruling(
        secondary_snr=lines["secondary_snr"],
        centroid_sigma=lines["centroid_sigma"],
        density_ratio=density_ratio,
        is_grazing=p4["impact_parameter_b"] > 0.9,
        tic_id="TIC 999999999",
        n_transits=n_transits_ref,
        secondary_depth_ratio=lines["secondary_depth_ratio"],
        coherent_evidence=lines["coherent_evidence"],
        alias_secondary_ratio=lines["alias_secondary_ratio"],
    )

    cvs = (0.97 * lines["s_periodicity"] + 0.83 * lines["s_depth"] +
           0.61 * lines["s_limb"] + 0.31 * lines["s_stellar"]) / \
          (0.97 + 0.83 + 0.61 + 0.31)

    fp_verdicts = [1 if t["verdict"] == "PASS" else 0 for t in fp["tests"]]
    if len(fp_verdicts) < 12:
        fp_verdicts.append(-1)  # FP-10 skipped

    ref = {
        "a_au": k1["a_au"],
        "k": g2["k"],
        "delta_ld_corrected": g2["delta_ld_corrected"],
        "rp_earth": g2["R_p_earth"],
        "rho_transit_gcc": d3["rho_transit_gcc"],
        "rho_tic_gcc": d3["rho_tic_gcc"],
        "density_ratio": density_ratio,
        "P_tr": p4["P_tr"],
        "impact_parameter_b": p4["impact_parameter_b"],
        "ingress_hrs": p4["T_ingress_hrs"],
        "i_min_deg": p4["i_min_deg"],
        "n_tests": fp["n_tests"],
        "n_pass": fp["n_pass"],
        "n_fail": fp["n_fail"],
        "n_critical": fp["n_critical"],
        "n_critical_pass": fp["n_critical_pass"],
        "overall_verdict": fp["overall_verdict"],
        "sovereign_verdict": fp["overall_verdict"],
        "cvs_verdict": "PLANET" if cvs > 0.80 else "NOT_PLANET",
        "n_transits": (n_transits_ref if n_transits_ref is not None else -1),
        "fp_verdicts": fp_verdicts,
    }
    return ref


def compare(flat_c, ref, tol_rel=2e-3, tol_abs=0.02, verbose=False):
    problems = []
    # per-test verdict comparison (FP-1..FP-11: 11 fixed + FP-10 slot)
    names = ["FP1", "FP2", "FP3", "FP4", "FP5", "FP5b", "FP5c",
             "FP6", "FP7", "FP8", "FP9", "FP10"]
    vc = flat_c.pop("fp_verdicts", None)
    vr = ref.pop("fp_verdicts", None)
    if vc is not None and vr is not None:
        for i, (a, b) in enumerate(zip(vc, vr)):
            if a != b:
                problems.append((f"fp_verdict[{names[i]}]", a, b, "verdict"))
    for f in FIELDS_NUM:
        vc, vr = flat_c[f], ref[f]
        if not math.isfinite(vc) or not math.isfinite(vr):
            if vc != vr:
                problems.append((f, vc, vr, "non-finite mismatch"))
            continue
        if vr == 0.0:
            if abs(vc) > tol_abs:
                problems.append((f, vc, vr, "abs"))
        else:
            rel = abs(vc - vr) / max(abs(vr), 1e-12)
            if rel > tol_rel and abs(vc - vr) > tol_abs:
                problems.append((f, vc, vr, f"rel={rel:.2e}"))
    for f in FIELDS_INT:
        if flat_c[f] != ref[f]:
            problems.append((f, flat_c[f], ref[f], "int"))
    for f in FIELDS_STR:
        if flat_c[f] != ref[f]:
            problems.append((f, flat_c[f], ref[f], "str"))
    return problems


def write_lc(lines, path):
    with open(path, "w") as fh:
        fh.write("time,flux,flux_err,model\n")
        for t, f in zip(lines["lc"]["t"], lines["lc"]["f"]):
            fh.write(f"{t},{f},0.001,1.0\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--seed", type=int, default=20260817)
    ap.add_argument("--lc", action="store_true", help="include light-curve (FP-10) cases")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    n_lc_cases = args.n // 2 if args.lc else 0

    n_pass = 0
    n_bad = 0
    for i in range(args.n):
        lines, extra = make_candidate(rng, lc=(args.lc and i < n_lc_cases))
        lc_path = None
        if lines.get("lc"):
            lc_path = os.path.join(HERE, "tmp", "lc.csv")
            os.makedirs(os.path.dirname(lc_path), exist_ok=True)
            write_lc(lines, lc_path)
        flat_c = run_c(lines, lc_path)
        ref = python_ref(lines, extra, lc_path)
        probs = compare(flat_c, ref, verbose=args.verbose)
        if probs:
            n_bad += 1
            print(f"[{i}] MISMATCH:")
            print(f"    C_full: {flat_c}")
            print(f"    py_full: {ref}")
            print(f"    cand={ {k: v for k, v in lines.items() if k != 'lc'} }")
            for f, vc, vr, kind in probs:
                print(f"    {f}: C={vc!r} py={vr!r} ({kind})")
        else:
            n_pass += 1
    print(f"parity: {n_pass}/{args.n} matched, {n_bad} mismatched")


if __name__ == "__main__":
    main()