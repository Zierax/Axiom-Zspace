"""Differential test: C99 zspace_card audit vs Python TransitAuditor."""
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
C99_DIR = os.path.abspath(os.path.join(HERE, ".."))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
BIN = os.path.join(C99_DIR, "build", "zspace_card")
LC = os.path.join(HERE, "tmp_lc.csv")

sys.path.insert(0, ROOT)
import numpy as np
from zspace_engine.auditors import TransitAuditor

rows = []
with open(LC) as fh:
    reader = csv.reader(fh)
    next(reader)
    for r in reader:
        rows.append((float(r[0]), float(r[1])))
t = np.asarray([r[0] for r in rows])
f = np.asarray([r[1] for r in rows])
print(f"LC: {len(rows)} points")

CASES = [
    (3.70846503, 3.63724088, 4.0, 0.0056, "clean planet"),
    (3.695262, 2459300.0 + 1.2, 4.0, 0.005, "pipeline-style"),
    (7.41693006, 3.63724088, 4.0, 0.0056, "alias period"),
]

aud = TransitAuditor(run_mcmc=False)
all_pass = True
for period, t0, dur_hrs, depth, label in CASES:
    duration = dur_hrs / 24.0

    # Python reference
    eo = aud.even_odd_test(t, f, period, t0, duration)
    dc = aud.depth_consistency_score(t, f, period, t0, duration)
    se = aud.secondary_eclipse_test(t, f, period, t0, duration)
    from zspace_engine.detectors import BLSDetector
    bp, bf, _ = BLSDetector.fold_and_bin(t, f, period, t0, n_bins=200)
    ie = aud.ingress_egress_test(bp, bf, period, duration, depth)

    out = subprocess.run(
        [BIN, "audit", LC, f"{period:.8f}", f"{t0:.8f}", f"{dur_hrs:.6f}", f"{depth:.8f}"],
        capture_output=True, text=True)
    c = json.loads(out.stdout)
    eo_c, dc_c, se_c, ie_c = c["even_odd"], c["depth_consistency"], c["secondary_eclipse"], c["ingress_egress"]

    import re as _re
    m = _re.search(r"χ²_red=([\d.]+)", dc.proof)
    py_chi2 = float(m.group(1)) if m else float("nan")

    tol = 1e-4
    ok = True
    checks = [
        ("even_odd delta_sigma", eo_c["delta_sigma"], eo.delta_sigma, 1e-3),
        ("even_odd n_even", eo_c["n_even"], eo.n_even, 0),
        ("even_odd n_odd", eo_c["n_odd"], eo.n_odd, 0),
        ("even_odd depth_even", eo_c["depth_even"], eo.depth_even, tol),
        ("even_odd depth_odd", eo_c["depth_odd"], eo.depth_odd, tol),
        ("even_odd is_eb", eo_c["is_eb_flag"], eo.is_eb_flag, None),
        ("depth n", dc_c["n_transits"], dc.depths.size, 0),
        ("depth mean", dc_c["mean_depth"], dc.mean_depth, 1e-6),
        ("depth cv", dc_c["cv"], dc.cv, 1e-3),
        ("depth chi2_red", dc_c["chi2_red"], py_chi2, 2e-3),
        ("depth s_depth", dc_c["s_depth"], dc.s_depth, 1e-6),
        ("sec primary", se_c["primary_depth"], se.primary_depth, tol),
        ("sec secondary", se_c["secondary_depth"], se.secondary_depth, 1e-7),
        ("sec ratio", se_c["secondary_ratio"], se.secondary_ratio, 1e-3),
        ("sec snr", se_c["secondary_snr"], se.secondary_snr, 0.1),
        ("ie ingress_frac", ie_c["ingress_fraction"], ie.ingress_fraction, 0.05),
        ("ie flat_frac", ie_c["flat_fraction"], ie.flat_fraction, 0.05),
        ("ie is_v_shape", ie_c["is_v_shape"], ie.is_v_shape, None),
        ("ie fp_risk", ie_c["fp_risk"], ie.fp_risk, None),
    ]
    print(f"\n[{label}] period={period:.6f}")
    for name, a, b, tol_v in checks:
        if tol_v is None:
            match = (a == b)
        elif tol_v == 0:
            match = (int(a) == int(b))
        else:
            match = abs(float(a) - float(b)) < tol_v
        ok &= match
        print(f"  {name:<24} C={a!s:<14} PY={b!s:<14} {'OK' if match else '**DIFF**'}")
    all_pass &= ok

print(f"\n{'ALL PASS' if all_pass else 'FAILURES PRESENT'}")
sys.exit(0 if all_pass else 1)