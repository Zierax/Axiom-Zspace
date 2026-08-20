"""Differential test: C99 zspace_card eph vs Python EphemerisResolver."""
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
C99_DIR = os.path.abspath(os.path.join(HERE, ".."))
BIN = os.path.join(C99_DIR, "build", "zspace_card")
LC = os.path.join(HERE, "tmp_lc.csv")

sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
from zspace_engine.ephemeris import EphemerisResolver

rows = []
with open(LC) as fh:
    reader = csv.reader(fh)
    header = next(reader)
    for r in reader:
        rows.append((float(r[0]), float(r[1])))
time = [r[0] for r in rows]
flux = [r[1] for r in rows]
import numpy as np
t_np = np.asarray(time)
f_np = np.asarray(flux)
print(f"LC: {len(rows)} points")

PMIN, PMAX = 0.5, 25.0
CASES = [
    (3.70846503, 3.63724088, "fundamental case"),
    (7.41693006, 3.63724088, "2P alias (P*2)"),
    (11.12539509, 3.63724088, "3P alias (P*3)"),
    (1.85423252, 3.63724088, "P/2 over-harmonic"),
    (1.23615501, 3.63724088, "P/3 over-harmonic"),
    (4.20000000, 3.63724088, "noise period"),
]

resolver = EphemerisResolver()
all_pass = True
for period, t0, label in CASES:
    py = resolver.resolve(t_np, f_np, period, t0, period_min=PMIN, period_max=PMAX)
    out = subprocess.run(
        [BIN, "eph", LC, f"{period:.8f}", f"{t0:.8f}", f"{PMIN}", f"{PMAX}"],
        capture_output=True, text=True)
    cj = json.loads(out.stdout)
    p_class = cj["classifier"]
    p_mult = cj["multiple"]
    p_phys = cj["physical_period"]
    p_conf = cj["confidence"]
    p_pat = cj["pattern"].replace("'", '"')
    try:
        p_pat = json.loads(p_pat)
    except Exception:
        pass
    ok = (p_class == py.classifier and p_mult == py.multiple
          and abs(p_phys - py.physical_period) < 1e-6
          and abs(p_conf - py.confidence) < 1e-3
          and str(p_pat) == str(py.pattern))
    all_pass &= ok
    print(f"\n[{label}] period={period:.8f}")
    print(f"  C99 : class={p_class:<16} mult={p_mult} phys={p_phys:.8f} "
          f"conf={p_conf:.3f} pattern={p_pat}")
    print(f"  PY  : class={py.classifier:<16} mult={py.multiple} phys={py.physical_period:.8f} "
          f"conf={py.confidence:.3f} pattern={py.pattern}")
    print(f"  MATCH: {ok}")
    if not ok:
        print(f"  C99 evidence: {cj['evidence']}")
        print(f"  PY  evidence: {py.evidence}")
        print(f"  C99 flags   : {cj['flags']}")
        print(f"  PY  flags   : {py.flags}")

print(f"\n{'ALL PASS' if all_pass else 'FAILURES PRESENT'}")
sys.exit(0 if all_pass else 1)
