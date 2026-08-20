"""Speed bench: C99 zspace_card eph vs Python EphemerisResolver (same LC)."""
import csv
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
C99_DIR = os.path.abspath(os.path.join(HERE, ".."))
BIN = os.path.join(C99_DIR, "build", "zspace_card")
LC = os.path.join(HERE, "tmp_lc.csv")

sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
from zspace_engine.ephemeris import EphemerisResolver

rows = []
with open(LC) as fh:
    reader = csv.reader(fh)
    next(reader)
    for r in reader:
        rows.append((float(r[0]), float(r[1])))
t_np = [r[0] for r in rows]
f_np = [r[1] for r in rows]
import numpy as np
t_np = np.asarray(t_np)
f_np = np.asarray(f_np)

PMIN, PMAX = 0.5, 25.0
CASES = [
    (3.70846503, 3.63724088),
    (7.41693006, 3.63724088),
    (1.85423252, 3.63724088),
]

resolver = EphemerisResolver()

# warmup
resolver.resolve(t_np, f_np, 3.70846503, 3.63724088, period_min=PMIN, period_max=PMAX)

py_times = []
for period, t0 in CASES:
    t0w = time.perf_counter()
    resolver.resolve(t_np, f_np, period, t0, period_min=PMIN, period_max=PMAX)
    py_times.append(time.perf_counter() - t0w)

c_times = []
for period, t0 in CASES:
    subprocess.run([BIN, "eph", LC, f"{period:.8f}", f"{t0:.8f}", f"{PMIN}", f"{PMAX}"],
                   capture_output=True, text=True)
    t0w = time.perf_counter()
    out = subprocess.run([BIN, "eph", LC, f"{period:.8f}", f"{t0:.8f}", f"{PMIN}", f"{PMAX}"],
                         capture_output=True, text=True)
    c_times.append(time.perf_counter() - t0w)
    json.loads(out.stdout)

print(f"LC: {len(rows)} points, {len(CASES)} resolve calls each")
print(f"Python EphemerisResolver: best {min(py_times)*1000:.1f} ms  "
      f"(mean {sum(py_times)/len(py_times)*1000:.1f} ms)")
print(f"C99   zspace_card eph   : best {min(c_times)*1000:.1f} ms  "
      f"(mean {sum(c_times)/len(c_times)*1000:.1f} ms)")
print(f"Speedup (best/best): {min(py_times)/min(c_times):.1f}x")