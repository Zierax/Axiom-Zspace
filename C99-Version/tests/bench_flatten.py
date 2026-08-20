import csv
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
import numpy as np

from c99_bridge import run_c99_flatten

rows = []
with open("C99-Version/tests/tmp_lc.csv") as fh:
    rd = csv.reader(fh)
    next(rd)
    for r in rd:
        rows.append((float(r[0]), float(r[1])))
t = np.asarray([r[0] for r in rows])
f = np.asarray([r[1] for r in rows])

t0 = time.perf_counter()
out = run_c99_flatten(t, f, 3.7)
dt = time.perf_counter() - t0
print(f"run_c99_flatten: {dt*1000:.0f}ms, n_out={len(out)}")

from zspace_engine.ingestion import LightCurveIngester
ing = LightCurveIngester(tic_id="x")
t0 = time.perf_counter()
ff, tr, w = ing._savgol_flatten(t, f, 3.7)
dt = time.perf_counter() - t0
print(f"scipy _savgol_flatten: {dt*1000:.0f}ms, w={w}")

d = np.abs(np.asarray(out) - ff)
print(f"max diff C vs scipy: {d.max():.3e}")
