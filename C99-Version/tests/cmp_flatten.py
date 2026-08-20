import csv
import json
import os
import subprocess

import numpy as np
from scipy.signal import savgol_filter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
BIN = os.path.join(ROOT, "C99-Version", "build", "zspace_card")
LC = os.path.join(HERE, "tmp_lc.csv")

rows = []
with open(LC) as fh:
    reader = csv.reader(fh)
    next(reader)
    for r in reader:
        rows.append((float(r[0]), float(r[1])))
t = np.asarray([r[0] for r in rows])
f = np.asarray([r[1] for r in rows])

def py_flatten(time, flux, period_days):
    order = np.argsort(time, kind="mergesort")
    time_sorted = time[order]
    flux_sorted = flux[order]
    cadence = float(np.median(np.diff(time_sorted)))
    if period_days is not None and period_days > 0:
        window_days = 0.75 * period_days
        window_pts = int(round(window_days / cadence))
    else:
        window_pts = int(round(3.0 / cadence))
    if window_pts < 51:
        window_pts = 51
    if window_pts % 2 == 0:
        window_pts += 1
    n_pts = flux_sorted.size
    window_pts = min(window_pts, n_pts if n_pts % 2 == 1 else n_pts - 1)
    if window_pts < 5 or n_pts < 5:
        return flux, 0
    polyorder = min(3, window_pts - 1)
    trend = savgol_filter(flux_sorted, window_length=window_pts, polyorder=polyorder)
    flux_flat_sorted = flux_sorted / trend
    inv_order = np.empty_like(order)
    inv_order[order] = np.arange(order.size)
    return flux_flat_sorted[inv_order], window_pts

all_pass = True
for P in (3.7, 0.5, 13.5, 0.0, 1.5):
    out_csv = os.path.join(HERE, "flat_out.csv")
    r = subprocess.run([BIN, "flatten", LC, out_csv, str(P)],
                       capture_output=True, text=True)
    meta = json.loads(r.stdout)
    c_flat = np.asarray([float(x[1]) for x in csv.reader(open(out_csv))
                         if x and x[0] != "time"])
    py_flat, py_w = py_flatten(t, f, P)
    diff = np.abs(c_flat - py_flat)
    rel = diff / np.maximum(np.abs(py_flat), 1e-12)
    worst = float(np.max(diff))
    worst_rel = float(np.max(rel))
    ok = worst < 1e-8 and meta["window_pts"] == py_w
    all_pass &= ok
    print(f"P={P}: window C={meta['window_pts']} PY={py_w}  max_abs={worst:.3e} "
          f"max_rel={worst_rel:.3e}  {'OK' if ok else '**DIFF**'}")

print("ALL PASS" if all_pass else "FAILURES PRESENT")