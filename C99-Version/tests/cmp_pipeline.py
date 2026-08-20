import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
import numpy as np
from scipy.signal import savgol_filter

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "build", "zspace_card")
LC = os.path.join(HERE, "tmp_lc.csv")

rows = []
with open(LC) as fh:
    rd = csv.reader(fh)
    next(rd)
    for r in rd:
        rows.append((float(r[0]), float(r[1])))
t = np.asarray([r[0] for r in rows])
f = np.asarray([r[1] for r in rows])

out = subprocess.run([BIN, "pipeline", LC, "0.5", "13.5", "3.7"],
                     capture_output=True, text=True)
d = json.loads(out.stdout)
period = d["period_days"]

med = np.median(f)
norm = f / med

order = np.argsort(t, kind="mergesort")
ts = t[order]
fs = norm[order]
cadence = float(np.median(np.diff(ts)))
w = int(round(0.75 * period / cadence))
if w < 51:
    w = 51
if w % 2 == 0:
    w += 1
n = fs.size
w = min(w, n if n % 2 == 1 else n - 1)
po = min(3, w - 1)
trend = savgol_filter(fs, window_length=w, polyorder=po)
ffs = fs / trend
inv = np.empty_like(order)
inv[order] = np.arange(order.size)
py_flat = ffs[inv]

c_flat = np.asarray(d["flattened"])
diff = np.abs(c_flat - py_flat)
print(f"window C={d['window_pts']} PY={w} poly C={d['polyorder']} PY={po}")
print(f"flattened max diff: {diff.max():.3e}  {'OK' if diff.max() < 1e-8 else '**DIFF**'}")

print("bls:", {k: d[k] for k in ("period_days", "snr", "fap", "duration_hrs", "t0_days", "depth")})
print("even_odd delta_sigma:", d["even_odd"]["delta_sigma"])
print("depth cv:", d["depth_consistency"]["cv"], "s_depth:", d["depth_consistency"]["s_depth"])
print("sec snr:", d["secondary_eclipse"]["secondary_snr"])
print("ie:", d["ingress_egress"]["ingress_fraction"], d["ingress_egress"]["flat_fraction"],
      d["ingress_egress"]["is_v_shape"], d["ingress_egress"]["fp_risk"])