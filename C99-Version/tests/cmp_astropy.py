"""Differential verification: C99 zspace_card bls vs astropy BoxLeastSquares.

Runs both on the same full-range periodogram (same freq grid, same duration
ladder) and compares the best-fit parameters point by point.
"""
import json
import os
import subprocess
import sys

import numpy as np
from astropy.timeseries import BoxLeastSquares
import astropy.units as u

HERE = os.path.dirname(os.path.abspath(__file__))
LC = os.path.join(HERE, "tmp_lc.csv")
BIN = os.path.join(HERE, "..", "build", "zspace_card")

arr = np.loadtxt(LC, skiprows=1, delimiter=",")
t_obs, flux = arr[:, 0], arr[:, 1]
print(f"LC: {len(t_obs)} points")

P_MIN, P_MAX = 0.5, 13.5

baseline = float(t_obs[-1] - t_obs[0])
freq_min = 1.0 / P_MAX
freq_max = 1.0 / P_MIN
df = 1.0 / (baseline * 10.0)
n_freqs = max(int((freq_max - freq_min) / df), 2000)
freq_grid = np.linspace(freq_min, freq_max, n_freqs)
period_grid = (1.0 / freq_grid[::-1]) * u.day

durations = np.array([0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,
                      2.5, 3, 4, 5, 6, 8, 12]) / 24.0
durations = durations[durations * 24.0 < max(P_MIN * 24.0, 0.5)]
durations = durations * u.day

model = BoxLeastSquares(t=t_obs * u.day, y=flux * u.dimensionless_unscaled)
pg = model.power(period=period_grid, duration=durations)
power = np.asarray(pg.power, dtype=np.float64).ravel()
period_arr = np.asarray(pg.period.to("d").value, dtype=np.float64).ravel()
depth_arr = np.asarray(pg.depth, dtype=np.float64).ravel()
depth_snr = np.asarray(pg.depth_snr, dtype=np.float64).ravel()
transit_time = np.asarray(pg.transit_time.to("d").value, dtype=np.float64).ravel()
duration_best = np.asarray(pg.duration.to("d").value, dtype=np.float64).ravel()
depth_err = np.asarray(pg.depth_err, dtype=np.float64).ravel()

idx = int(np.argmax(power))
print(f"astropy best: P={period_arr[idx]:.6f}  power={power[idx]:.6f}  "
      f"depth={depth_arr[idx]:.6f}  snr={depth_snr[idx]:.4f}  "
      f"dur={duration_best[idx]*24:.2f}h  t0={transit_time[idx]:.4f}")

res = subprocess.run([BIN, "bls", LC, str(P_MIN), str(P_MAX)],
                     capture_output=True, text=True)
c = json.loads(res.stdout)
print(f"C99     best: P={c['period_days']:.6f}  power={c['power']:.6f}  "
      f"depth={c['depth']:.6f}  snr={c['snr']:.4f}  "
      f"dur={c['duration_hrs']:.2f}h  t0={c['t0_days']:.4f}")

tols = dict(period=0.02, power=1e-3, depth=1e-3, dur=0.6, t0=0.15)
ok = True
for key, a, b in [("period", period_arr[idx], c["period_days"]),
                  ("power", power[idx], c["power"]),
                  ("depth", depth_arr[idx], c["depth"]),
                  ("dur", duration_best[idx] * 24, c["duration_hrs"]),
                  ("t0", transit_time[idx], c["t0_days"])]:
    d = abs(a - b)
    status = "OK " if d <= tols[key] else "FAIL"
    ok &= d <= tols[key]
    print(f"  {key:6s} astropy={a:.6f}  C99={b:.6f}  diff={d:.6f}  [{status}]")
# SNR is NOT comparable: astropy depth_snr uses ivar=1 (unitless), while C99
# reports the ZS matched-filter SNR (detectors.py): depth/(sigma_oot*sqrt(1/n)).
print(f"  snr    C99(ZS matched-filter)={c['snr']:.1f}  astropy(depth_snr)={depth_snr[idx]:.4f}  [different definitions]")

# periodogram shape agreement at top-10 frequencies
top = np.argsort(power)[-10:][::-1]
print("top-10 period agreement:", "OK" if any(
    abs(period_arr[i] - c["period_days"]) < 0.02 for i in top) else "FAIL")
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)