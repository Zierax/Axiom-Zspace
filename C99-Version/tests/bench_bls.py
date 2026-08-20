"""bench_bls.py  ·  BLS parity + speed test (C99 vs astropy)
Generates a synthetic light curve with a known transit period,
runs the C99 BLS engine and astropy BLS, compares recovered period.
"""
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "build", "zspace_card")

def make_lc(period=3.69526, depth=0.0042, dur_hrs=2.0, n=20000,
            baseline=30.0, rng=None):
    rng = rng or np.random.default_rng(20260817)
    t = np.sort(rng.uniform(0.0, baseline, n))
    phase = (t / period) % 1.0
    in_t = (phase < dur_hrs / 24.0 / period) | (phase > 1.0 - dur_hrs / 24.0 / period)
    flux = np.ones(n) - depth * in_t
    flux += rng.normal(0.0, 0.0004, n)
    return t, flux

def run_c(t, f, pmin, pmax):
    lc = os.path.join(HERE, "tmp_lc.csv")
    np.savetxt(lc, np.column_stack([t, f]), fmt="%.8f", delimiter=",",
               header="time,flux", comments="")
    p = subprocess.run([BIN, "bls", lc, str(pmin), str(pmax)],
                       capture_output=True, text=True, timeout=600)
    if p.returncode != 0:
        return None, p.stderr
    return json.loads(p.stdout), None

def run_py(t, f, pmin, pmax):
    from astropy.timeseries import BoxLeastSquares
    import astropy.units as u
    baseline = float(t[-1] - t[0])
    freq_min = 1.0 / pmax
    freq_max = 1.0 / pmin
    df = 1.0 / (baseline * 10.0)
    n_freqs = max(int((freq_max - freq_min) / df), 2000)
    freq_grid = np.linspace(freq_min, freq_max, n_freqs)
    period_grid = (1.0 / freq_grid[::-1]) * u.day
    durations = np.array([0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,
                          2.5, 3, 4, 5, 6, 8, 12]) / 24.0
    durations = durations[durations * 24.0 < max(pmin * 24.0, 0.5)]
    durations = durations * u.day
    model = BoxLeastSquares(t=t * u.day, y=f * u.dimensionless_unscaled)
    pg = model.power(period=period_grid, duration=durations)
    power = np.asarray(pg.power, dtype=np.float64).ravel()
    period_arr = np.asarray(pg.period.to("d").value, dtype=np.float64).ravel()
    idx = int(np.argmax(power))
    return {
        "period_days": float(period_arr[idx]),
        "power": float(power[idx]),
    }

def main():
    truth = 3.69526
    t, f = make_lc(period=truth)
    print(f"LC: n={len(t)} points, truth P={truth} d")

    best_tc = None
    cres = None
    for _ in range(3):
        t0 = time.perf_counter()
        r, cerr = run_c(t, f, 0.5, 13.5)
        el = time.perf_counter() - t0
        if cerr:
            print("C ERROR:", cerr)
            sys.exit(1)
        if best_tc is None or el < best_tc:
            best_tc, cres = el, r
    tc = best_tc

    t0 = time.perf_counter()
    pres = run_py(t, f, 0.5, 13.5)
    tp = time.perf_counter() - t0

    print(f"C99 : P={cres['period_days']:.6f} d  power={cres['power']:.4f}  "
          f"dur={cres['duration_hrs']:.2f} h  snr={cres['snr']:.1f}  [{tc*1000:.0f} ms]")
    print(f"Py  : P={pres['period_days']:.6f} d  power={pres['power']:.4f}  [{tp:.2f} s]")

    p_err = abs(cres['period_days'] - pres['period_days']) / pres['period_days'] * 100.0
    print(f"period match: {p_err:.4f}% off  |  speedup: {tp/tc:.0f}x")

if __name__ == "__main__":
    main()