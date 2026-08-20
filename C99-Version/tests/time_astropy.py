import os, time as _time
import numpy as np
from astropy.timeseries import BoxLeastSquares
import astropy.units as u

HERE = os.path.dirname(os.path.abspath(__file__))
arr = np.loadtxt(os.path.join(HERE, "tmp_lc.csv"), skiprows=1, delimiter=",")
t_obs, flux = arr[:, 0], arr[:, 1]
print(f"LC: {len(t_obs)} points")

durations = np.array([0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,
                      2.5, 3, 4, 5, 6, 8, 12]) / 24.0
durations = durations[durations * 24.0 < max(0.5 * 24.0, 0.5)]
durations = durations * u.day

model = BoxLeastSquares(t=t_obs * u.day, y=flux * u.dimensionless_unscaled)
period = np.linspace(0.5, 13.5, 2000) * u.day

t0 = _time.time()
pg = model.power(period=period, duration=durations)
el = _time.time() - t0
i = int(np.argmax(pg.power))
print(f"astropy elapsed: {el:.2f}s  best P={pg.period[i].value:.4f}  power={pg.power[i]:.4f}")