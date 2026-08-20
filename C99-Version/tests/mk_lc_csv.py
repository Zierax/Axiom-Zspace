import numpy as np

n_points = 3000
baseline = 60.0
period_true = 2.0
depth_true = 0.0053
dur_hrs_true = 2.5

rng = np.random.default_rng(1000)
t_arr = baseline * np.arange(n_points) / n_points
phase = ((t_arr - 1.0) / period_true) % 1.0
dur_frac = (dur_hrs_true / 24.0) / period_true
in_transit = (phase < dur_frac) | (phase > 1.0 - dur_frac)
signal = np.where(in_transit, -depth_true, 0.0)
noise = (rng.random(n_points) - 0.5) * 2.0 * 0.0015
flux = 1.0 + signal + noise

with open('/tmp/lc_syn2.csv', 'w') as f:
    f.write('time,flux\n')
    for t, fl in zip(t_arr, flux):
        f.write(f'{t:.8f},{fl:.10f}\n')
print('written /tmp/lc_syn2.csv')