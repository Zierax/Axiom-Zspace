"""Time each TransitAuditor audit on the synthetic pipeline inputs."""
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np
from zspace_engine.auditors import TransitAuditor
from zspace_engine.detectors import BLSDetector
from run_pipeline import generate_synthetic_transit

t, f = generate_synthetic_transit()
period, t0, dur = 3.695262, 2459300.0, 4.0 / 24.0

aud = TransitAuditor(run_mcmc=False)

t0w = time.perf_counter()
r1 = aud.secondary_eclipse_test(t, f, period, t0, dur)
dt1 = time.perf_counter() - t0w

t0w = time.perf_counter()
r2 = aud.even_odd_test(t, f, period, t0, dur)
dt2 = time.perf_counter() - t0w

t0w = time.perf_counter()
r3 = aud.depth_consistency_score(t, f, period, t0, dur)
dt3 = time.perf_counter() - t0w

t0w = time.perf_counter()
r4 = aud.limb_shape_score(period, dur, 1e-3, time=t, flux=f, t0=t0)
dt4 = time.perf_counter() - t0w

print(f"secondary_eclipse : {1000*dt1:7.1f} ms  sec_ratio={r1.secondary_ratio:.3f} snr={r1.secondary_snr:.1f}")
print(f"even_odd          : {1000*dt2:7.1f} ms  d_sigma={r2.delta_sigma:.2f} flag={r2.is_eb_flag}")
print(f"depth_consistency : {1000*dt3:7.1f} ms")
print(f"limb_shape        : {1000*dt4:7.1f} ms")
print(f"TOTAL             : {1000*(dt1+dt2+dt3+dt4):7.1f} ms")