"""Time savgol re-flatten vs remaining PHASE-2 pieces on pipeline inputs."""
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np
from run_pipeline import generate_synthetic_transit
from zspace_engine.ingestion import LightCurveIngester

t, f = generate_synthetic_transit()
period = 3.695262

ing = LightCurveIngester(tic_id="SYNTHETIC")
t0w = time.perf_counter()
f2, _, _ = ing._savgol_flatten(t, f, period)
dt = time.perf_counter() - t0w
print(f"savgol re-flatten : {1000*dt:7.1f} ms  ({len(t)} pts)")

t0w = time.perf_counter()
for _ in range(5):
    f2b, _, _ = ing._savgol_flatten(t, f, period)
print(f"savgol x5         : {1000*(time.perf_counter()-t0w):7.1f} ms")