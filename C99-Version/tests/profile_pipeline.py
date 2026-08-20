"""Profile where Python pipeline time goes (synthetic run)."""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
os.environ.setdefault("ZSPACE_SKIP_LOGS", "1")

from run_pipeline import run_synthetic_test

t0 = time.perf_counter()
res = run_synthetic_test(engine="python")
dt = time.perf_counter() - t0
print(f"\nTOTAL pipeline: {dt*1000:.0f} ms")
print(f"verdict: {res.get('verdict')}")
print(f"cvs: {res.get('cvs')}")
for k, v in res.items():
    if isinstance(v, (int, float, str)) and k not in ("verdict", "cvs"):
        print(f"  {k}: {v}")