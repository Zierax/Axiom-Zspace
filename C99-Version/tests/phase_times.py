"""Per-phase timing inside the pipeline process (logging handler)."""
import logging
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("ZSPACE_SKIP_LOGS", "1")

from run_pipeline import run_synthetic_test

marks = []


class TimestampHandler(logging.Handler):
    def emit(self, record):
        msg = record.getMessage()
        if "[PHASE" in msg or "TOTAL" in msg:
            marks.append((time.perf_counter(), msg.split("]")[0].split("[")[1]))


logging.getLogger().addHandler(TimestampHandler())
for name in list(logging.root.manager.loggerDict):
    logging.getLogger(name).addHandler(TimestampHandler())

t0 = time.perf_counter()
res = run_synthetic_test(engine="python")
dt = time.perf_counter() - t0

prev = t0
for t, label in marks:
    print(f"{label:<45} {1000*(t-prev):8.0f} ms")
    prev = t
print(f"{'TOTAL pipeline':<45} {1000*dt:8.0f} ms")
print(f"verdict: {res.get('verdict')}")