"""Per-audit timing via log-line handler inside the real pipeline."""
import logging
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("ZSPACE_SKIP_LOGS", "1")

from run_pipeline import run_synthetic_test

marks = []


class H(logging.Handler):
    def emit(self, record):
        msg = record.getMessage()
        for key in ("Even/Odd Delta-sigma", "Depth CV", "Shape ratio", "Ingress fraction",
                    "PHASE 1", "PHASE 2", "PHASE 2.5", "PHASE 3", "TIC source"):
            if key in msg:
                marks.append((time.perf_counter(), key))
                break


logging.getLogger().addHandler(H())
for name in list(logging.root.manager.loggerDict):
    logging.getLogger(name).addHandler(H())

t0 = time.perf_counter()
res = run_synthetic_test(engine="python")
dt = time.perf_counter() - t0

prev = marks[0][0]
for t, key in marks:
    print(f"{key:<32} {1000*(t-prev):8.0f} ms")
    prev = t
print(f"{'TOTAL pipeline':<32} {1000*dt:8.0f} ms")