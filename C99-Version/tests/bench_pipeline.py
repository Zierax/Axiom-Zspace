import statistics
import subprocess
import time

res = {}
for engine in ("python", "c99"):
    ts = []
    for _ in range(3):
        t0 = time.perf_counter()
        subprocess.run(["python3", "run_pipeline.py", "--synthetic", "--engine", engine],
                       capture_output=True, text=True, timeout=900)
        ts.append(time.perf_counter() - t0)
    res[engine] = ts
    print(f"{engine}: runs={['%.2f' % t for t in ts]} median={statistics.median(ts):.2f}s")