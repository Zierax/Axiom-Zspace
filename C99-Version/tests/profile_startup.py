import subprocess
import time

t0 = time.perf_counter()
r = subprocess.run(["python3", "run_pipeline.py", "--synthetic", "--engine", "c99"],
                   capture_output=True, text=True, timeout=900)
t1 = time.perf_counter()
print(f"total wall: {t1 - t0:.2f}s")
log = r.stdout + r.stderr
import re
first = None
last = None
for line in log.splitlines():
    m = re.search(r"(\d{2}:\d{2}:\d{2},\d{3})", line)
    if m:
        h, mi, s = m.group(1).split(",")[0].split(":")
        t = int(h) * 3600 + int(mi) * 60 + float(s) + float("0." + m.group(1).split(",")[1])
        if first is None:
            first = t
        last = t
print(f"logged span: {last - first:.2f}s  -> pre-log startup: {(t1 - t0) - (last - first):.2f}s")
