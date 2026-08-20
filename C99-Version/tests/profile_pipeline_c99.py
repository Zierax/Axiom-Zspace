import re
import subprocess

r = subprocess.run(["python3", "run_pipeline.py", "--synthetic", "--engine", "c99"],
                   capture_output=True, text=True, timeout=900)
log = r.stdout + r.stderr
times = {}
for line in log.splitlines():
    m = re.search(r"(\d{2}:\d{2}:\d{2},\d{3}).*?(INFO|WARNING|ERROR) - (.*)", line)
    if m:
        h, mi, s = m.group(1).split(",")[0].split(":")
        t = int(h) * 3600 + int(mi) * 60 + float(s) + float("0." + m.group(1).split(",")[1])
        msg = m.group(3).strip()
        times.setdefault(t, []).append(msg)
keys = sorted(times)
print("PHASE timing (c99 engine):")
prev = None
for t in keys:
    for msg in times[t]:
        if msg.startswith("[PHASE") or "Gate:" in msg or "Verdict" in msg or "Period:" in msg:
            dt = f"+{(t - prev):.3f}s" if prev is not None else "(start)"
            print(f"  {t:9.3f}  {dt:>8}  {msg[:80]}")
            prev = t
print(f"total: {keys[-1] - keys[0]:.2f}s")
