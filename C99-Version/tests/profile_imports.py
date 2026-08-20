import subprocess
import time

mods = [
    "run_pipeline.py",  # full module import (before main guard?)
]
code = """
import time
t0 = time.perf_counter()
import run_pipeline
print(f"import run_pipeline: {time.perf_counter() - t0:.2f}s")
"""
t0 = time.perf_counter()
r = subprocess.run(["python3", "-c", code], capture_output=True, text=True, timeout=300, cwd="/mnt/d/Axioms/Axiom-Zspace-CODE")
print(r.stdout)
print(r.stderr[-500:] if r.returncode else "")

code2 = """
import time
t0 = time.perf_counter()
import zspace_engine.detectors
t1 = time.perf_counter()
import zspace_engine.auditors
t2 = time.perf_counter()
import zspace_engine.context
t3 = time.perf_counter()
import zspace_engine.validator
t4 = time.perf_counter()
import zspace_engine.report
t5 = time.perf_counter()
import zspace_engine.ingestion
t6 = time.perf_counter()
print(f"detectors {t1-t0:.2f}s auditors {t2-t1:.2f}s context {t3-t2:.2f}s validator {t4-t3:.2f}s report {t5-t4:.2f}s ingestion {t6-t5:.2f}s")
"""
t0 = time.perf_counter()
r = subprocess.run(["python3", "-c", code2], capture_output=True, text=True, timeout=300, cwd="/mnt/d/Axioms/Axiom-Zspace-CODE")
print(r.stdout)
print(r.stderr[-500:] if r.returncode else "")