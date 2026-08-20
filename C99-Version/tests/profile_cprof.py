"""cProfile hotspot identification for run_synthetic_test."""
import cProfile
import os
import pstats
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("ZSPACE_SKIP_LOGS", "1")

from run_pipeline import run_synthetic_test

prof = cProfile.Profile()
prof.enable()
res = run_synthetic_test(engine="python")
prof.disable()

st = pstats.Stats(prof)
st.sort_stats("cumulative")
st.print_stats(18)