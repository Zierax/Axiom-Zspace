import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
print("HERE", HERE)
print("ROOT", ROOT)
print("exists", os.path.exists(os.path.join(ROOT, "run_pipeline.py")))
print("sys.path", sys.path[:3])
sys.path.insert(0, ROOT)
import run_pipeline
print("import ok")
