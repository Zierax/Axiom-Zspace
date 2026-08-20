import csv
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
import numpy as np

rows = []
with open("C99-Version/tests/tmp_lc.csv") as fh:
    reader = csv.reader(fh)
    next(reader)
    for r in reader:
        rows.append((float(r[0]), float(r[1])))
t = np.asarray([r[0] for r in rows])
f = np.asarray([r[1] for r in rows])
print("n", len(rows), "t0", t[0], "t1", t[-1])
from zspace_engine.auditors import TransitAuditor
aud = TransitAuditor(run_mcmc=False)
eo = aud.even_odd_test(t, f, 3.695262, 2459300.0 + 1.2, 4.0 / 24.0)
print("eo ok", eo.delta_sigma)