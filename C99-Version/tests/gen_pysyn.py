"""Generate Python synthetic targets (true planets) as CSVs + manifest for C batch."""
import sys
sys.path.insert(0, '/mnt/d/Axioms/Axiom-Zspace-CODE')
import numpy as np
from benchmarks_controlled.synthetic import generate_true_planet

n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
outdir = '/tmp/pysyn'
import os
os.makedirs(outdir, exist_ok=True)
manifest = []
for i in range(n):
    period = 2.0 + (i * 0.9) % 10.0
    target = generate_true_planet(i, period, target_snr=10.0, seed=20260816)
    path = f'{outdir}/syn_{i}.csv'
    np.savetxt(path, np.column_stack([target.time, target.flux]),
               delimiter=',', header='time,flux', comments='')
    manifest.append(path)
    print(f'{i}: P={target.label_period:.3f} depth={target.injected_depth:.5f} pts={len(target.time)}')

with open(f'{outdir}/manifest.txt', 'w') as f:
    f.write('\n'.join(manifest) + '\n')
print(f'manifest: {outdir}/manifest.txt')