# Axiom-ZSpace v1

**White-box exoplanet transit detection for TESS/Kepler light curves — with a
measured 0/80 contamination false-positive rate.**

This is the first public release of Axiom-ZSpace: a blind-search transit
detection pipeline built to be *audited*. No black box — every verdict ships
with a human-readable proof chain explaining exactly why a candidate was
certified or rejected, and every threshold in the pipeline is a measured,
benchmarked value.

## Why this is different from most existing solutions

| Most black-box detectors | Axiom-ZSpace v1 |
|---|---|
| Verdicts come from a learned model — you cannot see *why* | **White-box proof engine**: 10 open false-positive gates (FP-1…FP-10), each with a threshold, a direction, and a purpose; every card carries the full proof chain |
| FPR claims come from a paper's chosen test set | **Measured in-repo**: 0/80 (0.0%) contamination FPR on the controlled benchmark, with the per-target JSONs committed as evidence |
| Thresholds are scattered constants tuned ad hoc | **One catalog** (`config/production.yaml` + `zspace_engine/thresholds.py`), every value documented with measured evidence and pros/cons in the auto-generated THRESHOLDS_REPORT.md |
| Reproducibility is "trust us" | **Determinism is a tested contract**: same seed ⇒ same results, asserted by 7 dedicated tests |
| Science-grade sanity is optional | **Physical-invariant auditing** as a first-class stage: even/odd depth, U/V shape, ingress/egress, secondary eclipse, stellar-density consistency, circuit breaker |
| Numbers that disappear after publication | **Evidence is versioned with the code** — clones can recompute every metric themselves |

## Measured results (balanced profile, 2026-08-16)

| Benchmark | Result |
|---|---|
| Controlled: recall on 100-target suite | 52/100 |
| Controlled: contamination FPR | **0/80 (0.0%)** |
| Real Kepler: recall@target | 41.7% (5/12), ≤0.01% period accuracy when matched |
| Real Kepler: quiet-star certification (proxy FPR) | 33.3% (4/12) — reported honestly, with the full caveat list |
| Coherent override of the FAP firewall | OFF by default — probes measured up to 62.5% contamination FPR when ON |

Honest numbers, honest limits: the real-data values are a *measurement of the
pipeline*, not a tune; quiet stars are a proxy false set; recall is
precision-first (52/100 controlled), not a complete census. Read
[`docs/BENCHMARKS.md`](docs/BENCHMARKS.md) before quoting any of this.

## Quick start

```bash
pip install -r requirements.txt          # Python ≥ 3.10

python run_pipeline.py --synthetic       # offline self-test, no internet
python run_pipeline.py --tic 260128333   # single TESS target, end-to-end
python -m pytest tests/ -q               # 101 tests, offline-green

# reproduce the measured benchmarks (docs/BENCHMARKS.md):
python benchmarks_controlled/run_controlled.py --true 100 --false 80 \
    --out benchmarks_controlled/runs/MY_RUN --seed 20260814
python benchmarks_real/run_real.py --n-true 12 --n-false 12 \
    --out benchmarks_real/runs/MY_REAL
```

## What's in this release

- Blind BLS detection with a period-prior ladder (prior, not lock), harmonic
  rejection, and a self-calibrating FAP firewall.
- Ephemeris resolution + identity gate (wrong ephemerides are demoted, never
  certified).
- 10-gate false-positive ruling engine with a circuit-breaker invariant.
- CVS classification (4 tiers) with measured weights (0.97/0.83/0.61/0.31).
- Three threshold profiles; the `sensitive` profile is explicitly
  experimental and unmeasured.
- 101-test suite (determinism, calibration, physics audits, circuit breaker).
- Full documentation: architecture, benchmarks, testing, production,
  contributing, glossary, quickstart.
- Repository hygiene: secrets, `archive/`, run outputs and caches are
  git-ignored by hard rule; unverifiable legacy claims removed from the
  changelog.

## Links

- Docs hub: [`docs/README.md`](docs/README.md)
- Threshold reference: [`THRESHOLDS_REPORT.md`](THRESHOLDS_REPORT.md)
- Benchmark evidence (committed): `benchmarks_controlled/evidence/`,
  `benchmarks_real/evidence/`
- Changelog: [`CHANGELOG.txt`](CHANGELOG.txt)
