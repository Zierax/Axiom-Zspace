# Axiom-ZSpace Documentation

The documentation hub. Read in this order for a full picture: **Quickstart →
Architecture → Benchmarks → Testing → Production → Contributing**. Use the
Glossary for terminology, THRESHOLDS_REPORT.md for every tunable number, and
CHANGELOG.txt for release history (v2.x entries are legacy, not evidence).

## Reading map

| Doc | Audience | Answers |
|---|---|---|
| [`QUICKSTART.md`](QUICKSTART.md) | everyone | install, run a TIC, run the suite, benchmark commands, troubleshooting |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | developers | how the pipeline works stage-by-stage, module register, gate engine, CVS, invariants |
| [`BENCHMARKS.md`](BENCHMARKS.md) | everyone | the two benchmarks, definitions, measured numbers, re-measure procedure |
| [`TESTING.md`](TESTING.md) | developers | the 101-test suite, determinism contract, adding tests |
| [`PRODUCTION.md`](PRODUCTION.md) | operators | batch/sector operation, config knobs, measured-vs-unmeasured scope |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | contributors | the one rule (measure before you tune), claim discipline, review checklist |
| [`GLOSSARY.md`](GLOSSARY.md) | everyone | every acronym and metric defined |
| [`VALIDATION_IMPROVEMENTS.md`](VALIDATION_IMPROVEMENTS.md) | historians | pre-v1 validation-fix writeup (legacy history) |

## Generated / root documents

| Doc | Origin |
|---|---|
| `THRESHOLDS_REPORT.md` (root) | auto-generated from live config via `python -m zspace_engine.thresholds_report` — the reference for every threshold value and its measured evidence |
| `README.md` (root) | the public face — start there |
| `CHANGELOG.txt` (root) | release history; entries marked LEGACY are history, not evidence |

## Where the numbers live

```
evidence (versioned):    benchmarks_controlled/evidence/{OVERHARM_FIX2,PROBE_FPR68*}
                         benchmarks_real/evidence/REAL_FINAL
threshold reference:     THRESHOLDS_REPORT.md            (generated)
config source of truth:  config/production.yaml  →  thresholds.py
tests (executed truth):  tests/ (101 tests, offline-green)
snapshot cache:          benchmarks_real/data/*.json    (regenerate via scripts/fetch_nea_snapshot.py)
```

## Quick reference — the commands that matter

```bash
pip install -r requirements.txt
python run_pipeline.py --synthetic                    # offline self-test
python run_pipeline.py --tic 260128333                # single target
python -m pytest tests/ -q                            # 101 tests
python benchmarks_controlled/run_controlled.py --true 100 --false 80 \
    --out benchmarks_controlled/runs/MY --seed 20260814
python benchmarks_real/run_real.py --n-true 12 --n-false 12 \
    --out benchmarks_real/runs/MY
python -m zspace_engine.thresholds_report             # regenerate reference report
```