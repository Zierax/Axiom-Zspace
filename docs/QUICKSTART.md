# Quickstart

Everything you need to go from a fresh clone to a working pipeline, then to a
detection. For the full mental model see `docs/ARCHITECTURE.md`; for measured
numbers and how to reproduce them see `docs/BENCHMARKS.md`.

## 1. Requirements

| Item | Requirement |
|---|---|
| Python | ≥ 3.10 |
| Core stack | numpy, scipy, astropy |
| Optional (network features) | lightkurve, astroquery |
| Benchmark-only | batman-package (transit model), pandas/matplotlib |

`requirements.txt` covers all of it, including `pytest` for the test suite.

## 2. Install

```bash
python -m pip install -r requirements.txt
```

> Tip: use a virtual environment.
> On Windows consoles set `PYTHONIOENCODING=utf-8` if exotic characters
> (`σ`, `Δ`) ever render wrong; the CLI now forces UTF-8 itself.

## 3. Provenance of every number

Before doing anything: the threshold catalog in `config/production.yaml` and
`zspace_engine/thresholds.py` is the single source of truth; the reference
report `THRESHOLDS_REPORT.md` is generated from it:

```bash
python -m zspace_engine.thresholds          # show the active-profile table
python -m zspace_engine.thresholds --profile sensitive   # a specific profile
python -m zspace_engine.thresholds_report   # regenerate THRESHOLDS_REPORT.md
```

## 4. Run modes

### 4.1 Offline self-test (no internet, nothing to download)

```bash
python run_pipeline.py --synthetic
```

Builds pure-synthetic light curves with an injected transit and runs the full
blind pipeline end-to-end (detection → audits → validation → discovery card).
Seconds to a minute. Second sanity:

```bash
python -m pytest tests/test_reproducibility.py tests/test_synthetic.py -q
```

### 4.2 Single TIC (TESS)

```bash
python run_pipeline.py --tic 260128333
```

Downloads the light curve (MAST), caches it locally, runs detection + audits +
validation, writes the discovery card under the output directory (default
`axiom_output/`). Add optional rigor:

```bash
python run_pipeline.py --tic 260128333 --mcmc --tpf-centroids --multi-sector
```

(`--mcmc` = full posterior transit fit; `--tpf-centroids` = pixel-level
centroid; `--multi-sector` = cross-sector consistency. All three are
**off by default** — they are expensive and only `multi-sector` is cheap at
single-target scale.)

### 4.3 Whole sector / many sectors

```bash
python run_pipeline.py --sector 42
python run_pipeline.py --sectors 1-50 --max-targets 200   # capped, test-sized
```

Batch mode dispatches through `SectorProcessor`. **Heavy**: hours plus
GB-scale FITS caches. Batch mode is NOT part of the measured benchmark claims
(v1.0 measured single-target engine behavior only — see
`docs/PRODUCTION.md`).

### 4.4 Full CLI reference

```
--tic TIC            process one TIC ID
--synthetic          run the offline synthetic self-test
--sector N           process a whole TESS sector
--sectors A-B        process a range of sectors
--max-targets N      cap targets per sector (testing)
--output DIR         override the output base directory
--mcmc               enable MCMC posterior validation
--tpf-centroids      enable TPF centroid analysis
--multi-sector       enable multi-sector consistency
```

## 5. Test suite

```bash
python -m pytest tests/ -q        # 101 tests
```

- Green both offline and with network (network-dependent tests mock or
  degrade gracefully).
- Rough runtime: ~5 minutes.
- Map of what every test file covers: `docs/TESTING.md`.

## 6. Benchmarks (reproducibility in two commands)

```bash
# controlled — offline, ~minutes, deterministic with --seed
python benchmarks_controlled/run_controlled.py --true 100 --false 80 \
    --out benchmarks_controlled/runs/MY_RUN --seed 20260814

# real Kepler — network-bound first run (~2 h), cached afterwards
python benchmarks_real/run_real.py --n-true 12 --n-false 12 \
    --out benchmarks_real/runs/MY_REAL
```

Measured values (BIG400 800 targets: 41.2% recall, 4.25% contamination FPR;
fixed-seed suite: 52/100 recall, 0/80 FPR; real: 41.7% recall@target, 33.3%
quiet-star proxy FPR) and what they mean: `docs/BENCHMARKS.md`.

## 7. Reproducing the real benchmark sample

`benchmarks_real/data/*.json` pin the sample. To refresh after NEA changes:

```bash
python scripts/fetch_nea_snapshot.py --out benchmarks_real/data
```

This rewrites `real_ps_star.json`, `real_quiet.json`, and a
`_snapshot_meta.json` sidecar. `real_known_signals.json` (the truth keys the
runner caches) is managed by the runner itself.

## 8. Configuration

- `config/production.yaml` — the canonical config (threshold catalog +
  network/cache/output settings).
- Environment override: `AXIOM_CONFIG=/path/to/config.yaml`.
- Cache dir: `.cache/fits` (see `fits_cache:` section; size-capped in config).

## 9. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `ImportError: No module named 'zspace_engine'` | run from the repository root (the package is importable in place); check `pip install -r requirements.txt` |
| Tests fail on `astroquery`/SIMBAD logs | those are warnings, not failures; network-dependent tests pass offline by design |
| `--tic` run hangs on download | network-bound; retry or check `network:` timeouts in config |
| Extra `axiom_pipeline.log` appears at the repo root | created by the logging subsystem on first import — it is git-ignored, safe to delete |
| `thresholds --show` prints a `runpy RuntimeWarning` | benign duplicate-module warning of the CLI entry; ignore |
| Benchmark numbers drift from docs | you changed thresholds — re-measure and update the docs (the one rule, §6 above) |
| `UnicodeEncodeError` from CLI tools | the CLI forces UTF-8; on very old consoles use `PYTHONIOENCODING=utf-8` |

## 10. Repository layout (one glance)

```
run_pipeline.py          CLI entry (single TIC / synthetic / sector)
zspace_engine/           the engine (detectors, auditors, validator, …)
benchmarks_controlled/   synthetic benchmark + evidence runs
benchmarks_real/         real Kepler benchmark + evidence runs
scripts/                 snapshot fetch, sector-aggregation utilities
scripts/legacy/          pre-v1 tools (C engine, old bench harnesses) — NOT v1
tests/                   the 101-test suite
config/                  canonical configuration
docs/                    this documentation hub
archive/                 local-only pre-v1 history — NEVER pushed
```