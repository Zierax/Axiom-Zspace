<div align="center">

# Axiom-ZSpace

**Blind-search exoplanet transit detection for TESS/Kepler light curves**

BLS period search with a period-prior ladder · physical-invariant auditing ·
a false-positive ruling engine with a circuit breaker · a **measured,
single-source threshold catalog**

<sub>v1.0 — public release · 101-test suite · controlled + real benchmarks as
first-class evidence · [documentation hub →](docs/README.md)</sub>

</div>

---

## Overview

Axiom-ZSpace runs a light curve through a single deterministic pipeline:

```
ingestion → BLS detection (period-prior ladder, harmonic rejection, FAP
firewall) → ephemeris resolution → transit-physics audits → context/density
checks → false-positive gate engine (circuit breaker) → CVS classification →
discovery card with a full proof chain
```

Every tunable number lives in **one place** —
`config/production.yaml` + `zspace_engine/thresholds.py` — and every number is
**measured**: threshold changes must be re-benchmarked before shipping
(the one rule, enforced in [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md)).

## Features

- **Blind search with a prior, not a lock**: the period hint (benchmarks)
  is only *preferred* when its harmonic-ladder peak carries enough power;
  otherwise the global BLS peak wins.
- **Physical-invariant auditing**: even/odd depth (Welch t-test), U/V shape
  ratio, per-transit depth consistency, ingress/egress, secondary eclipse,
  stellar-density mismatch, centroid (optional), MCMC posterior (optional).
- **Evidence-bearing verdicts**: every card ships a human-readable proof
  chain; no gate result is decided silently.
- **Circuit breaker**: a failed critical gate makes `SOVEREIGN_PASS`
  impossible — asserted by tests.
- **One measured catalog**: three profiles
  (`conservative` / `balanced` / `sensitive`), auto-generated reference report
  (`THRESHOLDS_REPORT.md`).
- **Determinism as a contract**: same seed ⇒ same results, asserted by
  `tests/test_reproducibility.py`.
- **Dual sovereign engine**: the verdict engine also ships as a
  **dependency-free C99 binary** whose math kernels are machine-generated
  from a strict numpy subset by [**Purce**](https://github.com/Zierax/Purce),
  differentially verified against the Python reference (148/148 kernels,
  90/90 cards) — selectable with `--engine {python,c99}` (`python` default for reference, `c99` **recommended for batch**: `46 ms/TIC` light `604×` / `4.8 s/TIC` heavy `5.8×`, `docs/BENCHMARKS.md:151`).
- **Offline-first**: the engine and the test suite run without network; the
  real benchmark rehydrates from pinned snapshot JSONs + a one-time MAST
  download.

## Measured results (balanced profile, 2026-08-16/17)

| Benchmark | Result |
|---|---|
| **Controlled BIG400 (400 true + 400 false, 8 seeds)** | recall **41.2%** (165/400) · contamination FPR **4.25%** (17/400) · wrong-ephemeris 20 |
| Controlled, fixed seed 20260814 (100+80) | recall 52/100 · **0/80 (0%) contamination FPR** (one-sample — see note) |
| Real Kepler (12 true hosts / 12 quiet stars): recall@target | **41.7% (5/12)** — ≤0.01% period accuracy when matched |
| Real Kepler: quiet-star certification (proxy FPR) | 33.3% (4/12) — honest interpretation in report §2.5 |
| Coherent override (FP-2 bypass) | **OFF by default** — probes measured contamination FPR up to 62.5% when ON |

### Engine speed — Python vs C99 (measured 2026-08-21, 16 cores, `-O3 -march=native -flto -fopenmp`, `C99-Version/bin/zspace_card`)

| Dataset | n_points | Python `run_controlled` | C99 `bin/zspace_card batch` | Speedup | Verdict parity |
|---|---|---:|---:|---:|---|
| Controlled 100-light (syn 3k) | ~3k | 27.8 s / TIC | **42.8 ms / TIC** (`100 in 4.28s`) | **604×** | 400/400 identical |
| Heavy 90k (5-sector 2-min, ~87k) | ~87k | 27.8 s / TIC | **4.8 s / TIC** (`10 in 48s`) | **5.8×** | 8/10* |
| `verify_compare` | — | — | — | — | **148/148** kernels |
| `parity_card --n 90 --lc` | — | — | — | — | **90/90** cards |

\* Heavy 8/10 agreement — 2 marginal `FAP≈0.05` flips (syn_0 `fap 0.008→0.078`, syn_2 alias), both `O(n·n_freq)` bound. Light is the benchmark for `1000×` target. **Recommendation:** `python` default for reference/single, **`c99` recommended for batch** (`--engine c99`) — `docs/BENCHMARKS.md:151`, `docs/C99_ENGINE.md:18`.

> **FPR honesty:** 0/80 is the fixed-seed sample; across 8 fresh samples
> (BIG400, 800 targets) the measured contamination FPR is **4.25%**, dominated
> by high-SNR eclipsing binaries that pass every gate (pure noise: 0.6%).
> Both numbers are committed evidence — see `docs/BENCHMARKS.md` §1.6.

Evidence runs are versioned in the repo:
`benchmarks_controlled/evidence/{OVERHARM_FIX2,BIG400,PROBE_FPR68*}` and
`benchmarks_real/evidence/REAL_FINAL` — per-target JSONs, recomputable
metrics, and the reports that derive from them. How to reproduce everything
in two commands: [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md).

> The real numbers are a **measurement, not a tune** — re-run both benchmarks
> after any threshold edit.

## Quick start

```bash
pip install -r requirements.txt          # Python ≥ 3.10

python run_pipeline.py --synthetic       # offline self-test (no internet)
python run_pipeline.py --synthetic --engine c99   # same, C99 sovereign engine
python -m pytest tests/ -q               # 101 tests, ~5 min, offline-green

python run_pipeline.py --tic 260128333   # single TIC, end-to-end
python run_pipeline.py --sector 42       # whole sector (heavy; see PRODUCTION.md)
```

CLI: `--tic ID` · `--synthetic` · `--engine {python,c99}` · `--sector N` /
`--sectors A-B` · `--max-targets N` · `--output DIR` · `--mcmc` ·
`--tpf-centroids` · `--multi-sector` (full reference:
[`docs/QUICKSTART.md`](docs/QUICKSTART.md)).

## Documentation

| Doc | What it gives you |
|---|---|
| [`docs/README.md`](docs/README.md) | the hub + where every number lives |
| [`docs/QUICKSTART.md`](docs/QUICKSTART.md) | install, run modes, troubleshooting |
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | stage-by-stage design, gates, verdicts, CVS, invariants |
| [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md) | methodology + measured results + re-measure procedure |
| [`docs/TESTING.md`](docs/TESTING.md) | the 101-test suite, determinism contract |
| [`docs/PRODUCTION.md`](docs/PRODUCTION.md) | batch operation, config, measured-vs-unmeasured scope |
| [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) | the one rule + claim discipline + review checklist |
| [`docs/GLOSSARY.md`](docs/GLOSSARY.md) | every acronym and metric |
| [`docs/C99_ENGINE.md`](docs/C99_ENGINE.md) | the C99 sovereign engine (Purce, verification, build) |
| `THRESHOLDS_REPORT.md` (root) | every threshold value, its evidence, pros/cons — auto-generated |
| `CHANGELOG.txt` | release history (v2.x entries are labeled LEGACY — history, not evidence) |

## Repository layout

```
run_pipeline.py          CLI entry point (+ c99_bridge.py — Python↔C99 adapter)
zspace_engine/           the engine: ingestion, detectors, ephemeris,
                         auditors, context, validator, core, report,
                         thresholds, sector processor, output organizer
C99-Version/             the C99 sovereign engine (Purce-generated kernels,
                         zspace_core, zspace_card CLI, verification harnesses)
benchmarks_controlled/   synthetic benchmark (injector + runner) + evidence
benchmarks_real/         real Kepler benchmark (pinned NEA snapshots) + evidence
tests/                   the 101-test suite (+ validator demo helper)
scripts/                 snapshot fetcher, sector-aggregation utility
scripts/legacy/          pre-v1 tools (C engine, old harnesses) — NOT part of v1
config/                  canonical configuration (threshold catalog)
docs/                    this documentation set
archive/                 local-only pre-v1 history — git-ignored, never pushed
```

## Honesty & reproducibility

- **Evidence over claims.** Every number in the docs traces to a committed
  artifact; the CHANGELOG's v2.x entries are explicitly labeled legacy
  because their batch-scale numbers were not preserved. Anything else is not
  claimed.
- **Deterministic samples.** The controlled benchmark is seed-fixed
  (`--seed 20260814`); the real sample is pinned by offline NEA snapshot JSONs
  (`scripts/fetch_nea_snapshot.py`). Two identical commands produce identical
  results — asserted by tests, verifiable by any clone.
- **Privacy by construction.** `archive/`, run outputs, caches, discovery
  cards, credential-bearing scripts and keys are git-ignored by hard rule.

## Contributing & license

Read [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) first — especially the
"one rule": *no threshold change ships without both benchmarks re-measured*,
and *no number ships without an artifact*.

---

<sub>Built for blind, honest, reproducible transit searching. If the numbers
surprise you, the evidence tells you why — that is the point.</sub>