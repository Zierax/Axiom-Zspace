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
- **Offline-first**: the engine and the test suite run without network; the
  real benchmark rehydrates from pinned snapshot JSONs + a one-time MAST
  download.

## Measured results (balanced profile, 2026-08-16)

| Benchmark | Result |
|---|---|
| Controlled: recall, 100-target suite | **52/100** (first-80 prefix 40/80; 32/80 pre-calibration baseline, archived) |
| Controlled: contamination FPR | **0/80 (0%)** — all 80 false targets rejected |
| Controlled: wrong-ephemeris certs (true set) | 1 (counted separately) |
| Real Kepler (12 true hosts / 12 quiet stars): recall@target | **41.7% (5/12)** — ≤0.01% period accuracy when matched |
| Real Kepler: quiet-star certification (proxy FPR) | 33.3% (4/12) — honest interpretation in report §2.5 |
| Coherent override (FP-2 bypass) | **OFF by default** — probes measured contamination FPR up to 62.5% when ON |

Evidence runs are versioned in the repo:
`benchmarks_controlled/evidence/OVERHARM_FIX2`, `PROBE_FPR68*`, and
`benchmarks_real/evidence/REAL_FINAL` — per-target JSONs, recomputable
metrics, and the reports that derive from them. How to reproduce everything
in two commands: [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md).

> The real numbers are a **measurement, not a tune** — re-run both benchmarks
> after any threshold edit.

## Quick start

```bash
pip install -r requirements.txt          # Python ≥ 3.10

python run_pipeline.py --synthetic       # offline self-test (no internet)
python -m pytest tests/ -q               # 101 tests, ~5 min, offline-green

python run_pipeline.py --tic 260128333   # single TIC, end-to-end
python run_pipeline.py --sector 42       # whole sector (heavy; see PRODUCTION.md)
```

CLI: `--tic ID` · `--synthetic` · `--sector N` / `--sectors A-B` ·
`--max-targets N` · `--output DIR` · `--mcmc` · `--tpf-centroids` ·
`--multi-sector` (full reference: [`docs/QUICKSTART.md`](docs/QUICKSTART.md)).

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
| `THRESHOLDS_REPORT.md` (root) | every threshold value, its evidence, pros/cons — auto-generated |
| `CHANGELOG.txt` | release history (v2.x entries are labeled LEGACY — history, not evidence) |

## Repository layout

```
run_pipeline.py          CLI entry point
zspace_engine/           the engine: ingestion, detectors, ephemeris,
                         auditors, context, validator, core, report,
                         thresholds, sector processor, output organizer
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