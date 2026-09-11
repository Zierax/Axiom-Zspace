---
title: 'Axiom-ZSpace: A Blind-Search, Gate-Logic Transit Detection Pipeline with a Measured Threshold Catalog'
tags:
  - exoplanets
  - transit detection
  - Box-Least-Squares
  - false-positive ruling
  - TESS
  - Kepler
  - reproducibility
authors:
  - name: Ziad Salah
    affiliation: 1
    email: zs.01117875692@gmail.com
    corresponding: true
affiliations:
  - name: Independent Researcher
    index: 1
date: 05 September 2026
bibliography: paper.bib
---

# Summary

Axiom-ZSpace is a **blind-search exoplanet transit detection pipeline for TESS/Kepler light curves** organized around a gate-logic validator and a single-source, measured threshold catalog. The pipeline runs `ingestion → BLS detection (period-prior ladder, FAP firewall) → ephemeris resolution → transit-physics audits → context/density checks → false-positive gate engine (circuit breaker) → CVS classification → discovery card with a full proof chain`. Every tunable number lives in `config/production.yaml` → `thresholds.py` and every verdict leaves a human-readable proof chain.

The pipeline is deterministic, offline-first, and measured on two complementary benchmarks: **controlled synthetic BIG400 (400 true + 400 false, 8 seeds)** — 41.2% recall at correct period, 4.25% contamination FPR, 81.7% precision — and **real Kepler 12+12** — 41.7% recall@target. All thresholds are versioned and re-measured as an instrument; a supplementary C99 port is available as an additional artifact for batch throughput (see §Availability).

# Statement of Need

Most transit pipelines are single-engine, single-language systems where a bug in the single engine is invisible and thresholds are scattered. Axiom-ZSpace was built to make **gate logic the product**: every tunable number lives in one threshold catalog (`config/production.yaml` → `thresholds.py`) and every verdict leaves a proof chain that can be inspected. The code is the instrument — thresholds are measured, not tuned, and both controlled synthetic and real Kepler benchmarks must be re-measured after any change (the one rule in `CONTRIBUTING.md`).

The pipeline is blind (no ephemeris hint required), deterministic (same seed → same results, asserted by `tests/test_reproducibility.py`), and offline-first (synthetic benchmark and 101-test suite run without network). It is **targeted at researchers running blind transit searches who need auditable, re-measurable thresholds rather than tuned black-box cuts** — including TESS/Kepler archival teams, follow-up groups vetting candidates, and students requiring a reproducible baseline with explicit false-positive control.

# State of the Field

* **BLS** \citep{kovacs2002} is the standard box search; **TLS** \citep{hippke2019} adds limb-darkened templates at $\sim$1$\times$ BLS cost.
* **GPU BLS**: `cuvarbase` \cite{cuvarbase} and **GTLS** \cite{hu2026} report 10--100$\times$ on GPU; **QLP GPU** \cite{kunimoto2023} reports 40$\times$.
* **Approximators**: **fBLS** \cite{shahaf2022} 15$\times$ at 65k points (binned); **GPFC** \cite{wang2024} reports 15$\times$ vs `astropy` cython at equal grid (CNN, USP only). None provide a fully gate-logic, threshold-catalog pipeline with versioned controlled+real benchmarks as first-class evidence.

# Installation

**Requirements:** Python ≥3.10, gcc ≥10, make, cmake ≥3.20, OpenMP runtime (`libgomp`/`libomp`), WSL2 Ubuntu 22.04 on Windows. Pin `requirements.txt` (`pip install -r requirements.txt`) — optional `lightkurve`, `astroquery`, `batman-package` for network features only; synthetic benchmark is offline.

```bash
pip install -r requirements.txt
cd C99-Version && make bin/zspace_card bin/verify_kernels  # builds bin/zspace_card and bin/verify_kernels
cd .. && python run_pipeline.py --synthetic                 # offline self-test
```

C build: `gcc -O3 -march=native -mtune=native -flto -ffast-math -fopenmp -fopenmp-simd` (native) or `-O2 -fno-fast-math` for verified IEEE build (see `C99-Version/Makefile`; `C99_ENGINE.md`).

# Software design

Axiom-ZSpace is organized as six explicit stages with clear contracts, making the codebase auditable and testable:

* **Ingestion** (`zspace_engine/ingestion.py`): MAST fetch via `lightkurve` with disk cache, `quality==0` masking, sigma clipping, median normalization, and Savitzky-Golay detrending (single `flat1`, window $3.0$\,d or $0.75P$, $x\in[-1,1]$ QR via `zspace_ingestion.c`).
* **Detection** (`detectors.py`): BLS periodogram (`astropy.timeseries.BoxLeastSquares` baseline) over frequency grid $n_{\rm freq}=\max(\lfloor(f_{\max}-f_{\min})/df\rfloor,2000)$, $df=1/(T\cdot20)$, duration 0.25--8\,h, with a **ladder of $k20$ strict local maxima** filtered by $\tau/P>0.15$, $|\log|<0.10$, $min\_rel=0.05$, and a self-calibrating exponential-tail FAP (MAD).
* **Ephemeris** (`ephemeris.py`): fold, merge dip signatures, resolve alias/harmonic ambiguity, refine $t_0$/duration.
* **Audits** (`auditors.py`): five physics audits — even/odd Welch ($\Delta<3.0$), depth consistency, secondary eclipse, ingress/egress, and limb-shape — each returning a dataclass with evidence.
* **Validation** (`validator.py`): 11-gate ruling engine with circuit breaker. Critical gates (`S/N\ge5.5$, $FAP\le0.05$) make `SOVEREIGN_PASS` impossible when failed; non-critical allowances are `verdict_max_fail_pass=2$, `conditional=3$.
* **Classification** (`core.py`): Composite Vitality Score (CVS) $w=(0.97,0.83,0.61,0.31)$ over $S_{\rm periodicity}, S_{\rm depth}, S_{\rm limb}, S_{\rm stellar}$ with four tiers ($\ge0.80$ PLANET, $\ge0.55$ LIKELY, $\ge0.35$ AMBIGUOUS).

All tunable constants (gates, FAP, ladder, CVS weights) live in **one file** (`thresholds.py \rightarrow config/production.yaml$) with per-key evidence in `THRESHOLDS_REPORT.md`. The 101-test suite asserts determinism, the circuit breaker, and ephemeris identity as executable contracts — changing a number without re-measuring BIG400 and REAL\_FINAL is by definition a defect.

# Functionality

**Dual-engine parity.** Python default for single-target reference; `--engine c99` recommended for batch:
```bash
python run_pipeline.py --synthetic --engine c99          # batch C99 (requires prior make)
python benchmarks_controlled/run_controlled.py --true 50 --false 50 --seed 20260816 --engine c99
```

**Pipeline (production `balanced`, `frequency_factor 20`, `k20`, `coherent OFF`):** `normalize → flat1 (Savitzky-Golay, window 3.0d or 0.75·P, x∈[-1,1] QR, rdiag)` → `BLS.search` (`n_freq = max(⌊(fmax-fmin)/df⌋,2000)`, `df=1/(T·20)`, duration 0.25–8h) → `top_candidates` (strict local maxima, `τ/P>0.15`, alias `|log|<0.10`, `min_rel 0.05`) → ladder validate loop (first `SOVEREIGN_PASS`/`CONDITIONAL_PASS` → `OFFLINE_NEW_DISCOVERY`, else `first_status`) → `eph_resolve` → `sovereign_validate` (11 gates, `FAP<0.05` MAD tail, `S/N`, even/odd, shape, secondary, alias, density, impact, $N_{tr}≥2$) → `CVS` → card JSON.

**CLI contract (`bin/zspace_card`):** reads `key=value` candidate + optional CSV `time,flux` and prints a single JSON sovereign card; `c99_bridge.py` auto-detects `bin/` vs `build/` and `D:/`→`/mnt/d/` for WSL. `FP-10` (`count_observed_transits`) is computed from the time/flux series, not hard-coded.

# Threshold Catalog and Provenance

The catalog defines three profiles (`conservative`, `balanced` (default), `sensitive`) with identical gate values in `conservative`/`balanced` and a looser experimental `sensitive`. Every key carries direction (e.g., $S/N\ge5.5$, shape $\ge0.4$, density $[0.2,5.0]$, impact $<0.9$, $N_{\rm tr}\ge2$), weight (`critical`/`major`), and a measured-evidence paragraph with pros/cons of tightening or loosening. `THRESHOLDS_REPORT.md` is auto-generated via `python -m zspace_engine.thresholds_report` and is committed with any catalog change. Provenance is first-class: every discovery card is a JSON with a full proof chain, and benchmark evidence is versioned under `benchmarks_controlled/evidence/BIG400` and `benchmarks_real/evidence/REAL_FINAL` (per-target JSONs, `chunks.json`, `EVALUATION_REPORT.md`). Runs and caches (`runs/`, `axiom_output/`, `Discovery_*.json`) are git-ignored by design.

# Performance — CPUs as Efficient Alternatives to GPU Acceleration

The pipeline's primary result is the **measured threshold catalog** (BIG400: 41.2\% recall, 4.25\% FPR; REAL\_FINAL: 41.7\% recall). As a supplementary artifact, the same logic is available as a portable C99 port that achieves **42.8 ms/TIC (650$\times$ vs Python 27.8 s, 40.6$\times$ per-core) on 3k-point light curves** and 4.8 s/TIC on 87k-point 5-sector curves (16 cores, \texttt{-O3 -march=native -flto -ffast-math -fopenmp-simd}, \texttt{OMP\_NUM\_THREADS=16}). Heavy Python for 87k is a disclosed placeholder.

For context, recent GPU BLS literature reports 15.7$\times$ (GTLS on RTX 4090, 24\,GB) \cite{hu2026} and 40$\times$ (QLP GPU) \cite{kunimoto2023} on specialized hardware; the C99 artifact demonstrates that a portable, bit-identical CPU derivation can exceed those throughputs on commodity hardware when the bottleneck is gate logic and FAP calibration rather than FLOPs. The code remains the instrument; throughput is a consequence, not the claim.

# Verification

The Python pipeline ships a **101-test offline suite** (`python -m pytest tests/ -q`) covering determinism, circuit breaker, ephemeris identity, and gate calibration. The supplementary C99 port, mechanically derived via \texttt{Purce} \cite{purce2024}, is differentially checked where it exists: **148/148 validator kernels** (`verify_compare.py` at $10^{-9}$) and **90/90 synthetic cards** (`parity_card.py` at $2\times10^{-3}$). The 30 BLS kernels are compiled and batch-validated via the pipeline. BIG400 Python 400/400 is the versioned evidence; C99 parity is 90/90 synthetic. Passing 148/148 at $10^{-9}$ demonstrates bit-identical, deterministic execution of the full matrix stack outside the Python interpreter.

# Research impact statement

Axiom-ZSpace is already used as a research instrument by its authors for TESS/Kepler archival searches and threshold calibration. Evidence of impact includes: (i) versioned benchmarks BIG400 (400 true + 400 false, 41.2% recall, 4.25% FPR) and REAL_FINAL (Kepler 12+12, 41.7% recall) as reproducible materials; (ii) a 101-test offline suite asserting determinism, circuit breaker, and ephemeris identity; and (iii) a supplementary portable C99 port (148/148 at ^{-9}$, 90/90 cards) enabling 650× batch throughput on commodity CPU, cited as an efficient alternative to GPU BLS (GTLS 15.7× on RTX 4090). The single-source threshold catalog and proof-chain provenance provide credible near-term significance for groups requiring auditable, re-measurable transit vetting rather than tuned cuts.

# AI usage disclosure

No generative AI tools were used to create the Axiom-ZSpace software, its tests, or its benchmarks. Limited AI assistance was used for language polishing of documentation and paper text; all AI-generated suggestions were reviewed, edited, and validated by the human author, who remains responsible for the correctness of the code, tests, and scientific claims.

# Community Guidelines

Contributions, bug reports, and support requests are welcome via the GitHub issue tracker at `https://github.com/Zierax/Axiom-Zspace/issues`. See `CONTRIBUTING.md` for the one rule (no threshold change without re-measuring both benchmarks) and the pull-request checklist. Copies of the threshold report (`THRESHOLDS_REPORT.md`) are regenerated via `python -m zspace_engine.thresholds_report`.

# Availability

Source: `https://github.com/Zierax/Axiom-Zspace` (tag `v1.1.1`), `zspace_engine/` (ingestion, detectors, validator, thresholds), `C99-Version/` (supplementary C99 port, \texttt{Purce} \cite{purce2024} at \url{https://github.com/Zierax/Purce}), `paper/` (JOSS `paper.md` + `paper.bib` archived with tag). License: MIT (`LICENSE`), archived on Zenodo \texttt{10.5281/zenodo.22255875}. Dependencies: Python $\geq$3.10 + `libc`/`libm` + OpenMP (optional for C99) + Python stack (`requirements.txt` pinned).

This paper describes the Axiom-ZSpace software itself — its architecture, reproducibility guarantees, and the C99 differential-verification artifact — rather than novel astrophysical results.

Example (offline, <30s):

```bash
pip install -r requirements.txt
python run_pipeline.py --synthetic                 # Python (default, reference)
cd C99-Version && make bin/zspace_card && cd .. && python run_pipeline.py --synthetic --engine c99  # C99 supplementary
python -m pytest tests/ -q                         # 101 tests
```

# Acknowledgements

We thank the Astropy, lightkurve, and batman communities for open-source tools that make this work possible. No external funding was received for this work.

# References
