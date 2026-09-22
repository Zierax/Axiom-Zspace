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
date: 5 September 2026
bibliography: paper.bib
---

# Summary

Axiom-ZSpace searches TESS and Kepler light curves for transits without being told where to look. It chains ingestion, BLS detection over a period-prior ladder, ephemeris resolution, physics audits, an 11-gate ruling engine with a circuit breaker, and CVS scoring into one deterministic pipeline. Every check writes its reasoning into a proof chain on the output card, and every tunable number lives in a single catalog (`config/production.yaml`, read through `thresholds.py`).

Two benchmarks pin the numbers. The synthetic BIG400 set (400 injected planets, 400 contaminants, 8 seeds) gives 41.2% recall at the right period, 81.7% precision, and 4.25% contamination FPR. Twelve Kepler hosts with confirmed planets plus twelve quiet stars give 41.7% recall. A C99 port of the same logic ships alongside for batch runs.

# Statement of Need

Most transit codes mix detection, vetting, and thresholds in one place, so a silent bug or a hand-tuned cut is hard to catch. Axiom-ZSpace separates them: detection proposes, eleven gates dispose, and the catalog records every threshold with the measurement behind it. Change a number and both benchmarks must be re-run. That rule lives in `CONTRIBUTING.md` and the test suite enforces it.

Nothing here needs a hint or a network connection. The same seed gives the same result (covered by `tests/test_reproducibility.py`), and the synthetic suite plus all 101 tests run offline. This is for people running blind searches who want to inspect why a candidate passed: archive teams, follow-up groups vetting candidates, and students who need a baseline they can reproduce.

# State of the Field

* **BLS** [@kovacs2002] is the standard box search; **TLS** [@hippke2019] adds limb-darkened templates at $\sim$1$\times$ BLS cost.
* **GPU BLS**: `cuvarbase` [@cuvarbase] and **GTLS** [@hu2026] report 10--100$\times$ on GPU; **QLP GPU** [@kunimoto2023] reports 40$\times$.
* **Approximators**: **fBLS** [@shahaf2022] 15$\times$ at 65k points (binned); **GPFC** [@wang2024] reports 15$\times$ vs `astropy` cython at equal grid (CNN, USP only). None ship the full chain (gates, catalog, proof chains) with versioned synthetic and real benchmarks next to the code.

# Installation

**Requirements:** Python ≥3.10, gcc ≥10, make, cmake ≥3.20, OpenMP runtime (`libgomp`/`libomp`), WSL2 Ubuntu 22.04 on Windows. Pin `requirements.txt` (`pip install -r requirements.txt`) — optional `lightkurve`, `astroquery`, `batman-package` for network features only; synthetic benchmark is offline.

```bash
pip install -r requirements.txt
cd C99-Version && make bin/zspace_card bin/verify_kernels  # builds bin/zspace_card and bin/verify_kernels
cd .. && python run_pipeline.py --synthetic                 # offline self-test
```

C build: `gcc -O3 -march=native -mtune=native -flto -ffast-math -fopenmp -fopenmp-simd` (native) or `-O2 -fno-fast-math` for verified IEEE build (see `C99-Version/Makefile`; `C99_ENGINE.md`).

# Software design

Six stages, each with explicit inputs and outputs:

* **Ingestion** (`zspace_engine/ingestion.py`): MAST fetch via `lightkurve` with disk cache, `quality==0` masking, sigma clipping, median normalization, and Savitzky-Golay detrending (single `flat1`, window 3.0 d or 0.75P, x in [-1,1] QR via `zspace_ingestion.c`).
* **Detection** (`detectors.py`): BLS periodogram (`astropy.timeseries.BoxLeastSquares` baseline) over frequency grid $n_{\mathrm{freq}}=\max(\lfloor(f_{\max}-f_{\min})/df\rfloor,2000)$, $df=1/(T \cdot 20)$, duration 0.25--8 h, with a **ladder of $k20$ strict local maxima** filtered by $\tau/P>0.15$, $|\log|<0.10$, $min\_rel=0.05$, and a self-calibrating exponential-tail FAP (MAD).
* **Ephemeris** (`ephemeris.py`): fold, merge dip signatures, resolve alias/harmonic ambiguity, refine $t_0$/duration.
* **Audits** (`auditors.py`): five physics audits — even/odd Welch ($\Delta<3.0$), depth consistency, secondary eclipse, ingress/egress, and limb-shape — each returning a dataclass with evidence.
* **Validation** (`validator.py`): 11-gate ruling engine with circuit breaker. Critical gates ($S/N \ge 5.5$, $FAP \le 0.05$) make `SOVEREIGN_PASS` impossible when failed; non-critical allowances are `verdict_max_fail_pass=2`, `conditional=3`.
* **Classification** (`core.py`): Composite Vitality Score (CVS) $w=(0.97,0.83,0.61,0.31)$ over $S_{\mathrm{periodicity}}, S_{\mathrm{depth}}, S_{\mathrm{limb}}, S_{\mathrm{stellar}}$ with four tiers ($\ge0.80$ PLANET, $\ge0.55$ LIKELY, $\ge0.35$ AMBIGUOUS).

All tunable constants (gates, FAP, ladder, CVS weights) live in one file (`thresholds.py` to `config/production.yaml`) with per-key evidence in `THRESHOLDS_REPORT.md`. The 101-test suite checks determinism, the circuit breaker, and ephemeris identity as runnable contracts. A threshold change without re-measuring BIG400 and REAL_FINAL is rejected.

# Functionality

**Dual-engine parity.** Python default for single-target reference; `--engine c99` recommended for batch:
```bash
python run_pipeline.py --synthetic --engine c99          # batch C99 (requires prior make)
python benchmarks_controlled/run_controlled.py --true 50 --false 50 --seed 20260816 --engine c99
```

**Pipeline (production `balanced`, `frequency_factor 20`, `k20`, `coherent OFF`).** Normalize, then Savitzky-Golay `flat1` (3.0 d window). `BLS.search` over a frequency-duration grid, then `top_candidates` keeps strict local maxima. Each candidate goes through fold, audits, `eph_resolve`, and the 11-gate `sovereign_validate`. The first `SOVEREIGN_PASS` or `CONDITIONAL_PASS` certifies `OFFLINE_NEW_DISCOVERY`; otherwise the first status stands. Output is a JSON card with the full chain.

**CLI contract (`bin/zspace_card`):** reads `key=value` candidate + optional CSV `time,flux` and prints a single JSON sovereign card; `c99_bridge.py` auto-detects `bin/` vs `build/` and `D:/`→`/mnt/d/` for WSL. `FP-10` (`count_observed_transits`) is computed from the time/flux series, not hard-coded.

# Threshold Catalog and Provenance

The catalog defines three profiles: `conservative`, `balanced` (default), and a looser experimental `sensitive`. Every key records its direction (e.g., $S/N\ge5.5$, shape $\ge0.4$, density $[0.2,5.0]$, impact $<0.9$, $N_{\mathrm{tr}}\ge2$), its weight (`critical`/`major`), and a short evidence note on what tightening or loosening costs. `THRESHOLDS_REPORT.md` regenerates via `python -m zspace_engine.thresholds_report` and ships with every catalog change. Every discovery card is JSON with its proof chain, and benchmark evidence sits versioned under `benchmarks_controlled/evidence/BIG400` and `benchmarks_real/evidence/REAL_FINAL`. Run outputs and caches stay git-ignored.

# Performance — CPUs as Efficient Alternatives to GPU Acceleration

The headline numbers belong to the catalog: 41.2% recall and 4.25% FPR on BIG400, 41.7% recall on REAL_FINAL. The C99 port exists for batch work. It does 42.8 ms per target on 3k-point curves (650$\times$ vs single-thread Python at 27.8 s, 40.6$\times$ per core) and 4.8 s per target on 87k-point 5-sector curves (16 cores, `-O3 -march=native -flto -ffast-math -fopenmp-simd`). The 87k Python baseline was not re-measured, so treat that ratio as provisional.

For reference, GTLS reports 15.7$\times$ on an RTX 4090 [@hu2026] and QLP GPU reports 40$\times$ [@kunimoto2023], both on dedicated hardware. The C99 port is bit-identical to Python and needs no GPU. Where the cost is gate logic and FAP calibration rather than raw FLOPs, a plain CPU build can win.

# Verification

101 tests ship with the Python pipeline and run offline (`python -m pytest tests/ -q`): determinism, circuit breaker, ephemeris identity, gate calibration. The C99 port, derived mechanically via `Purce` [@purce2024], is checked against Python where the harness covers it: 148/148 validator kernels at $10^{-9}$ and 90/90 synthetic cards at $2\times10^{-3}$. The 30 BLS kernels compile and get validated through batch runs instead. Versioned evidence is BIG400 at 400/400 in Python; C99 parity on that sample is 90/90 synthetic.

# Research impact statement

The authors use this pipeline for TESS/Kepler archival searches and threshold calibration. What backs that up: versioned BIG400 and REAL_FINAL benchmarks anyone can re-run, 101 offline tests, and the C99 port (148/148 at $10^{-9}$, 90/90 cards) that makes batch scans practical at 650× on ordinary CPUs. Groups that need to show their work on vetting, rather than trust a cut, get a catalog and proof chains they can audit.

# AI usage disclosure

No generative AI tools were used to create the Axiom-ZSpace software, its tests, or its benchmarks. Limited AI assistance was used for language polishing of documentation and paper text; all AI-generated suggestions were reviewed, edited, and validated by the human author, who remains responsible for the correctness of the code, tests, and scientific claims.

# Community Guidelines

Contributions, bug reports, and support requests are welcome via the GitHub issue tracker at `https://github.com/Zierax/Axiom-Zspace/issues`. See `CONTRIBUTING.md` for the one rule (no threshold change without re-measuring both benchmarks) and the pull-request checklist. Copies of the threshold report (`THRESHOLDS_REPORT.md`) are regenerated via `python -m zspace_engine.thresholds_report`.

# Availability

Source: `https://github.com/Zierax/Axiom-Zspace` (tag `v1.1.2`), `zspace_engine/` (ingestion, detectors, validator, thresholds), `C99-Version/` (supplementary C99 port, `Purce` [@purce2024] at <https://github.com/Zierax/Purce>), `paper/` (JOSS `paper.md` + `paper.bib` archived with tag). License: MIT (`LICENSE`), archived on Zenodo `10.5281/zenodo.22255875`. Dependencies: Python $>=3.10$ + `libc`/`libm` + OpenMP (optional for C99) + Python stack (`requirements.txt` pinned).

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
