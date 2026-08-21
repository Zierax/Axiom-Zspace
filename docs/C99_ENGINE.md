# The C99 Engine

Axiom-ZSpace ships **two sovereign engines** with identical physics and
identical verdicts:

| Engine | Code | Purpose |
|---|---|---|
| `python` | `zspace_engine/validator.py` (reference) | default; the measured benchmark baseline |
| `c99` | `C99-Version/` (see below) | a self-contained C99 re-implementation, dependency-free |

The C99 engine is **not a wrapper** around the Python code: its math kernels
are machine-generated from a strict numpy subset by
[**Purce**](https://github.com/Zierax/Purce), compiled with a stock C99
compiler, and **differentially verified** against the Python reference before
any use. It links **no external libraries** — only `libc` (`<stdio.h>`,
`<stdlib.h>`, `<string.h>`, `<math.h>`).

## Quick reference

```bash
# run the pipeline with the C99 sovereign engine
python run_pipeline.py --synthetic --engine c99
python run_pipeline.py --tic 260128333 --engine c99

# controlled benchmark with the C99 engine (same sample as python)
python benchmarks_controlled/run_controlled.py --true 50 --false 50 \
    --seed 20260816 --engine c99 --out benchmarks_controlled/runs/MY_C99

# build + verify the C99 engine (requires WSL with gcc on Windows)
cd C99-Version
make bin/zspace_card                       # build the CLI binary
make bin/verify_kernels                    # build the differential verifier
python tests/verify_compare.py             # 148/148 kernels vs Python
python tests/parity_card.py                # 90/90 full-card parity
```

## Repository layout (`C99-Version/`)

```
C99-Version/
├── Makefile                 # bin/verify_kernels, bin/zspace_card (gcc -O3 -march=native -flto -fopenmp -std=c99, 7 TUs)
├── CMakeLists.txt           # Release -O3 -march=native -fopenmp (build/zspace_card)
├── purce_src/
│   └── zspace_kernels.py    # the strict numpy subset that Purce translates
├── generated/               # Purce output (committed): 148 kernels + headers
│   ├── purce_src.gen.h      # prototypes regenerated from kernel bodies (the reliable header)
│   └── purce_src_*.c        # one C file per kernel + *.prov.json provenance
├── src/
│   ├── zspace_core.c/.h     # §1–§6 sovereign logic, FP gate engine, CVS
│   ├── zspace_bls.c/.h      # BLS periodogram, ladder top_candidates (k20), FAP
│   ├── zspace_eph.c/.h      # ephemeris/fold, density
│   ├── zspace_audit.c/.h    # even/odd, depth_cons, secondary, ingress
│   ├── zspace_ingestion.c/.h# SG flatten (QR, x∈[-1,1], rdiag)
│   ├── zspace_batch.c       # batch: single flat1, BLS 20x20, ladder validate
│   └── zspace_card.c        # CLI: batch/manifest, bls, card; key=value + CSV
├── tests/
│   ├── gen_verify_kernels.c # generates verify_kernels.c from the kernel list
│   ├── verify_compare.py    # differential kernel verification (148/148 PASS)
│   └── parity_card.py       # full-card parity (90/90 PASS, tol 2e-3)
└── bin/ | build/            # built binaries (git-ignored, c99_bridge auto-detects both)
```

## How the engines stay identical

1. **One math source.** Every numeric formula of the sovereign validator is
   written once as a strict numpy kernel in `purce_src/zspace_kernels.py`.
   Purce translates each kernel to a standalone C function. The rule that
   keeps translation honest: **one line = one numpy call** — no arithmetic
   inside numpy-call arguments, no nested calls. (Purce v0.1.0 does not
   reliably translate composite expressions or `BinOp` nodes inside call
   arguments; known broken forms are documented in the kernel file header.)
2. **Differential verification.** `make bin/verify_kernels` produces a C
   program that runs every kernel on randomized inputs; `verify_compare.py`
   runs the same inputs through the Python originals and compares (relative
   tolerance 1e-9). **148/148 kernels PASS** on every rebuild.
3. **Full-card parity.** `tests/parity_card.py` runs the actual `zspace_card`
   binary against the Python `ProofEngine` on synthetic candidates and light
   curves (with and without a light curve), comparing every numeric field
   (relative tol 2e-3) plus verdicts and per-test FP gate results:
   **90/90 PASS**.
4. **Identical gate logic.** `src/zspace_core.c` reimplements the §1–§6
    decision chain (Kepler consistency, geometric consistency, density
    constraint, transit probability, the 11-test false-positive ruling engine,
    the chi² fit and the sovereign verdict assembly) with the same thresholds
    (`zspace_engine/thresholds.py`), read as a fixed catalog at compile time.
    FP-10 (`count_observed_transits`) is computed from the raw time/flux
    series in C, not hardcoded.
    > **Engine divergence (documented):** `frequency_factor 20.0` catalog (`config/production.yaml:20`, `zspace_batch.c:176`) vs `zspace_engine/detectors.py:223` default `10.0` (sector path `sector_processor.py:635` still `10.0`); `ladder k20 min_rel 0.05` hard-coded in C `zspace_batch.c:184`/`zspace_bls.c:584` (Python `thresholds`); `n_freqs` `int()` vs `floor()` `zspace_bls.c:83` byte-identical for positive range; `single flat1` reuse `zspace_batch.c:153` vs `zspace_card.c:688` double flatten; `coherent 0` hard-coded in C `zspace_batch.c:271` (Python profile-dependent, `balanced OFF`).
5. **Controlled-benchmark parity.** Running the same 50+50 sample
    (seed 20260816) through both engines produces identical recall/FPR/
    precision/F1 — see `docs/BENCHMARKS.md` §1.7.

## The CLI contract (`zspace_card`)

The binary reads a candidate file (`key=value` lines) and an optional CSV
light curve (header `time,flux[,flux_err,model_flux]` — two columns suffice;
the chi² section is only computed when all four columns are present), and
prints a single JSON card to stdout:

```json
{
  "schema": "Axiom-ZSpace Sovereign ...",
  "sovereign_verdict": "SOVEREIGN_PASS | CONDITIONAL_PASS | FALSE_POSITIVE",
  "cvs": 0.89242684,
  "cvs_verdict": "PLANET",
  "all_sections_pass": true,
  "n_transits": 2,
  "section_1_kepler": { "a_au": ..., "residual_si_pct": ..., ... },
  "section_2_geometry": { ... },
  "section_3_density": { "density_ratio": ..., ... },
  "section_4_probability": { "P_tr": ..., "impact_parameter_b": ..., ... },
  "section_5_fp_ruling": { "n_pass": 11, "n_tests": 12,
                           "fp_verdicts": [...], "overall_verdict": ... },
  "inputs": { ... }
}
```

`c99_bridge.py` (repo root) is the thin Python adapter used by
`run_pipeline.py` and the controlled benchmark: it writes the candidate file
+ CSV to a temp dir, invokes the binary (directly when an `.exe` exists,
otherwise through WSL with automatic `D:/` → `/mnt/d/` path conversion), and
parses the card. On malformed JSON it raises with the raw stdout dumped to
`%TEMP%\zspace_c99_stdout.txt` for diagnosis.

## Build & platform notes

- **Linux/macOS (Makefile, release):**
  ```make
  CC=gcc
  CFLAGS=-std=c99 -O3 -march=native -mtune=native -flto -ffast-math -fno-math-errno -funroll-loops -DNDEBUG -Wall -Wextra -fopenmp -fopenmp-simd -Igenerated -Isrc
  bin/zspace_card: src/zspace_card.c src/zspace_core.c src/zspace_bls.c src/zspace_eph.c src/zspace_audit.c src/zspace_ingestion.c src/zspace_batch.c $(GEN)
  ```
  Output `C99-Version/bin/zspace_card` (Makefile) and `C99-Version/build/zspace_card` (CMake `Release -O3 -march=native -fopenmp`). `c99_bridge.py` auto-detects both and converts `D:/`→`/mnt/d/` for WSL.
- **Windows:** the produced binary is an ELF — build and run it under WSL
  (`wsl bash -lc "cd /mnt/d/.../C99-Version && make bin/zspace_card"`). The
  bridge and the parity harness handle the path conversion automatically.
- **Re-generating kernels** (after editing `zspace_kernels.py`):
  `purce compile <SRC> -o C99-Version/generated` — then **regenerate the
  header** from the kernel bodies (regex `^void (\w+)\(([^)]*)\)\s*\{` over
  the `.c` files; the Purce-generated header is not reliable), then rebuild
  and re-run both verifications.

## Why a second engine at all?

Independent implementation, one source of truth. Two engines that agree on
90/90 cards and 148/148 kernels make single-language bugs visible: a
disagreement between the engines is a bug in at least one of them, and the
C99 binary provides a path for deployment targets where a Python interpreter
is not wanted.

---

*Engine status: differential verification 148/148 PASS, card parity 90/90
PASS, controlled-benchmark agreement measured — `docs/BENCHMARKS.md` §1.7.*
