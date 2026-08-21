# Testing

The verification suite is the contract of this repository: **101 tests,
offline-green, covering determinism, calibration, physics audits, and the
circuit-breaker invariant**. This document maps every test file to what it
proves, explains how determinism is enforced, and how to add tests.

## 1. Running

```bash
python -m pytest tests/ -q                 # everything, ~5 min
python -m pytest tests/test_reproducibility.py -q   # determinism only
python -m pytest tests/test_synthetic.py -q        # synthetic contract
python -m pytest -k "fp4 or fp7" -q               # subset by keyword
```

Run from the repository root. Network-dependent tests (catalog cross-checks)
mock or degrade gracefully, so the suite is green offline **and** online —
this is a design requirement, tested explicitly.

## 2. Suite map

| Test file | Proves / covers |
|---|---|
| `test_reproducibility.py` | **Determinism contract** (the 7 tests): same seed ⇒ same sample ⇒ same results for the controlled benchmark; asserted against the fixed seed used for evidence runs |
| `test_synthetic.py` | Synthetic-target contract: generated targets carry every field the pipeline needs (`cvs_verdict`, periods, transit parameters); injector shapes are stable |
| `test_period_prior_selector.py` | Period-prior selector: explicit prior-power floor + harmonic ladder ratio (asserts the (3000/5000)² ≈ 0.35 floor logic) — the "prior, not lock" guarantee |
| `test_physics_audit_fixes.py` | Physics audits under noisy fixtures (batman model + noise); a true transit still certifies end-to-end through the full blind path |
| `test_ephemeris_identity_gate.py` | **Wrong-ephemeris gate**: candidates whose recovered period contradicts archive truth are demoted, never certified as SOVEREIGN_PASS |
| `test_catalog_integration.py` | Threshold-catalog integration into the engine (config values reach detectors/validator) |
| `test_requirements_verification.py` | Runtime requirements & metadata verification (package importability); source of the benign SIMBAD warnings |
| `test_check_external_catalogs.py` | `check_external_catalogs()` error paths: invalid TIC IDs, timeouts, offline behavior — must never break a verdict |
| `test_task_2_2_verification.py` | Per-target verification records (candidate → result → verdict provenance) |
| `test_task_3_2_demo.py`, `test_task_3_2_fp9_integration.py` | End-to-end demo path + integration of the FP-rejection chain (strict-gate demotion) |
| `test_task_4_1_fp4_shape_ratio.py` | **FP-4 calibration**: U/V shape ratio ≥ 0.4 separates true transits from V-shaped artifacts on synthetic distributions |
| `test_task_4_2_combined_shape_density.py` | **FP-7 + shape combination**: density-band [0.2, 5.0] and shape jointly reject false classes that pass either gate alone |
| `test_task_5_1_circuit_breaker.py` | **Circuit breaker**: a failed critical gate makes SOVEREIGN_PASS impossible, regardless of other scores |
| `test_task_7_1_conflict_detection.py`, `test_task_7_1_integration.py` | **Conflict detection**: high-SNR/contradicting-gate conflicts are recorded as explicit proof-chain entries |
| `test_task_7_2_conflict_logging.py` | Conflicts are surfaced in the per-target log/report (closure of §7.1) |
| `test_task_8_1_tic_id_passing.py` | TIC ID passthrough end-to-end (identifier survives detection → validation → card) |

Helpers (not collected by pytest, run explicitly):

| Helper | Use |
|---|---|
| `tests/run_validator.py` | Manual validator demo: `--synthetic`, `--tic <id>`, `--force-known` — inspect a single verdict + proof chain |
| `tests/demo_catalog_check.py` | Interactive catalog sanity demo |

## 3. How determinism is enforced

`tests/test_reproducibility.py` runs the controlled benchmark **twice** with
the same seed and asserts bit-identical results, plus asserts the fixed-seed
snapshot numbers. The contract lives in `benchmarks_controlled/runner`
(seed-routed RNGs for injector and engine) — any code change that makes two
same-seed runs diverge fails the suite *before* any human reads benchmark
numbers.

## 4. Adding tests (the rules)

1. **No test may require network to pass.** If you need external data, mock
   it (pattern: `test_check_external_catalogs.py`) or write an offline branch.
2. **Assert measured truth, not wishes.** Thresholds in tests must match the
   **calibrated** profile (e.g. `fp4_shape_min == 0.4`, `fp7_density_band ==
   [0.2, 5.0]`). Calibration expectations are part of the regression contract.
3. **Name the contract.** Prefer explicit names
   (`test_circuit_breaker_no_critical_fail_pass`) over `test_thing_works`.
4. **Keep runtime sane.** Target < 60 s per file where possible; the full
   suite should stay ~5 min so it is cheap to run before every push.
5. **Update `docs/TESTING.md`** (this map) when you add a file.

## 5. The C99 engine verification (separate from pytest)

The C99 sovereign engine (`C99-Version/`, see
[`docs/C99_ENGINE.md`](C99_ENGINE.md)) is guarded by its own differential
harnesses — they are not part of the 101-test pytest suite because they
require a C toolchain (WSL on Windows):

```bash
cd C99-Version
make bin/verify_kernels && python tests/verify_compare.py   # 148/148 kernels vs Python
python tests/parity_card.py                                  # 90/90 full-card parity
```

- `tests/gen_verify_kernels.py` + `verify_compare.py` run **every generated
  kernel** on randomized inputs against its Python original (rel. tol 1e-9).
- `tests/parity_card.py` runs the actual `zspace_card` binary against the
  Python `ProofEngine` on synthetic candidates and light curves, comparing
  every numeric field and per-test FP verdict (rel. tol 2e-3, abs 0.02).

A rebuild of the kernels (Purce re-run) must re-pass both before any
`--engine c99` result is quoted.

## 6. Before pushing code

```bash
python -m pytest tests/ -q          # 101 passed, 0 failed
python -m zspace_engine.thresholds_report  # if thresholds touched
```

Anything that changes detection/validation *behavior* (not just docs/tests)
requires re-running **both** benchmarks first — see `docs/CONTRIBUTING.md`.