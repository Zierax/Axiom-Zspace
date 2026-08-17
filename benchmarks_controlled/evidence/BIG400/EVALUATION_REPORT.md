# Axiom-ZSpace Aggregated Controlled Evaluation

**Generated:** 2026-08-17 02:06:38 UTC
**Suite:** BIG400 (aggregated from 8 chunk runs)
**Composition:** 400 true + 400 false targets; chunks of 50+50 with seeds 20260816..20260823 (seed = suite_seed + chunk index).
**Method:** full local pipeline, BLIND search (no period prior, no ground-truth
expected period), validator archive query stubbed OFFLINE so the sovereign
verdict reflects ONLY intrinsic physics/QA gating. Metrics computed with the
canonical compute_metrics() of benchmarks_controlled/run_controlled.py.

## Key Metrics

| Metric | Value |
|--------|-------|
| True (injected) planets | 400 |
| Contamination targets | 400 |
| **Recall@correct period (TP)** | 41.2%  (165/400) |
| **Detection recall (any period)** | 46.2%  (185/400) |
| Detected-at-any-period (raw BLS) | 400/400 |
| **FPR (contamination certified as planet)** | 4.25%  (17 FP / 400 ) |
| Wrong-ephemeris certs on TRUE set | 20 (certified at a period outside the 5% window) |
| **Precision (certified & correct-period)** | 81.7% |
| **F1 (period-level)** | 0.548 |
| Confusion | TP=165 FP=37 FN=215 TN=383 |

## Recall by injected SNR

| SNR tier | recovered/n | recall |
|---|---|---|
| 5.5-8 | 3/136 | 2.2% |
| 8-14 | 56/136 | 41.2% |
| 14-30 | 106/128 | 82.8% |

## Per-subkind contamination results

| kind | count | certified-FP |
|---|---|---|
| eb | 40 | 8 |
| grazing_eb | 40 | 6 |
| noise | 160 | 1 |
| rotation | 80 | 2 |
| single_event | 80 | 0 |

## Provenance (chunks)

| chunk | true | false | seed |
|---|---|---|---|
| c0 | 50 | 50 | 20260816 |
| c1 | 50 | 50 | 20260817 |
| c2 | 50 | 50 | 20260818 |
| c3 | 50 | 50 | 20260819 |
| c4 | 50 | 50 | 20260820 |
| c5 | 50 | 50 | 20260821 |
| c6 | 50 | 50 | 20260822 |
| c7 | 50 | 50 | 20260823 |

Each chunk is independently reproducible from its seed with run_controlled.py.

## Reproduce this suite

```bash
python benchmarks_controlled/run_controlled.py --true 50 --false 50 \
    --out benchmarks_controlled/runs/<cN> --seed <suite_seed + N>    # every chunk N in 0..7
python scripts/aggregate_benchmark_runs.py --chunks-dir benchmarks_controlled/runs/<dir> \
    --out benchmarks_controlled\evidence\BIG400 --suite-name BIG400 --suite-seed 20260816
```

## Files

- `results_true.json` — merged per-target results for injected planets
- `results_false.json` — merged per-target results for contamination
- `chunks.json` — chunk/seed provenance
- `EVALUATION_REPORT.md` — this report

**Note on honesty:** wrong-ephemeris certifications on the true set are counted
separately from contamination false-certifications (the FPR row is contamination-only).
