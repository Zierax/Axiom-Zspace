# Axiom-ZSpace Controlled-Honesty Evaluation

**Generated:** 2026-08-16 23:50:38 UTC
**Run:** OVERHARM_FIX2
**Started:** 2026-08-16 02:47:37 UTC
**Method:** full local pipeline, BLIND search (no period prior, no ground-truth
expected period), validator archive query stubbed OFFLINE so the sovereign
verdict reflects ONLY intrinsic physics/QA gating.

## Key Metrics

| Metric | Value |
|--------|-------|
| True (injected) planets | 100 |
| Contamination targets | 80 |
| **Recall@correct period (TP)** | 52.0%  (52/100) |
| **Detection recall (any period)** | 53.0%  (53/100) |
| Detected-at-any-period (raw BLS) | 100/100 |
| **FPR (contamination certified as planet)** | 0.00%  (0 FP / 80 ) |
| Wrong-ephemeris certs on TRUE set | 1 (certified at a period outside the 5% window) |
| **Precision (certified & correct-period)** | 98.1% |
| **F1 (period-level)** | 0.680 |
| Confusion | TP=52 FP=1 FN=47 TN=80 |

## Recall by injected SNR

| SNR tier | recovered/n | recall |
|---|---|---|
| 5.5-8 | 2/34 | 5.9% |
| 8-14 | 19/33 | 57.6% |
| 14-30 | 31/33 | 93.9% |

## Per-subkind contamination results

| kind | count | certified-FP |
|---|---|---|
| eb | 8 | 0 |
| grazing_eb | 8 | 0 |
| noise | 32 | 0 |
| rotation | 16 | 0 |
| single_event | 16 | 0 |


## Files

- `results_true.json` — per-target results for injected planets
- `results_false.json` — per-target results for contamination
- `validation/` — sovereign validation cards written by the pipeline
- `EVALUATION_REPORT.md` — this report

**Note on honesty:** raw "detected" flags exclude nothing; wrong-ephemeris
certifications on the true set are counted separately from contamination
false-certifications (the FPR row above is contamination-only).
