# Axiom-ZSpace Controlled-Honesty Evaluation

**Generated:** 2026-08-16 23:50:38 UTC
**Run:** PROBE_FPR68d
**Started:** 2026-08-14 23:01:45 UTC
**Method:** full local pipeline, BLIND search (no period prior, no ground-truth
expected period), validator archive query stubbed OFFLINE so the sovereign
verdict reflects ONLY intrinsic physics/QA gating.

## Key Metrics

| Metric | Value |
|--------|-------|
| True (injected) planets | 6 |
| Contamination targets | 8 |
| **Recall@correct period (TP)** | 0.0%  (0/6) |
| **Detection recall (any period)** | 0.0%  (0/6) |
| Detected-at-any-period (raw BLS) | 6/6 |
| **FPR (contamination certified as planet)** | 0.00%  (0 FP / 8 ) |
| Wrong-ephemeris certs on TRUE set | 0 (certified at a period outside the 5% window) |
| **Precision (certified & correct-period)** | 0.0% |
| **F1 (period-level)** | 0.000 |
| Confusion | TP=0 FP=0 FN=6 TN=8 |

## Recall by injected SNR

| SNR tier | recovered/n | recall |
|---|---|---|
| 5.5-8 | 0/2 | 0.0% |
| 8-14 | 0/2 | 0.0% |
| 14-30 | 0/2 | 0.0% |

## Per-subkind contamination results

| kind | count | certified-FP |
|---|---|---|
| eb | 1 | 0 |
| grazing_eb | 1 | 0 |
| noise | 2 | 0 |
| rotation | 2 | 0 |
| single_event | 2 | 0 |


## Files

- `results_true.json` — per-target results for injected planets
- `results_false.json` — per-target results for contamination
- `validation/` — sovereign validation cards written by the pipeline
- `EVALUATION_REPORT.md` — this report

**Note on honesty:** raw "detected" flags exclude nothing; wrong-ephemeris
certifications on the true set are counted separately from contamination
false-certifications (the FPR row above is contamination-only).
