# Axiom-ZSpace Pipeline Benchmark Report

**Generated:** 2026-04-26 19:57:44 UTC  
**Benchmark ID:** 20260426_152109  
**Test Dataset:** NASA Exoplanet Archive Confirmed Planets (TESS)

---

## Executive Summary

This benchmark tests the Axiom-ZSpace exoplanet detection pipeline against 619 confirmed planets from the NASA Exoplanet Archive to measure detection accuracy, false negative rate, and validation performance.

### Key Metrics

| Metric | Value | Percentage |
|--------|-------|------------|
| **Total Planets Tested** | 619 | 100.0% |
| **Successfully Detected** | 574 | 92.7% |
| **Missed (False Negatives)** | 3 | 0.5% |
| **Processing Failures** | 42 | 6.8% |

### Detection Performance

| Category | Count | Percentage |
|----------|-------|------------|
| **True Positives** | 574 | 92.7% |
| **False Negatives** | 3 | 0.5% |
| **SOVEREIGN_PASS** | 61 | 9.9% |
| **CONDITIONAL_PASS** | 6 | 1.0% |
| **FALSE_POSITIVE Verdict** | 3 | 0.5% |

### Average Performance Metrics

| Metric | Value |
|--------|-------|
| **Average SNR** | 60.2σ |
| **Average CVS Score** | 0.496 |
| **Average Period Error** | 71.43% |
| **Average Processing Time** | 26.8 seconds |

---

## Analysis

### Detection Rate: 92.7%

**GOOD** - The pipeline has strong detection performance with room for improvement.

### False Negative Rate: 0.5%

**EXCELLENT** - Very low false negative rate, minimal risk of missing real planets.

---

## Validation Verdict Distribution

The sovereign verdict determines whether a detected signal is classified as a planet candidate or false positive:

- **SOVEREIGN_PASS**: 61 (9.9%) - High confidence planet candidates
- **CONDITIONAL_PASS**: 6 (1.0%) - Moderate confidence, flagged for review
- **FALSE_POSITIVE**: 3 (0.5%) - Known planets incorrectly rejected


### ⚠️ False Positive Verdicts on Known Planets

3 known planets were incorrectly classified as FALSE_POSITIVE. This indicates the validation thresholds may be too strict. Review the detailed results to identify which tests are causing false rejections.

---

## Recommendations

### 3. Improve Robustness
- Add better error handling for edge cases
- Improve lightcurve quality checks
- Handle missing stellar parameters gracefully

---

## Files Generated

- `results.json` - Raw benchmark results
- `BENCHMARK_REPORT.md` - This summary report
- `detailed_analysis.md` - Detailed per-planet results
- `charts/` - Visualization charts and diagrams

---

**End of Report**
