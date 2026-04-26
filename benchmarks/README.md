# Axiom-ZSpace Benchmark System

This directory contains benchmark results for the Axiom-ZSpace exoplanet detection pipeline tested against confirmed planets from the NASA Exoplanet Archive.

## Purpose

The benchmark system validates the pipeline's performance by:
- Testing detection accuracy on known planets
- Measuring false negative rate (missed planets)
- Analyzing validation logic effectiveness
- Identifying areas for improvement

## Running a Benchmark

### Basic Usage

```bash
# Test against 500 known planets (default)
python benchmark_known_planets.py

# Test against 100 planets
python benchmark_known_planets.py --limit 100

# Custom output directory
python benchmark_known_planets.py --output my_benchmarks
```

### Requirements

```bash
pip install matplotlib seaborn astroquery astropy numpy
```

## Benchmark Output

Each benchmark run creates a timestamped directory with:

```
benchmarks/YYYYMMDD_HHMMSS/
├── BENCHMARK_REPORT.md          # Executive summary with key metrics
├── detailed_analysis.md         # Per-planet results and analysis
├── results.json                 # Raw benchmark data
├── charts/                      # Visualization charts
│   ├── detection_status.png     # Detection vs missed vs failed
│   ├── verdict_distribution.png # Sovereign verdict breakdown
│   ├── snr_distribution.png     # SNR histogram
│   ├── cvs_distribution.png     # CVS score histogram
│   ├── period_error.png         # Period recovery accuracy
│   └── processing_time.png      # Performance metrics
└── validation/                  # Individual validation cards
    └── Discovery_*.json
```

## Key Metrics

### Detection Rate
Percentage of known planets successfully detected by BLS and passing initial gates.

**Target:** >95% (Excellent), >85% (Good), >70% (Acceptable)

### False Negative Rate
Percentage of known planets missed or incorrectly rejected as false positives.

**Target:** <5% (Excellent), <15% (Acceptable), >15% (Needs Improvement)

### Sovereign Verdict Distribution
- **SOVEREIGN_PASS**: High confidence planet candidates
- **CONDITIONAL_PASS**: Moderate confidence, flagged for review
- **FALSE_POSITIVE**: Incorrectly rejected known planets (should be minimal)

## Interpreting Results

### High False Negative Rate
If many known planets are missed:
1. Check BLS detection thresholds (SNR, FAP)
2. Review lightcurve quality filters
3. Investigate common failure patterns

### False Positive Verdicts on Known Planets
If known planets are classified as FALSE_POSITIVE:
1. Review validation thresholds (density ratio, shape ratio)
2. Consider making strict tests non-critical
3. Add uncertainty propagation to physical tests
4. Check for systematic biases (e.g., short periods, grazing transits)

### Processing Failures
If many planets fail to process:
1. Improve error handling for edge cases
2. Add better lightcurve quality checks
3. Handle missing stellar parameters gracefully

## Benchmark History

Each benchmark run is preserved with a timestamp. Compare results across runs to track improvements:

```bash
# List all benchmark runs
ls -lt benchmarks/

# Compare two runs
diff benchmarks/20260418_013000/BENCHMARK_REPORT.md \
     benchmarks/20260418_020000/BENCHMARK_REPORT.md
```

## Best Practices

1. **Run benchmarks after major changes** to validation logic or thresholds
2. **Test with different planet populations** (hot Jupiters, Earth-sized, long periods)
3. **Document threshold changes** and their impact on metrics
4. **Track false positive verdicts** to avoid over-tuning
5. **Monitor processing time** to ensure scalability

## Troubleshooting

### "No planets fetched"
- Check internet connection
- Verify astroquery is installed: `pip install astroquery`
- NASA Exoplanet Archive may be temporarily unavailable

### "matplotlib not installed"
- Charts will be skipped but reports will still generate
- Install with: `pip install matplotlib seaborn`

### Slow processing
- Reduce `--limit` for faster testing
- Enable FITS caching in `config/production.yaml`
- Use `--max-targets` for quick validation

## Contributing

When modifying validation logic:
1. Run benchmark before changes (baseline)
2. Make your changes
3. Run benchmark after changes
4. Compare results and document impact
5. Ensure false negative rate doesn't increase significantly

---

**Last Updated:** 2026-04-18
