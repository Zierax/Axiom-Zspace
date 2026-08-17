# Production use

How to operate Axiom-ZSpace at scale, what is measured, and what is
**not** — the difference matters. The measured claims of v1.0 are the
single-target engine behavior (benchmarks). Batch/sector operation reuses the
same engine but has its own operational envelope that v1.0 did not
re-measure.

## 1. What is measured (v1.0)

| Scope | Status |
|---|---|
| Single-target blind pipeline (TIC-level: detection → audits → validation → card) | **Measured** — controlled BIG400 (400+400, 8 seeds) + fixed-seed 100+80 + real 12+12 benches |
| Certified-candidate FPR on the controlled false set | **Measured** — 4.25% across 800 targets (0/80 on the fixed-seed suite), balanced profile, override OFF |
| Threshold catalog | **Measured** — values carry evidence strings in THRESHOLDS_REPORT.md |

## 2. What is NOT measured (operational caveats — read before trusting batch output)

| Item | Status |
|---|---|
| Whole-sector batch certification rates | **Not re-measured in v1.0.** Batch mode exists (`--sector`, `--sectors`); its per-sector statistics were not part of the v1 evidence. Do not quote batch recall/FPR from any source except fresh runs. |
| v2.x-era sector aggregates (discoveries.json per sector, historical counts) | **Not preserved** — CHANGELOG documents them as legacy history only |
| `coherent_override_enabled: true` (sensitive profile) | **EXPERIMENTAL** — measured FPR explosion (12.5–62.5% on probes) when ON |
| Optional filters `--mcmc`, `--tpf-centroids`, `--multi-sector` | Off by default; their effect on FPR/recall is not separately benchmarked |

## 3. Operating the batch path

```bash
python run_pipeline.py --sector 42                    # full sector
python run_pipeline.py --sectors 1-50 --max-targets 200   # capped trial
```

- Dispatches through `zspace_engine/sector_processor.py`, one engine instance
  per target; per-sector `summary.json` + `discoveries.json` land under
  `output.base_directory` (default `axiom_output/`), routed per-verdict
  (`rejected/` for FALSE_POSITIVE) by `output_organizer.py`.
- **Cost**: hours per sector, GB-scale FITS cache (`fits_cache.cache_dir`,
  size-capped by `max_cache_size_gb`). Run on a machine with stable storage
  and network; resume-friendly because cached light curves are reused.
- **Determinism**: single-target runs with the same inputs are deterministic;
  batch runs are reproducible run-to-run given the same cache/catalog state.

## 4. Configuration knobs (non-threshold)

| Setting | Default | Meaning |
|---|---|---|
| `network.api_timeout_seconds` / `max_retry_attempts` | 30 / 2 | MAST/API resilience |
| `fits_cache.enabled` / `cache_dir` / `max_cache_size_gb` | true / `.cache/fits` / 50 | Light-curve caching (limit enforcement deferred — monitor disk) |
| `output.base_directory` / `preserve_timestamps` | `axiom_output` / true | Output routing |
| `logging.level` / `file` / `console` | INFO / `axiom_pipeline.log` / true | Central logging |
| `fp_filters.run_mcmc` | false | MCMC off for speed |
| env `AXIOM_CONFIG` | — | Point the pipeline at a different config file |

## 5. Known error classes to watch (from the measured benches)

1. **Period alias/harmonics** — the leading certification-error class on real
   data (wrong-ephemeris certs: Kepler-37 d at 39.79 d etc.). The
   `ladder`/FP-5c machinery reduces but does not eliminate it; human sanity
   check of certified periods against the search band is recommended.
2. **Quiet-star certifications** — a ~33% proxy-FPR was measured on 12 quiet
   stars; treat unconfirmed candidates as candidates, not planets, until
   follow-up.
3. **Shallow-signal misses** — recall is 41.2% controlled (BIG400; 52/100 on the fixed-seed suite), 41.7% real; the
   pipeline is calibrated for *precision-first* science at SNR ≥ ~5.5, not
   completeness.
4. **Override temptation** — enabling `coherent_override_enabled` "to catch
   more" is a measured FPR explosion. If you must, the only defensible path
   is: re-measure both benchmarks with the new profile and publish the
   numbers.

## 6. Operational checklist before a large scan

```
□ config/production.yaml reviewed (profile explicit, override OFF)
□ disk free ≥ 2× expected cache size
□ python -m pytest tests/ -q                 (101 passing)
□ snapshot refreshed? python scripts/fetch_nea_snapshot.py --out benchmarks_real/data
□ trial run: --sector N --max-targets 5      (sanity: cards look right)
□ log level + output dir set; discovery cards git-ignored (they are)
□ follow-up channel for candidates decided before the run
```

## 7. Monitoring during a run

- Watch `axiom_pipeline.log` for per-target completion and gate conflicts
  (§6.5 of ARCHITECTURE).
- Sample discovery cards early — verdict + proof chain tell you *why*, every
  time.
- Watch disk: caches accumulate; `max_cache_size_gb` enforcement is deferred
  (documented in CHANGELOG legacy section).