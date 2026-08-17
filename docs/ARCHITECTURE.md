# Architecture

This is the reference document for how Axiom-ZSpace works: the stages a light
curve goes through, the modules and classes that implement each stage, the
decision rules, and the data contracts between them. It is written to be
readable top-to-bottom; `docs/GLOSSARY.md` defines every acronym used here.

---

## 1. Purpose and scope

Axiom-ZSpace is a **blind-search transit detection pipeline** for TESS/Kepler
long-cadence light curves. "Blind" means: the pipeline searches for periodic
dips without being pointed at a known ephemeris. A period-PRIOR (not period
lock) exists for benchmark use — it *prefers* a BLS peak near a hint
(fundamental + harmonics) only when the hinted peak carries enough power, and
falls back to the global peak otherwise.

The pipeline's outputs for one target:

| Output | Format | Where |
|---|---|---|
| Detection + validation proof | JSON discovery card (`Discovery_ZS-T-….json`) | `axiom_output/…` (git-ignored) |
| Human-readable report | MD rendered from the card | `axiom_output/…` |
| Batch summary | `summary.json` per sector | `axiom_output/sector_N/` |
| Threshold reference | `THRESHOLDS_REPORT.md` (root) | auto-generated |

---

## 2. Pipeline stages at a glance

```
light curve (TESS/Kepler FITS, or synthetic, or raw arrays)
        │
        ▼
[1] ingestion.py        fetch / cache / quality mask / sigma-clip /
                        normalize / Savitzky-Golay flatten
        │
        ▼
[2] detectors.py        BLS periodogram (frequency + duration grid), FAP
                        under a self-calibrating exponential-tail noise model,
                        period-prior peak selection, harmonic rejection
        │
        ▼
[3] ephemeris.py        fold-based period/t0/duration refinement; alias
                        comparator; EphemerisResolver with harmonics ladder
        │
        ▼
[4] auditors.py         physical-invariant audits on the folded transit:
                        even/odd depth (Welch t-test), per-transit depth
                        consistency, U/V limb shape ratio, ingress/egress
                        duration ratio, secondary-eclipse search, MCMC
                        posterior fit (optional), stellar-density consistency
        │
        ▼
[5] context.py          stellar metadata (TIC), centroid-shift test (TPF,
                        optional), secondary-eclipse context, multi-sector
                        consistency (optional)
        │
        ▼
[6] validator.py        ProofEngine: let the FP gates vote; circuit breaker;
                        external-catalog cross-check (SIMBAD/Gaia); ephemeris
                        identity gate against archive truth; conflict
                        detection & logging; final verdict
        │
        ▼
[7] core.py / report.py CVS (Composite Vitality Score) classification into 4
                        tiers; enrichment (orbital mechanics, equilibrium
                        temperature, planet radius); discovery card JSON
        │
        ▼
    output_organizer.py route by verdict → axiom_output/… (kept / rejected)
```

Stages 1–7 all run inside a single `AxiomValidator.process_*` call
(`run_pipeline.py` orchestrates it), or batch-wise via `SectorProcessor`
(which loops the same engine per target).

---

## 3. Module register

| Module | Key public API | Responsibility |
|---|---|---|
| `ingestion.py` | `FITSCacheManager`, `LightCurveIngester` (`fetch`, `process`, `from_fits`, `from_arrays`), `LightCurveProduct` | MAST download with disk cache; quality-mask rejection; sigma clipping; normalization; flattening; physical consistency assertion on products |
| `detectors.py` | `BLSDetector` (`run`), `BLSResult`, `FAPValidator` | Box-least-squares search over frequency × duration grid; period-prior ladder (`k=20` harmonics, min relative SNR 5%); FAP from self-calibrating exponential-tail model (red-noise conservative) |
| `ephemeris.py` | `EphemerisResolver`, `FoldSignature`, `ResolvedEphemeris` | Fold at trial period, merge dip signatures, resolve alias ambiguity (harmonics ladder + spacing+equality heuristics), refine t0/duration |
| `auditors.py` | `TransitAuditor` (`even_odd_test`, `depth_consistency_score`, `limb_shape_score`, `secondary_eclipse_test`, `ingress_egress_test`), `MCMCValidator`, result dataclasses | Physics audits; per-audit result dataclasses feeding the proof chain |
| `context.py` | `TICMetadataFetcher`, `StellarContextAuditor`, `StellarDensityFilter`, `CentroidShiftTest`, `SecondaryEclipseSearch`, `MultiSectorConsistencyCheck` | Stellar parameters, density-ratio gate, centroid (TPF), context flags |
| `validator.py` | `ProofEngine`, `AxiomValidator`, `ArchiveQueryEngine`, `PeriodComparator`, `check_external_catalogs` | The gate engine (see §6); verdict assembly; archive cross-checks; per-target `ValidationResult` |
| `core.py` | `CompositeVitalityScore`, `VitalityMatrix`, `apply_hard_filters`, `OrbitalMechanics`, planet radius / equilibrium-temperature helpers | CVS scoring + 4-tier classification (§7); hard physical filters (SNR, depth, radius sanity) |
| `chi_squared.py` | `chi2`-based goodness-of-fit analysis | Used by auditors for model-vs-data quality assessment |
| `report.py` | `TruthimaticsReport` | Discovery card JSON with full proof chain |
| `thresholds.py` | `threshold(key)`, `profile_values(profile)`, `catalog()`, `summary()`, CLI (`--show/--profile/--set`) | Single source of truth for every tunable number (see §8) |
| `thresholds_report.py` | CLI module | Regenerates `THRESHOLDS_REPORT.md` from live config — the human reference for every value + its measured evidence |
| `sector_processor.py` | `SectorProcessor` | Batch mode: loop targets of a TESS sector; per-sector `summary.json` + `discoveries.json`; compact per-planet verdict records |
| `output_organizer.py` | `OutputOrganizer` | Route outputs by verdict to per-sector directories |
| `config.py` | `load_config`, `cvs_planet_threshold()`, `fap_threshold()`, `snr_threshold()`, `snr_ref()` | Config resolution (env override `AXIOM_CONFIG`, then `config/production.yaml`, then bundled defaults) |
| `logging_config.py` | `get_logger`, `suppress_astroquery_logger` | Central logging (`axiom_pipeline.log` + console) |
| `constants.py` | `G_SI`, `M_SUN`, `R_SUN`, `L_SUN`, `SIGMA_SB`, `AU`, `R_EARTH_SOLAR` | IAU 2015 Resolution B3 physical constants |

---

## 4. Stage detail

### 4.1 Ingestion (`ingestion.py`)

- `LightCurveIngester.fetch()` downloads TESS/Kepler light curves via
  `lightkurve`/MAST. `FITSCacheManager` keeps a disk cache
  (`fits_cache.cache_dir`, default `.cache/fits`, git-ignored):
  `is_cached → get_cached_fits → cache_fits → download_and_cache`.
- Quality flags are applied (`quality == 0` retained), flux is sigma-clipped,
  normalized to ~1.0, and long-term systematics removed with a
  Savitzky-Golay flatten (window configurable via `process()`).
- `LightCurveProduct` asserts physical consistency (finite flux, sane time
  spacing, non-degenerate baseline) before the pipeline proceeds.

### 4.2 Detection (`detectors.py`)

- `BLSDetector` builds a frequency grid (period band `[0.5, 13.5]` d, grid
  resolution governed by `detector_frequency_factor = 20` samples per
  fundamental frequency cell) and a duration grid, then computes the BLS
  statistic at each cell.
- **Period-prior selector** (benchmark / targeted search): if
  `period_prior_days` is given, the candidate ladder is the prior ± k
  harmonics with `ladder_k = 20`; a ladder peak is chosen over the global peak
  only if its SNR is at least `ladder_min_relative_snr = 0.05` × the global
  peak SNR. Otherwise the global peak wins — the detector stays blind.
- Harmonic rejection removes peaks that are integer multiples of a stronger
  fundamental (the “opening harmonics” where BLS energy leaks).
- `FAPValidator` computes the false-alarm probability of the winning peak
  against a **self-calibrating exponential-tail noise model** — this is the
  FP-2 firewall and the single most contested number in the pipeline (see
  THRESHOLDS_REPORT.md §3, and the override discussion in §6.4).

### 4.3 Ephemeris (`ephemeris.py`)

- Folds the light curve at the candidate period, identifies dip signatures
  (`FoldDip`), checks for equal-depth alternating events (EB tell-tale) and
  spacing presence; `EphemerisResolver` compares the fundamental against the
  harmonic ladder with an alias comparator and resolves period/t0/duration
  ambiguities. Output is a `ResolvedEphemeris`.
- **Ephemeris identity gate** (`validator.py`): if an archive match exists
  (NEA), the reported period/epoch are compared; a mismatch demotes the
  candidate and is logged as "wrong-ephemeris" rather than re-passing with
  the archive number — this kills the class of wrong-ephemeris
  `SOVEREIGN_PASS` bugs fixed in v2.x lineage.

### 4.4 Audits (`auditors.py`)

All audits run on the folded, phase-folded transit and return dataclasses that
the proof engine consumes:

| Audit | Result type | What it catches |
|---|---|---|
| Even/odd depth | `EvenOddResult` | EB whose primary/secondary alternate (Δσ > 3.0 → gate FP-3) |
| Per-transit depth consistency | `DepthConsistencyResult` | Depth scatter from rotation/spot mismodeling |
| U/V limb shape | `LimbShapeResult` | V-shaped grazing events vs U-shaped planet transits (ratio < 0.4 → FP-4) |
| Ingress/egress | `IngressEgressResult` | Duration-ratio anomalies |
| Secondary eclipse | `SecondaryEclipseResult` | Phase-0.5 eclipse → EB (FP-5, FP-5b/5c) |
| MCMC (optional, `--mcmc`) | `MCMCResult` | Full posterior transit fit (expensive; off by default) |
| Density consistency | via `context.py` | Transit density vs TIC stellar density mismatch (FP-7 band [0.2, 5.0]) |

### 4.5 Context (`context.py`)

`TICMetadataFetcher` retrieves stellar parameters (radius, mass, Teff, logg)
needed by the density gate and the CVS “stellar” component. Optional layers:
`CentroidShiftTest` (TPF pixel centroid displacement, `--tpf-centroids`),
`MultiSectorConsistencyCheck` (`--multi-sector`), all guarded by config flags
(`use_tpf_centroids`, `check_multi_sector`) that default to **false** in
production for cost and determinism.

---

## 5. Data contracts between stages

- **Detection → Ephemeris**: `BLSResult(period, snr, fap, duration, t0, power)`.
- **Ephemeris → Auditors**: `ResolvedEphemeris(period, t0, duration, alias_confidence)`.
- **Auditors → Validator**: result dataclasses (`EvenOddResult`, …), each with
  `passed: bool` and an evidence string.
- **Validator → Report**: `ValidationResult(gate_results: dict, verdict,
  proof_chain: list[str], conflicts: list[str])` — every verdict leaves a
  proof chain; nothing is decided silently.

---

## 6. The proof engine (validator)

### 6.1 Verdict tiers

| Verdict | Meaning |
|---|---|
| `SOVEREIGN_PASS` | All critical gates passed and ≤ `verdict_max_fail_pass` (2) non-critical fails |
| `CONDITIONAL_PASS` | All critical gates passed, ≤ `verdict_max_fail_conditional` (3) non-critical fails |
| `FALSE_POSITIVE` | Any critical gate failed (circuit breaker) OR more fails than allowed |
| `NO_DETECTION` | No candidate exceeded `no_detection_snr_floor` (5.5) |

### 6.2 The gates (threshold catalog keys → measured profile values)

| Gate | Key | Direction | Value (balanced) | Role |
|---|---|---|---|---|
| FP-1 | `fp1_snr_min` | SNR ≥ | 5.5 σ | Weak peaks never certify |
| FP-2 | `fp2_fap_max` | FAP ≤ | 0.05 | **Critical firewall** (see §6.4) |
| FP-3 | `fp3_eo_sigma_max` | Δσ ≤ | 3.0 | Even/odd EB killer |
| FP-4 | `fp4_shape_min` | shape ≥ | 0.4 | V-shaped artifact killer |
| FP-5 | `fp5_secondary_snr_max` | ≤ | 3.0 | Phase-0.5 eclipse |
| FP-5b | `fp5b_secondary_ratio_max` | ≤ | 0.30 | Eclipse depth ratio |
| FP-5c | `fp5c_alias_band` + `fp5c_alias_min_snr` | band | [0.20, 0.90] @ ≥15σ | Alias-eclipse sub-harmonics |
| FP-7 | `fp7_density_band` | band | [0.2, 5.0] | Transit/TIC density mismatch |
| FP-8 | `fp8_impact_max` | b ≤ | 0.9 | Grazing geometry |
| FP-10 | `fp10_min_transits` | ≥ | 2 | Single-event dips never certify |

Every gate has a weight (`critical` / `major`), a direction, a purpose string,
and **measured-evidence text with pros/cons of tightening/loosening** — all
rendered by `THRESHOLDS_REPORT.md`.

### 6.3 Circuit breaker

`SOVEREIGN_PASS` is literally impossible while any critical gate reports FAIL
(`validator.py` — asserted by `tests/test_task_5_1_circuit_breaker.py`). A
critical fail pins the verdict to `FALSE_POSITIVE` regardless of other scores.

### 6.4 The FP-2 firewall and the coherent override

- FP-2's FAP **saturates at ≈1.0 for everything under red noise** — truths and
  noise alike. On its own it cannot separate; the pipeline's measured
  separation comes from combining it with the other gates and the override.
- `coherent_override_enabled` allows repeated-observation coherent evidence
  (min SNR 6.5, min 3 transits, dip-fraction floor 0.6) to override an FP-2
  fail. **Measured result: it is OFF in the shipping profiles.** The
  `PROBE_FPR68` series measured contamination FPR of 62.5% (5/8), 12.5% (1/8)
  and 37.5% (3/8) on tiny 8-noise-target probes when the override is ON —
  noise folds were passing the firewall. See `docs/BENCHMARKS.md` and
  THRESHOLDS_REPORT.md §3.

### 6.5 Conflict detection & external catalogs

- **Conflict detection/logging**: when gates disagree with high SNR
  (e.g. strong peak vs density mismatch vs V-shape), a conflict record is
  appended to the proof chain and logged (`tests/test_task_7_1_*`,
  `test_task_7_2_conflict_logging.py`).
- `check_external_catalogs(tic_id)` queries SIMBAD (multiplicity keywords:
  Double, Multiple, EB*, V*, Binary, SB*, El*) and Gaia DR3 — best-effort,
  network tolerant, never required for a verdict offline.

---

## 7. CVS scoring and classification (`core.py`)

`CompositeVitalityScore` mixes four component scores:

```
CVS = (w_P·S_P + w_depth·S_depth + w_limb·S_limb + w_stellar·S_stellar) / Σw
```

with measured weights `w = (0.97, 0.83, 0.61, 0.31)` and the stellar
component scored from context metadata. Classification (4 tiers):

| CVS | Verdict string |
|---|---|
| ≥ 0.80 | `PLANET CANDIDATE` |
| ≥ 0.55 | `LIKELY PLANET CANDIDATE` |
| ≥ 0.35 | `AMBIGUOUS / REQUIRES FOLLOW-UP` |
| < 0.35 | `FALSE POSITIVE` |

A critical-gate veto caps CVS below 0.35 so a vetoed candidate can never
classify above `FALSE POSITIVE` (proof-chain entry: `VETO: …`).

`apply_hard_filters` and `OrbitalMechanics` add physical bounds: SNR floor,
minimum depth, sane radius range, semi-major axis, equilibrium temperature —
computed with IAU 2015 constants and used for the enrichment fields in the
card (radius in Earth units, etc.).

---

## 8. The threshold catalog (single source of truth)

- **Canonical storage**: `config/production.yaml` (section `thresholds:`,
  three profiles) — parsed by `zspace_engine/thresholds.py` at import;
  `zspace_engine/config.py` holds only fallback defaults.
- **The rule**: every tunable number of detection + validation lives in the
  catalog. Dev width like “9-gate engine”, “0.45 V-shape ingress fraction” —
  all catalog keys.
- **Profiles**:

| Profile | Philosophy | Status |
|---|---|---|
| `conservative` | FPR=0 priority, identical values to balanced | measured (same run) |
| `balanced` | **default** — shipping values | measured (OVERHARM_FIX2) |
| `sensitive` | recall priority (looser FAP/shape/density, override ON) | **EXPERIMENTAL — NOT measured** |

- `THRESHOLDS_REPORT.md` (root) is auto-generated by
  `python -m zspace_engine.thresholds_report` — it carries per-key evidence
  statements, the generated profile tables, and the measured metrics. It must
  be regenerated after any catalog edit (same command, committed together).

---

## 9. Outputs and routing

- `SectorProcessor` runs the same engine over every target of a sector,
  writing `summary.json`/`discoveries.json` per sector under
  `output.base_directory` (default `axiom_output/`).
- `OutputOrganizer` routes per-verdict: certified candidates → discovery
  cards; `FALSE_POSITIVE` → `rejected/` folders.
- Cards are `Discovery_ZS-T-….json` (git-ignored — the output tree is
  regenerated per run; only `benchmarks_*/evidence/` is versioned).

---

## 10. Design invariants (asserted by tests)

1. Thresholds are data (`tests/test_period_prior_selector.py` calibrates the
   selector with explicit prior-floor assertions; catalog tests pin values).
2. Determinism: same seed ⇒ same sample ⇒ same results
   (`tests/test_reproducibility.py`, 7 tests).
3. Circuit breaker: no critical-fail SOVEREIGN_PASS (`test_task_5_1_*`).
4. Ephemeris identity: wrong ephemeris never certified (`tests/test_ephemeris_identity_gate.py`).
5. Calibrated gates: FP-4 floor 0.4 and FP-7 band [0.2, 5.0] are exercised on
   synthetic true/false distributions (`test_task_4_1_*`, `test_task_4_2_*`).
6. No tuning without measurement (`docs/CONTRIBUTING.md` — the one rule).