# Glossary

Terms used across the pipeline, docs, and benchmarks — in the order you will
meet them in the ARCHITECTURE walkthrough.

| Term | Meaning |
|---|---|
| **TESS / Kepler** | Space telescopes producing the long-cadence light curves the pipeline ingests (TESS ~30-min, Kepler ~30-min LC). |
| **FITS** | File format of telescope light curves (time, flux, quality columns). |
| **Quality flag** | Per-cadence data-quality bitmask; the pipeline retains only `quality == 0` cadences. |
| **MAST** | Mikulski Archive for Space Telescopes — where light curves are downloaded from (via `lightkurve`). |
| **NEA** | NASA Exoplanet Archive — the truth source for the real benchmark and for cross-checks. |
| **BLS** | Box Least-Squares — period-search algorithm that fits a box-shaped dip at each trial period; the core detection statistic. |
| **SNR** | Signal-to-noise ratio of the detected dip (transit depth ÷ noise). Gate FP-1 requires ≥ 5.5. |
| **FAP** | False-Alarm Probability of the winning BLS peak under the self-calibrating exponential-tail noise model (red-noise conservative). The FP-2 firewall; saturates near 1.0 under red noise for truths and noise alike. |
| **Red noise** | Correlated (non-white) noise — stellar variability/systematics; the dominant noise regime in these light curves, which the FAP model must handle conservatively. |
| **Harmonic / alias** | A period that is an integer multiple/divisor of the true period where BLS power leaks (e.g. true 10 d also shows power at 5 d, 20 d…). The main class of wrong-ephemeris errors. |
| **Period prior / ladder** | `period_prior_days` hint searched over `k` harmonics ("ladder"); used only when the hinted peak carries ≥ 5% of the global peak SNR — prior, not lock. |
| **Fold** | Phase-folding the light curve at a trial period to stack transits. |
| **Ephemeris** | (period, epoch/t0, duration) description of a transit signal; the identity gate compares the recovered ephemeris to archive truth. |
| **t0** | Transit midtime (epoch) in days (BKJD). |
| **Gate** | A falsification test with a threshold, weight, and evidence string — FP-1 … FP-10 (see ARCHITECTURE §6.2). |
| **FP-1 … FP-10** | The false-positive gate chain: SNR floor, FAP ceiling, even/odd Δσ, shape ratio, secondary eclipse (+ratio/alias), density band, impact parameter, minimum transits. |
| **Circuit breaker** | Rule: any *critical* gate FAIL ⇒ verdict is FALSE_POSITIVE; SOVEREIGN_PASS is impossible (tested invariant). |
| **Veto** | A critical-gate failure that additionally caps the CVS below 0.35 so classification cannot rise above FALSE POSITIVE. |
| **CVS** | Composite Vitality Score — weighted blend of periodicity/depth/limb/stellar scores (w = 0.97/0.83/0.61/0.31). |
| **Classify tiers** | CVS ≥ 0.80 PLANET CANDIDATE · ≥ 0.55 LIKELY · ≥ 0.35 AMBIGUOUS · else FALSE POSITIVE. |
| **SOVEREIGN_PASS** | Highest verdict: all critical gates passed and ≤ 2 non-critical fails. |
| **CONDITIONAL_PASS** | All critical gates passed and ≤ 3 non-critical fails. |
| **FALSE_POSITIVE** | A critical gate failed, or fails exceeded the allowance. |
| **NO_DETECTION** | No candidate exceeded the SNR floor — the target was not detected. |
| **Proof chain** | The accumulating, human-readable list of every gate/audit result that produced a verdict; every card ships one. |
| **Conflict** | Recorded disagreement between gates at high SNR (e.g. strong peak vs density mismatch); logged explicitly. |
| **EB** | Eclipsing binary — the main false-positive class (primary/secondary eclipse alternation). |
| **Grazing EB** | EB whose sub-harmonics can mimic shallow planet-like transits; targeted by FP-5c and even/odd tests. |
| **U / V shape** | Transit ingress-egress morphology: planets ≈ U (flat-bottom), grazing events ≈ V. Ratio < 0.4 → FP-4 fail. |
| **Density (FP-7)** | Transit-derived stellar density ÷ TIC stellar density; planets fall in band [0.2, 5.0]. |
| **Impact parameter b** | Transit chord geometry; b ≤ 0.9 (FP-8) excludes grazing. |
| **TPF** | Target Pixel File (TESS/Kepler pixels) — used by the optional centroid test. |
| **logg** | Stellar surface gravity (cgs) — context metadata for density checks. |
| **SIMBAD / Gaia DR3** | External catalogs consulted by `check_external_catalogs` for multiplicity/duplicity flags (best-effort, offline-tolerant). |
| **DR25** | Kepler Data Release 25 stellar delivery — source of the quiet-star list (`nkoi=0` AND `nconfp=0`). |
| **TIC / KIC** | TESS Input Catalog / Kepler Input Catalog star identifiers. |
| **Discovery card** | The output JSON (`Discovery_…json`) with detection params, audit results, proof chain, verdict, and enrichment (radius, Teq, a/R★…). |
| **zspace_id** | Internal candidate identifier `ZS-T-‹TIC›-NN` (TIC + planet ordinal) glued to a target across outputs. |
| **Contamination FPR** | Fraction of false targets wrongly certified (measured 0/80). |
| **Recall@target** | Fraction of true targets certified with the period matching truth within the acceptance window (5% real, tighter controlled). |
| **Proxy FPR** | Quiet-star certification rate used as a stand-in FPR (measured 33.3% — quiet stars are not guaranteed planet-free). |
| **Wrong-ephemeris** | Certification at the right target but a *wrong* period (real signal, wrong harmonic); counted separately from contamination. |
| **Override** | `coherent_override_enabled` — lets coherent repeated transits override an FP-2 fail; **OFF by default** (measured FPR explosion when ON). |
| **Architect / reviewer terms** | "Measured" = backed by an evidence run in `benchmarks_*/evidence/` · "Legacy" = v2.x history without preserved artifacts (CHANGELOG §LEGACY). |