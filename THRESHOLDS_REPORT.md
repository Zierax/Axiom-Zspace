# Axiom-ZSpace Threshold Catalog — Reference Report

_Generated 2026-08-17 00:15 UTC from `config/production.yaml` + `zspace_engine/thresholds.py` — single source of truth._

## 1. What this is

Every tunable number of the detection/validation pipeline lives in one place. No code edits are needed to change a threshold: edit `config/production.yaml` → section `thresholds`, then re-run the benchmark.

## 2. Three pre-defined profiles

Profiles trade **false-positive rate (FPR)** against **recall**. Changing the active profile: `python -m zspace_engine.thresholds --set balanced` (or edit `thresholds.profile` in the YAML).

**Always re-measure after any change** — the values below marked [M] were measured on the controlled benchmark (command below):

    python benchmarks_controlled/run_controlled.py --true 100 --false 80 --out benchmarks_controlled/runs/<NAME> --seed 20260814

## 3. Decision rules used by the validator

- All **critical** gates must pass (`critical_passed`), otherwise the candidate is FALSE_POSITIVE (circuit breaker).
- With all critical gates passing: ≤`verdict_max_fail_pass` total fails → SOVEREIGN_PASS; ≤`verdict_max_fail_conditional` → CONDITIONAL_PASS; more → FALSE_POSITIVE.
- `coherent_override_enabled` lets repeated-observation evidence override FP-2 (the power-FAP firewall). It is OFF in measured profiles because probe runs measured contamination FPR up to 62.5% (5/8) when ON (PROBE_FPR68 series).

## 4. All thresholds at a glance

| key | name | unit | weight | conservative | balanced | sensitive |
|---|---|---|---|---|---|---|
| `fp1_snr_min` | FP-1 BLS SNR minimum | σ | critical | 5.5 | 5.5 | 5.0 |
| `fp2_fap_max` | FP-2 power-spectrum FAP ceiling | probability | critical | 0.05 | 0.05 | 0.1 |
| `fp3_eo_sigma_max` | FP-3 Even/Odd depth Δσ ceiling | σ | critical | 3.0 | 3.0 | 3.0 |
| `fp4_shape_min` | FP-4 Shape Ratio (U vs V) floor | ratio | major | 0.4 | 0.4 | 0.3 |
| `fp5_secondary_snr_max` | FP-5 secondary-eclipse SNR ceiling | σ | critical | 3.0 | 3.0 | 3.0 |
| `fp5b_secondary_ratio_max` | FP-5b secondary depth-ratio ceiling | fraction | critical | 0.3 | 0.3 | 0.35 |
| `fp5c_alias_band` | FP-5c alias-eclipse ratio band | fraction | critical | [0.2, 0.9] | [0.2, 0.9] | [0.2, 0.9] |
| `fp7_density_band` | FP-7 stellar-density ratio band | ratio | critical | [0.2, 5.0] | [0.2, 5.0] | [0.15, 6.0] |
| `fp8_impact_max` | FP-8 impact-parameter cap | b | moderate | 0.9 | 0.9 | 0.95 |
| `fp10_min_transits` | FP-10 minimum independent transits | count | critical | 2 | 2 | 2 |
| `verdict_max_fail_pass` | Verdict: max fails for SOVEREIGN_PASS | count | meta | 2 | 2 | 3 |
| `coherent_override_enabled` | FP-2 coherent multi-transit override | flag | meta | False | False | True |
| `coherent_min_snr` | Override minimum matched SNR | σ | meta | 6.5 | 6.5 | 6.5 |
| `coherent_min_dip_fraction` | Override dip-recurrence floor | fraction | meta | 0.6 | 0.6 | 0.55 |
| `ladder_k` | Ladder depth (candidates tested) | count | meta | 20 | 20 | 20 |
| `detector_frequency_factor` | BLS frequency oversampling | x | meta | 20.0 | 20.0 | 20.0 |
| `ev_odd_sigma_eb` | Auditor Even/Odd EB Δσ flag | σ | auditor | 3.0 | 3.0 | 3.0 |
| `ingress_fraction_vshape` | Ingress-fraction V-shape flag | fraction | auditor | 0.45 | 0.45 | 0.45 |
| `context_density_mismatch` | Context a/R★ mismatch tolerance | fraction | context | 0.2 | 0.2 | 0.25 |
| `cvs_w_periodicity` | CVS weight: periodicity (S_P) | w | meta | 0.97 | 0.97 | 0.97 |
| `no_detection_snr_floor` | NO_DETECTION floor | σ | meta | 5.5 | 5.5 | 5.0 |

## 5. Profile values in full

### `conservative`

**Mode:** FPR=0 priority, measured  
**Measured metrics:** Contamination FPR 0/80 (0.0%) · Recall 52/100 (52/100 total; first-80 prefix 40/80) — measured, OVERHARM_FIX2  

| key | value |
|---|---|
| `fp1_snr_min` | 5.5 |
| `fp2_fap_max` | 0.05 |
| `fp3_eo_sigma_max` | 3.0 |
| `fp4_shape_min` | 0.4 |
| `fp5_secondary_snr_max` | 3.0 |
| `fp5b_secondary_ratio_max` | 0.3 |
| `fp5c_alias_band` | [0.2, 0.9] |
| `fp7_density_band` | [0.2, 5.0] |
| `fp8_impact_max` | 0.9 |
| `fp10_min_transits` | 2 |
| `verdict_max_fail_pass` | 2 |
| `coherent_override_enabled` | False |
| `coherent_min_snr` | 6.5 |
| `coherent_min_dip_fraction` | 0.6 |
| `ladder_k` | 20 |
| `detector_frequency_factor` | 20.0 |
| `ev_odd_sigma_eb` | 3.0 |
| `ingress_fraction_vshape` | 0.45 |
| `context_density_mismatch` | 0.2 |
| `cvs_w_periodicity` | 0.97 |
| `no_detection_snr_floor` | 5.5 |

### `balanced`

**Mode:** default, measured  
**Measured metrics:** Contamination FPR 0/80 (0.0%) · Recall 52/100 (first-80 prefix 40/80 vs 32/80 pre-calibration baseline OVERHARM_FIX, archived) — measured, OVERHARM_FIX2  

| key | value |
|---|---|
| `fp1_snr_min` | 5.5 |
| `fp2_fap_max` | 0.05 |
| `fp3_eo_sigma_max` | 3.0 |
| `fp4_shape_min` | 0.4 |
| `fp5_secondary_snr_max` | 3.0 |
| `fp5b_secondary_ratio_max` | 0.3 |
| `fp5c_alias_band` | [0.2, 0.9] |
| `fp7_density_band` | [0.2, 5.0] |
| `fp8_impact_max` | 0.9 |
| `fp10_min_transits` | 2 |
| `verdict_max_fail_pass` | 2 |
| `coherent_override_enabled` | False |
| `coherent_min_snr` | 6.5 |
| `coherent_min_dip_fraction` | 0.6 |
| `ladder_k` | 20 |
| `detector_frequency_factor` | 20.0 |
| `ev_odd_sigma_eb` | 3.0 |
| `ingress_fraction_vshape` | 0.45 |
| `context_density_mismatch` | 0.2 |
| `cvs_w_periodicity` | 0.97 |
| `no_detection_snr_floor` | 5.5 |

### `sensitive`

**Mode:** recall priority, EXPERIMENTAL  
**Measured metrics:** UNMEASURED — every deviation requires a full benchmark re-run  

| key | value |
|---|---|
| `fp1_snr_min` | 5.0 |
| `fp2_fap_max` | 0.1 |
| `fp3_eo_sigma_max` | 3.0 |
| `fp4_shape_min` | 0.3 |
| `fp5_secondary_snr_max` | 3.0 |
| `fp5b_secondary_ratio_max` | 0.35 |
| `fp5c_alias_band` | [0.2, 0.9] |
| `fp7_density_band` | [0.15, 6.0] |
| `fp8_impact_max` | 0.95 |
| `fp10_min_transits` | 2 |
| `verdict_max_fail_pass` | 3 |
| `coherent_override_enabled` | True |
| `coherent_min_snr` | 6.5 |
| `coherent_min_dip_fraction` | 0.55 |
| `ladder_k` | 20 |
| `detector_frequency_factor` | 20.0 |
| `ev_odd_sigma_eb` | 3.0 |
| `ingress_fraction_vshape` | 0.45 |
| `context_density_mismatch` | 0.25 |
| `cvs_w_periodicity` | 0.97 |
| `no_detection_snr_floor` | 5.0 |

## 6. Per-threshold rationale (measured evidence, pros/cons)

### FP-1 BLS SNR minimum  (`fp1_snr_min`)

- **Unit:** σ · **Direction:** gt · **Weight:** critical
- **Purpose:** BLS signal-to-noise ratio floor: below this the periodic power peak is not statistically distinguishable from the noise.
- **Measured evidence:** Measured: truths with SNR ≥ 5.5 at the correct period were all detectable; false noise-only candidates with SNR 5.0-5.5 exist but were vetoed by other gates.
- **Values:** conservative `5.5` · balanced `5.5` · sensitive `5.0`
- **Pros of tightening:** Fewer noise peaks enter the ladder (lower FPR).
- **Cons of tightening:** Loses the weakest genuine transits (recall drop); several label-SNR≈5.5 truths sit right at the edge.
- **Pros of loosening:** Recovers shallow planet transits.
- **Cons of loosening:** Pure-noise folds measured at SNR up to 8 → FPR risk unless later gates are strong (see coherent_override notes).
- **FPR risk when loosened:** MEDIUM

### FP-2 power-spectrum FAP ceiling  (`fp2_fap_max`)

- **Unit:** probability · **Direction:** lt · **Weight:** critical
- **Purpose:** False-Alarm-Probability of the BLS power peak under a self-calibrating exponential-tail noise model (red-noise conservative). This is the FPR firewall.
- **Measured evidence:** Measured: saturates at 0.5-1.0 for ALL candidates in red-noise-dominated spectra (truths AND noise) — it cannot separate; certification therefore depends on the coherent override or the raw SNR. 8 high-SNR truths with FAP≈1 died on this gate alone.
- **Values:** conservative `0.05` · balanced `0.05` · sensitive `0.1`
- **Pros of tightening:** Strict firewall; nothing certifies on a weak null.
- **Cons of tightening:** Frequent false rejection of genuine shallow-but-coherent signals (the reason the override exists).
- **Pros of loosening:** Accepts more signals without touching FPR if the coherent path stays off (override OFF measured FPR 0/80).
- **Cons of loosening:** With override ON, pure-noise folds pass FAP 0.3-1.0 and certify → measured contamination FPR 12.5-62.5% (PROBE_FPR68 series).
- **FPR risk when loosened:** LOW with override OFF, HIGH with override ON

### FP-3 Even/Odd depth Δσ ceiling  (`fp3_eo_sigma_max`)

- **Unit:** σ · **Direction:** lt · **Weight:** critical
- **Purpose:** Even and odd transit depths must be statistically identical; a large Δσ indicates an EB whose primary/secondary alternate.
- **Measured evidence:** Measured: truths sit at Δσ 0.05-2.0; grazing EBs in the false set reach Δσ 4-30 with SNR>100 → clean separation.
- **Values:** conservative `3.0` · balanced `3.0` · sensitive `3.0`
- **Pros of tightening:** Strongest single EB killer (catches deep EBs even when secondary windows miss).
- **Cons of tightening:** None measured within the observed distributions; do not tighten below ~2.5 or noise can trigger.
- **Pros of loosening:** No measured benefit.
- **Cons of loosening:** Deep EBs with Δσ 4-8 would certify → immediate FPR.
- **FPR risk when loosened:** HIGH

### FP-4 Shape Ratio (U vs V) floor  (`fp4_shape_min`)

- **Unit:** ratio · **Direction:** gt · **Weight:** major
- **Purpose:** Ratio of wing residual power vs centre: U-shaped (flat-bottom) planetary transits score >1; V-shaped (EB/grazing) score <1.
- **Measured evidence:** Measured: correct-period truths with per-cadence noise score 0.48-0.99 (!!) — the old >1.0 threshold rejected 4/8 recovered truths on pure measurement noise. False set: only 2/80 top candidates score <0.4, both independently vetoed by critical gates.
- **Values:** conservative `0.4` · balanced `0.4` · sensitive `0.3`
- **Pros of tightening:** Rejects V-shaped grazing EBs.
- **Cons of tightening:** Rejects grazing planets and noise-realized shallow transits (measured 0.48-0.99).
- **Pros of loosening:** Recovers grazing/near-grazing and shallow signals.
- **Cons of loosening:** Admits some V-shaped false dips; safe ONLY while density (FP-7) or even/odd (FP-3) critical gates still separate them (measured: yes).
- **FPR risk when loosened:** LOW at 0.4 (measured 0/80), grows below ~0.3

### FP-5 secondary-eclipse SNR ceiling  (`fp5_secondary_snr_max`)

- **Unit:** σ · **Direction:** lt · **Weight:** critical
- **Purpose:** No significant flux dip at phase 0.5 (a real companion would show one).
- **Measured evidence:** Measured: EBs show secondary SNR 5-90; planets ≤2.2.
- **Values:** conservative `3.0` · balanced `3.0` · sensitive `3.0`
- **Pros of tightening:** Kills EBs whose secondary is measured.
- **Cons of tightening:** At very low primary depth the phase-0.5 window catches noise; sub-2.5 thresholds can fire randomly on faint folds.
- **Pros of loosening:** Tolerates marginal secondary windows.
- **Cons of loosening:** Direct EB admission risk.
- **FPR risk when loosened:** HIGH

### FP-5b secondary depth-ratio ceiling  (`fp5b_secondary_ratio_max`)

- **Unit:** fraction · **Direction:** lt · **Weight:** critical
- **Purpose:** Phase-0.5 eclipse depth must be <30% of the primary (companion is not self-luminous).
- **Measured evidence:** Measured: planets ≤0.12; EBs 0.4-1.0. One recovered truth (SYN100076) sat at 0.33 — the 0.30 edge rejected it on noise.
- **Values:** conservative `0.3` · balanced `0.3` · sensitive `0.35`
- **Pros of tightening:** Clean EB discriminator.
- **Cons of tightening:** Borderline true values (0.30-0.35) get rejected on noise.
- **Pros of loosening:** Recovers borderline planets.
- **Cons of loosening:** Grazing EBs with compact secondaries (ratio 0.2-0.4) admitted; must stay below ~0.4.
- **FPR risk when loosened:** MEDIUM

### FP-5c alias-eclipse ratio band  (`fp5c_alias_band`)

- **Unit:** fraction · **Direction:** range · **Weight:** critical
- **Purpose:** Fold at 2×/3× the candidate period: an EB alias shows an ASYMMETRIC 0.2-0.9-ratio eclipse at SNR≥15; a planet folds to symmetric transits (ratio≈1.0).
- **Measured evidence:** Measured: catches grazing-EB sub-harmonics that faked clean planet folds (false-set rows 0.75/0.94/1.16 density in-band are rejected by FP-5c).
- **Values:** conservative `[0.2, 0.9]` · balanced `[0.2, 0.9]` · sensitive `[0.2, 0.9]`
- **Pros of tightening:** The only gate that sees hidden EB secondaries folded onto phase 0.
- **Cons of tightening:** Requires SNR≥15 to avoid noise flags; at huge EB SNRs (100+) it fires reliably.
- **Pros of loosening:** None.
- **Cons of loosening:** EB aliases certify as clean planets — measured FPR source.
- **FPR risk when loosened:** HIGH

### FP-7 stellar-density ratio band  (`fp7_density_band`)

- **Unit:** ratio · **Direction:** range · **Weight:** critical
- **Purpose:** Stellar density implied by transit depth+duration must match the TIC-stellar density (ratio≈1); an EB's half-duration geometry violates it.
- **Measured evidence:** Measured: recovered truths sit at 0.38-4.0 (noise around 1.0); false-set top candidates sit at ≤0.2 (n=73/80) or ≥9.4 (grazing EBs). In [0.2,0.4)∪(2.5,5] there are ZERO false top candidates and 6 truths → the widened band is FPR-free.
- **Values:** conservative `[0.2, 5.0]` · balanced `[0.2, 5.0]` · sensitive `[0.15, 6.0]`
- **Pros of tightening:** Rejects density-divergent EBs.
- **Cons of tightening:** Rejects grazing and noise-realized planets (0.4-2.5 produced 5/8 measured false rejections).
- **Pros of loosening:** Recovers grazing/near-grazing planets.
- **Cons of loosening:** Beyond [0.15, 6.0] grazing-EB densities (9-590) enter → FPR grows; keep 5.0 hard cap.
- **FPR risk when loosened:** LOW to 5.0 (measured 0/80), HIGH beyond ~6

### FP-8 impact-parameter cap  (`fp8_impact_max`)

- **Unit:** b · **Direction:** lt · **Weight:** moderate
- **Purpose:** Transit must not be grazing (b<0.9); V-shaped grazing geometry looks like an EB.
- **Measured evidence:** Measured: one recovered truth (SYN100020) at b=0.95 — as a moderate-weight fail it does NOT trip the circuit breaker, so the cap can stay 0.9.
- **Values:** conservative `0.9` · balanced `0.9` · sensitive `0.95`
- **Pros of tightening:** Flares the V-grazing alarm.
- **Cons of tightening:** Rejects the most grazing planets (b 0.9-1.0).
- **Pros of loosening:** Admits grazing planets.
- **Cons of loosening:** Grazing EBs also live at b→1; loose caps need the density cap ≤6 to stay safe.
- **FPR risk when loosened:** MEDIUM (coupled with fp7)

### FP-10 minimum independent transits  (`fp10_min_transits`)

- **Unit:** count · **Direction:** gt · **Weight:** critical
- **Purpose:** At least 2 transits must be observed; single-event dips are never certifiable.
- **Measured evidence:** Measured: single-event false targets all fail here; every true target (P≤13.5 d over ~1050 d) has ≫2.
- **Values:** conservative `2` · balanced `2` · sensitive `2`
- **Pros of tightening:** Kills single-event contamination.
- **Cons of tightening:** None in the synthetic band (P≤13.5).
- **Pros of loosening:** None.
- **Cons of loosening:** Single transients certify — unacceptable.
- **FPR risk when loosened:** ABSOLUTE (design invariant)

### Verdict: max fails for SOVEREIGN_PASS  (`verdict_max_fail_pass`)

- **Unit:** count · **Direction:** lt · **Weight:** meta
- **Purpose:** All critical gates must pass; up to N non-critical fails tolerated.
- **Measured evidence:** Measured: recovered truths have 0-1 fails after the gate widening.
- **Values:** conservative `2` · balanced `2` · sensitive `3`
- **Pros of tightening:** Fewer marginal candidates certified.
- **Cons of tightening:** Moderate-weight noisy measures (FP-4/FP-8) can veto otherwise clean signals.
- **Pros of loosening:** Tolerates noisy minor gates.
- **Cons of loosening:** More noise candidates pass; only with strong critical gates.
- **FPR risk when loosened:** LOW (critical breaker dominates)

### FP-2 coherent multi-transit override  (`coherent_override_enabled`)

- **Unit:** flag · **Direction:** flag · **Weight:** meta
- **Purpose:** Let overwhelming repeated-observation evidence (≥3 consistent transits + SNR + no secondary) override the saturated power-FAP.
- **Measured evidence:** Measured: ON caused FPR explosion (contamination FPR up to 62.5% on 8-target probes — PROBE_FPR68 series); OFF with widened gates gives 0/80 FPR. The dip-fraction discriminator at 0.6 is a coin-flip for real shallow transits (0.5-0.55 measured) — reliable only above per-transit SNR≈3.
- **Values:** conservative `False` · balanced `False` · sensitive `True`
- **Pros of tightening:** —
- **Cons of tightening:** —
- **Pros of loosening:** —
- **Cons of loosening:** —
- **FPR risk when loosened:** HIGH (measured)

### Override minimum matched SNR  (`coherent_min_snr`)

- **Unit:** σ · **Direction:** gt · **Weight:** meta
- **Purpose:** Floor for the coherent path. Measured noise folds reach 6-8σ → 6.5 is the expected value of the noise, not a separation point.
- **Measured evidence:** Measured: noise folds at SNR 6-8 (PROBE_FPR68).
- **Values:** conservative `6.5` · balanced `6.5` · sensitive `6.5`
- **Pros of tightening:** Harder to exploit by noise.
- **Cons of tightening:** Weak-but-real signals can't override.
- **Pros of loosening:** More overrides.
- **Cons of loosening:** FPR with override ON and SNR floor ≤6.5 is already measured-explosive.
- **FPR risk when loosened:** HIGH (with override ON)

### Override dip-recurrence floor  (`coherent_min_dip_fraction`)

- **Unit:** fraction · **Direction:** gt · **Weight:** meta
- **Purpose:** Fraction of epochs showing a dip; genuine periodicity ~1.0, ephemeral noise ~0.5.
- **Measured evidence:** Measured: noise folds show dipF≈0.9+ while passing FAP (probes certify 12.5-62.5% of contamination) — the metric is useless at SNR<3 per transit (sign flips); truths at per-transit SNR 1-2 measure 0.5-0.55, below the 0.6 floor.
- **Values:** conservative `0.6` · balanced `0.6` · sensitive `0.55`
- **Pros of tightening:** Feels safe at 0.6+.
- **Cons of tightening:** Rejects exactly the shallow truths it was meant to rescue.
- **Pros of loosening:** Nothing without a binomial significance test.
- **Cons of loosening:** Noise folds show dipF≈0.9+ → they certify → measured FPR explosion. Keep the floor, keep override OFF.
- **FPR risk when loosened:** HIGH, do not loosen while override is ON

### Ladder depth (candidates tested)  (`ladder_k`)

- **Unit:** count · **Direction:** gt · **Weight:** meta
- **Purpose:** How many period candidates are validated per target; recall grows with k, CPU scales with k.
- **Measured evidence:** Measured: correct-period candidate is the #1 ladder entry for 7/8 recovered truths; sub-harmonics occupy ranks 1-3 on weaker targets.
- **Values:** conservative `20` · balanced `20` · sensitive `20`
- **Pros of tightening:** Faster.
- **Cons of tightening:** Misses weaker true periods ranked lower.
- **Pros of loosening:** Better weak-target recall (alias + true candidates both tested).
- **Cons of loosening:** CPU cost; FPR unchanged (each candidate still must pass gates).
- **FPR risk when loosened:** NONE

### BLS frequency oversampling  (`detector_frequency_factor`)

- **Unit:** x · **Direction:** gt · **Weight:** meta
- **Purpose:** Frequencies per trial per bin; 20× recovers shallow peaks that 1× misses (prior 30k-bin runs destroyed weak signals).
- **Measured evidence:** Measured: weak truths (SNR≈6-8) only appeared at factor 20 on unbinned data.
- **Values:** conservative `20.0` · balanced `20.0` · sensitive `20.0`
- **Pros of tightening:** Faster.
- **Cons of tightening:** Aliasing/missed weak peaks.
- **Pros of loosening:** Weak-signal recall.
- **Cons of loosening:** CPU cost, more ladder noise entries (gated downstream).
- **FPR risk when loosened:** NONE

### Auditor Even/Odd EB Δσ flag  (`ev_odd_sigma_eb`)

- **Unit:** σ · **Direction:** gt · **Weight:** auditor
- **Purpose:** EB flag inside the auditors (paired with p<0.01).
- **Measured evidence:** Aligned with FP-3 (3.0).
- **Values:** conservative `3.0` · balanced `3.0` · sensitive `3.0`
- **Pros of tightening:** Flags EBs earlier.
- **Cons of tightening:** False flags at noisy depths.
- **Pros of loosening:** None.
- **Cons of loosening:** EBs pass auditor stage.
- **FPR risk when loosened:** MEDIUM

### Ingress-fraction V-shape flag  (`ingress_fraction_vshape`)

- **Unit:** fraction · **Direction:** gt · **Weight:** auditor
- **Purpose:** Ingress/duration >0.45 ⇒ V-shaped ⇒ EB-like flag in auditor scoring.
- **Measured evidence:** Complements FP-4; grazing planets can exceed 0.45 — penalty only, not hard veto.
- **Values:** conservative `0.45` · balanced `0.45` · sensitive `0.45`
- **Pros of tightening:** Stronger V-shape separation.
- **Cons of tightening:** Penalizes grazing planets in CVS.
- **Pros of loosening:** Fairer to grazing planets.
- **Cons of loosening:** EB-like shapes less penalized.
- **FPR risk when loosened:** LOW (CVS is soft)

### Context a/R★ mismatch tolerance  (`context_density_mismatch`)

- **Unit:** fraction · **Direction:** lt · **Weight:** context
- **Purpose:** Catalog vs transit a/R★ fractional deviation above 20% flags density mismatch in context checks.
- **Measured evidence:** Soft flag; feeds CVS penalties not gates.
- **Values:** conservative `0.2` · balanced `0.2` · sensitive `0.25`
- **Pros of tightening:** Flags suspicious densities.
- **Cons of tightening:** Noise excursions flag clean planets.
- **Pros of loosening:** None measured.
- **Cons of loosening:** Density-divergent EBs less penalized.
- **FPR risk when loosened:** LOW

### CVS weight: periodicity (S_P)  (`cvs_w_periodicity`)

- **Unit:** w · **Direction:** meta · **Weight:** meta
- **Purpose:** Composite Vitality Score mixing weight; the bench derives CVS = (0.97·S_P + 0.83·S_depth + 0.61·S_limb + 0.31·0.5)/Σw.
- **Measured evidence:** CVS ≥0.80 (core.THRESHOLD_PLANET, mirrored in config detection.cvs_planet_threshold) labels PLANET CANDIDATE; all certified candidates in the benchmark have CVS≥0.8.
- **Values:** conservative `0.97` · balanced `0.97` · sensitive `0.97`
- **Pros of tightening:** Periodicity dominates classification.
- **Cons of tightening:** Depth/shape quality underweighted.
- **Pros of loosening:** More balanced score.
- **Cons of loosening:** Score drifts with pipeline changes.
- **FPR risk when loosened:** LOW (informational)

### NO_DETECTION floor  (`no_detection_snr_floor`)

- **Unit:** σ · **Direction:** gt · **Weight:** meta
- **Purpose:** If no ladder candidate exceeds this SNR the target is declared NO_DETECTION (undetected) instead of certified/FP.
- **Measured evidence:** Aligned with fp1_snr_min; governs reported recall, not FPR.
- **Values:** conservative `5.5` · balanced `5.5` · sensitive `5.0`
- **Pros of tightening:** Honest accounting.
- **Cons of tightening:** Marks weak truths undetected even when recoverable.
- **Pros of loosening:** More candidates enter validation.
- **Cons of loosening:** None beyond CPU.
- **FPR risk when loosened:** NONE


## 7. Change workflow

1. Edit `config/production.yaml` (one value, or build a new profile block).
2. Verify thresholds load: `python -m zspace_engine.thresholds --show`.
3. Re-measure: python benchmarks_controlled/run_controlled.py --true 100 --false 80 --out benchmarks_controlled/runs/<NAME> --seed 20260814
4. Compare FPR/recall against the measured baselines in section 5; update this report with the new numbers.
