# Axiom-ZSpace Real-Data Honesty Evaluation (Kepler / NEA)

**Generated:** 2026-08-16 08:46:40 UTC
**Run:** REAL_FINAL
**Started:** finalize
**Method:** identical blind pipeline to the controlled bench (BLS +
ladder + auditors + sovereign validator, archive stub OFFLINE); input =
real Kepler long-cadence light curves (Q4-Q9 ~540 d baseline).

## Key Metrics

| Metric | Value |
|---|---|
| Confirmed planets (NEA truth, 1.0-13.5 d) | 12 |
| Quiet stars (zero KOI proxy) | 12 |
| **Recall@target period** (certified & P within 5% of NEA period) | 41.7% (5/12) |
| **Recall@any known planet of host** (certified & P matches another known planet) | 50.0% (6/12) |
| **Quiet-star certification rate** (proxy FPR) | 33.3% (4/12) |
| **Total FP rate** (incl. wrong-ephemeris certs on true hosts) | 46.7% (7 FP / 15) |
| **Precision (target-level)** | 41.7% |
| Confusion | TP(target)=5 TP(other known)=1 FP(quiet)=4 FP(unknown)=3 TN=8 |

*Controlled-bench comparison (for context): recall 52%, FPR 0% on
  synthetic noise/contamination targets. Real data is harder:
  long cadence, real systematics, and multiple real signals per
  system (e.g. Kepler-37's outer planets alias into the band).*

## Target list (true)

| KIC | star | NEA P (d) | R_p (Re) | found P | perr% | match | status |
|---|---|---|---|---|---|---|---|
| 8478994 | Kepler-37 b / Kepler-37 c / Kepler-37 d / Kepler-37 e | 13.36702000000 | 0.30980000 | 7.958303378864127 | 40.5 | none | OFFLINE_NEW_DISCOVERY |
| 8073705 | Kepler-1445 b | 10.60052251000 | 0.96000000 | 11.558105963268083 | 9.0 | none | OFFLINE_NEW_DISCOVERY |
| 4814502 | Kepler-1933 b | 4.94329023400 | 1.17984089 | 6.139448350637403 | 24.2 | none | FALSE_POSITIVE |
| 9285568 | Kepler-1575 b | 2.55314213000 | 1.31000000 | 2.5531890382965083 | 0.0 | TARGET | OFFLINE_NEW_DISCOVERY |
| 9650989 | Kepler-1905 b | 3.42392993000 | 1.44607340 | 9.084245050393813 | 165.3 | none | FALSE_POSITIVE |
| 9649706 | Kepler-1072 b | 1.56906650200 | 1.58000000 | 1.569123509271052 | 0.0 | TARGET | OFFLINE_NEW_DISCOVERY |
| 7269974 | Kepler-160 b / Kepler-160 c | 4.30939700000 | 1.71500000 | 4.5665468179551425 | 6.0 | none | OFFLINE_NEW_DISCOVERY |
| 5088400 | Kepler-1556 b | 8.82713457000 | 1.96000000 | 9.804789382573572 | 11.1 | none | FALSE_POSITIVE |
| 4736569 | Kepler-1042 b | 10.13202575000 | 2.19000000 | 10.131848857147352 | 0.0 | TARGET | OFFLINE_NEW_DISCOVERY |
| 11134879 | Kepler-570 b | 4.30166200300 | 2.45000000 | 4.301903477048307 | 0.0 | TARGET | OFFLINE_NEW_DISCOVERY |
| 5351250 | Kepler-150 b / Kepler-150 c / Kepler-150 d / Kepler-150 e | 12.56093000000 | 2.79000000 | 7.381876896739553 | 41.2 | planet | OFFLINE_NEW_DISCOVERY |
| 7447200 | Kepler-210 b / Kepler-210 c | 7.97251300000 | 3.62000000 | 7.972644185313584 | 0.0 | TARGET | OFFLINE_NEW_DISCOVERY |

## Target list (false / quiet)

| KIC | found P | snr | status |
|---|---|---|---|
| 7581037 | 12.73990263555481 | 10.637457429817669 | FALSE_POSITIVE |
| 7592369 | 10.266214514560785 | 7.810349202362609 | FALSE_POSITIVE |
| 7602179 | 6.978110882956879 | 12.593032831686044 | FALSE_POSITIVE |
| 7618003 | 5.879875109553024 | 11.048526437308324 | OFFLINE_NEW_DISCOVERY |
| 7658768 | 7.300303611835133 | 10.65318998930664 | FALSE_POSITIVE |
| 7668232 | 11.401969264035833 | 24.532231152523625 | FALSE_POSITIVE |
| 7677951 | 11.831072865026659 | 13.73379677710345 | FALSE_POSITIVE |
| 7691290 | 2.039460813294779 | 8.305160104924296 | OFFLINE_NEW_DISCOVERY |
| 7700503 | 8.413435299295775 | 9.857249462490392 | OFFLINE_NEW_DISCOVERY |
| 7733138 | 3.1410432395332877 | 8.280039735904582 | FALSE_POSITIVE |
| 7744036 | 11.92846445453269 | 11.682773483362336 | OFFLINE_NEW_DISCOVERY |
| 7757922 | 0.515106019058996 | 7.985197346212704 | FALSE_POSITIVE |

## Known signals per true host (NEA cumulative)

| KIC | P (d) | disp | name |
|---|---|---|---|
| 8478994 | 39.79220077 | CANDIDATE | Kepler-37 d |
| 8478994 | 21.30181863 | CANDIDATE | Kepler-37 c |
| 8478994 | 13.3669309 | CANDIDATE | Kepler-37 b |
| 8478994 | 51.20690303 | FALSE POSITIVE | Kepler-37 e |
| 8073705 | 10.60052083 | CANDIDATE | Kepler-1445 b |
| 4814502 | 4.94330596 | CANDIDATE | Kepler-1933 b |
| 9285568 | 2.55315031 | CANDIDATE | Kepler-1575 b |
| 9650989 | 3.42390809 | CANDIDATE | Kepler-1905 b |
| 9649706 | 1.569066598 | CANDIDATE | Kepler-1072 b |
| 7269974 | 13.69942266 | CANDIDATE | Kepler-160 c |
| 7269974 | 4.309381978 | CANDIDATE | Kepler-160 b |
| 5088400 | 8.82715783 | CANDIDATE | Kepler-1556 b |
| 4736569 | 10.13200555 | CANDIDATE | Kepler-1042 b |
| 4736569 | 7.07394528 | CANDIDATE |  |
| 11134879 | 4.30165843 | CANDIDATE | Kepler-570 b |
| 5351250 | 7.381980206 | CANDIDATE | Kepler-150 c |
| 5351250 | 12.56096618 | CANDIDATE | Kepler-150 d |
| 5351250 | 30.8261043 | CANDIDATE | Kepler-150 e |
| 5351250 | 3.42806306 | CANDIDATE | Kepler-150 b |
| 5351250 | 93.803766 | FALSE POSITIVE |  |
| 7447200 | 7.972510087 | CANDIDATE | Kepler-210 c |
| 7447200 | 2.453236139 | CANDIDATE | Kepler-210 b |

## Honesty caveats (REQUIRED READING)

1. **The quiet-star FPR is a proxy.** 'Quiet' = zero KOIs / zero
   confirmed planets in DR25. Undiscovered planets or undetected
   eclipsing binaries below the Kepler detection limit may still
   exist; a certification on such a star is flagged FP but is not
   proven to be noise.
2. **Truth = NEA period (5% window).** A certification at a period
   matching ANOTHER known planet of the same host is scored as
   'other known' (real signal, different planet), not as a false
   positive. A certification matching a cataloged FALSE POSITIVE
   (e.g. an EB) is scored as FP — correctly.
3. **Band aliases:** signals with true period >13.5 d can appear as
   in-band sub-harmonics (e.g. Kepler-37 d: P=39.8 d, 5th harmonic
   at 7.96 d certified). The alias resolver tests 2x/3x only.
4. Long cadence (29.4 min) + quarter gaps + real systematics make
   real-data recall lower than the controlled synthetic recall.
5. Periods were restricted to the pipeline search band (0.5-13.5 d);
   the 5% window is an acceptance rule, not a precision claim.
