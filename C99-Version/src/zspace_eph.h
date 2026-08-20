/*
 * zspace_eph.h  ·  Harmonics / Alias Ephemeris Resolution (C99)
 * =============================================================
 * Mirrors zspace_engine/ephemeris.py EphemerisResolver:
 *   - fold & bin (BLSDetector.fold_and_bin), robust median/MAD noise
 *   - two-pass dip detection (3-MAD cores, wing expansion, relative floor)
 *   - own/2P/3P fold classification into FUNDAMENTAL / P_TRUE/2 / P_TRUE/3
 *   - over-harmonic branch (fold down at P/n, require single dip)
 */
#ifndef ZSPACE_EPH_H
#define ZSPACE_EPH_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double phase_center;   /* weighted phase of dip centre in [-0.5, 0.5] */
    double phase_lo;       /* phase of the first dip bin  */
    double phase_hi;       /* phase of the last dip bin   */
    double depth;          /* baseline - mean(binned flux inside group) >= 0 */
    double depth_snr;      /* depth * sqrt(n_bins_group) / per_bin_noise */
    double width_phase;    /* (phase_hi - phase_lo) + bin width */
} ZSEphDip;

typedef struct {
    double period;
    int n_bins;
    ZSEphDip *dips;        /* caller frees */
    int n_dips;
    double baseline;
    double per_bin_noise;
    double covered_fraction;
} ZSEphSignature;

typedef struct {
    double candidate_period;
    double physical_period;
    int multiple;          /* 1, 2 or 3 */
    char classifier[32];   /* FUNDAMENTAL | P_TRUE/2_ALIAS | P_TRUE/3_ALIAS | OVER_HARMONIC */
    double confidence;
    char evidence[512];
    char pattern[96];      /* {'own':n,'at_2p':n,'at_3p':n} compact string */
    int n_flags;
    char flags[4][40];
} ZSEphResult;

/* Fold and bin (mirrors BLSDetector.fold_and_bin). Bins with < 3 points are
   NaN (represent as `valid` flags). Returns 0 on success. */
int zs_eph_fold_and_bin(const double *time, size_t n_time,
                        const double *flux,
                        double period, double t0, int n_bins,
                        double *bin_phase,      /* n_bins */
                        double *bin_flux,       /* n_bins, NaN-free when valid */
                        unsigned char *valid);  /* n_bins */

/* Full dip signature (mirrors _fold_dip_signature). Caller frees sig->dips. */
int zs_eph_dip_signature(const double *time, size_t n_time,
                         const double *flux,
                         double period, double t0, int n_bins,
                         ZSEphSignature *sig);

/* Free a signature returned by zs_eph_dip_signature. */
void zs_eph_signature_free(ZSEphSignature *sig);

/* Resolve a candidate's ephemeris (mirrors EphemerisResolver.resolve). */
int zs_eph_resolve(const double *time, size_t n_time,
                   const double *flux,
                   double period_best, double t0,
                   double period_min, double period_max,
                   ZSEphResult *out);

#ifdef __cplusplus
}
#endif

#endif /* ZSPACE_EPH_H */