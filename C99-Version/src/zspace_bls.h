/*
 * zspace_bls.h  ·  BLS Periodogram Search (C99)
 * Mirrors zspace_engine/detectors.py BLSDetector
 */
#ifndef ZSPACE_BLS_H
#define ZSPACE_BLS_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Configuration ─────────────────────────────────────────────────────────── */
typedef struct {
    double period_min_days;     /* e.g., 0.5 */
    double period_max_days;     /* e.g., 13.5 */
    double duration_min_hrs;    /* e.g., 0.5 */
    double duration_max_hrs;    /* e.g., 12.0 */
    int n_freq;                 /* frequency grid points (0 = auto from baseline) */
    int n_dur;                  /* duration grid points (0 = fixed Python grid) */
    double frequency_factor;    /* oversampling factor, e.g., 10.0 */
    int oversample;             /* BLS bin oversample (astropy: 10) */
    double fap_threshold;       /* FAP threshold for detection */
    double snr_threshold;       /* SNR threshold for detection */
    double period_prior_days;   /* optional period prior (0 = no prior) */
    int ladder_k;               /* harmonics for prior ladder, e.g., 20 */
    double ladder_min_rel_snr;  /* min relative SNR for prior peak, e.g., 0.05 */
    int store_grid;             /* 1 = keep full power_grid (ladder/prior), 0 = spectrum only */
} ZSBLSConfig;

/* ── BLS Result ────────────────────────────────────────────────────────────── */
typedef struct {
    double best_period_days;
    double best_power;
    double best_snr;
    double best_fap;
    double best_duration_hrs;
    double best_t0_days;
    double best_depth;
    int has_detection;
    /* FAP decomposition (detectors.py _finalize_result) */
    double fap_power;          /* power-spectrum exponential-tail FAP */
    double fap_snr;            /* matched-filter SNR trial-corrected FAP */
    double s_periodicity;      /* periodicity score S_P in [0,1] */
    double noise_floor;        /* median of power spectrum (FAP loc) */
    double noise_rms;          /* std of power spectrum */
    double fap_loc;            /* exponential-tail location */
    double fap_scale;          /* exponential-tail scale (MAD-normalised) */
    int    n_independent;      /* independent trial periods */
    /* Optional: full grid for ladder/prior selection */
    double *period_grid;        /* size n_periods */
    double *power_spectrum;     /* per-period max power, size n_periods */
    double *duration_spectrum;  /* per-period best duration (days), size n_periods */
    double *t0_spectrum;        /* per-period best t0 (days), size n_periods */
    double *depth_spectrum;     /* per-period best depth, size n_periods */
    double *power_grid;         /* size n_periods x n_dur (flattened) */
    int n_periods;
    int n_durations;
} ZSBLSResult;

/* ── Candidate from ladder (for prior search) ──────────────────────────────── */
typedef struct {
    double period_days;
    double power;
    double snr;
    double fap;
    double duration_hrs;
    double t0_days;
    double depth;
} ZSBLSCandidate;

/* ── Default config (mirrors zspace_engine/detectors.py defaults) ───────────── */
extern const ZSBLSConfig ZS_BLS_CONFIG_DEFAULT;

/* ── Core functions ────────────────────────────────────────────────────────── */

/* Run full BLS periodogram on light curve.
 * Returns 0 on success, -1 on error.
 * Result must be freed with zs_bls_result_free(). */
int zs_bls_search(const double *time, size_t n_time,
                  const double *flux, const double *flux_err,
                  const ZSBLSConfig *config,
                  ZSBLSResult *out);

/* Run BLS at a specific period (for prior ladder refinement).
 * Returns 0 on success. */
int zs_bls_run_at_period(const double *time, size_t n_time,
                         const double *flux, const double *flux_err,
                         double target_period_days, double duration_days,
                         ZSBLSCandidate *out);

/* Generate period-prior ladder candidates from a global BLS result.
 * Returns number of ladder candidates written to array (max ladder_k). */
int zs_bls_ladder_candidates(const ZSBLSResult *global,
                             double period_prior_days,
                             int ladder_k,
                             double ladder_min_rel_snr,
                             ZSBLSCandidate *ladder_out);

/* Candidate ladder: local-maxima of the power spectrum (mirrors
 * zspace_engine/detectors.py BLSDetector.top_candidates).
 * time/flux are needed to compute matched-filter SNR per candidate.
 * Returns number of candidates written (<= max_out, <= k). */
int zs_bls_top_candidates(const ZSBLSResult *global,
                          const double *time, size_t n_time,
                          const double *flux,
                          int k,
                          double min_relative_snr,
                          ZSBLSCandidate *out,
                          int max_out);

/* Free result resources */
void zs_bls_result_free(ZSBLSResult *res);

/* Free candidate resources */
void zs_bls_candidate_free(ZSBLSCandidate *cand);

#ifdef __cplusplus
}
#endif

#endif /* ZSPACE_BLS_H */