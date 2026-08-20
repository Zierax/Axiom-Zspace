/*
 * zspace_bls.c  ·  BLS Periodogram Search (C99)
 * Mirrors astropy.timeseries.BoxLeastSquares (fast Cython method)
 * + zspace_engine/detectors.py BLSDetector (SNR/FAP/S_P finalization)
 *
 * Algorithm (per period):
 *   - bin data into n_bins phase bins (bin_duration = min_dur / oversample)
 *   - prefix sums over bins -> O(1) window sums
 *   - scan all window offsets x durations, maximize log-likelihood
 *     LL = 0.5 * ivar_in * (y_out - y_in)^2   [likelihood objective]
 *   - require y_out >= y_in (dip only)
 */
#include "zspace_bls.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const ZSBLSConfig ZS_BLS_CONFIG_DEFAULT = {
    .period_min_days = 0.5,
    .period_max_days = 13.5,
    .duration_min_hrs = 0.5,
    .duration_max_hrs = 12.0,
    .n_freq = 0,                 /* auto: derived from baseline */
    .n_dur = 0,                  /* auto: fixed Python-style grid */
    .frequency_factor = 20.0,
    .oversample = 10,
    .fap_threshold = 1e-4,
    .snr_threshold = 5.5,
    .period_prior_days = 0.0,
    .ladder_k = 20,
    .ladder_min_rel_snr = 0.05,
    .store_grid = 0,
};

/* Python reference duration grid (hours): [0.25 ... 8] (12 trimmed by period_min) */
static const double ZS_DUR_HRS_GRID[] = {
    0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0
};
#define ZS_N_DUR_GRID 14

/* ── Internal helpers ──────────────────────────────────────────────────────── */

static inline double wrap_into(double x, double period) {
    return x - period * floor(x / period);
}

static int cmp_double_asc(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* ── Full periodogram (astropy-mirrored) ───────────────────────────────────── */

int zs_bls_search(const double *time, size_t n_time,
                  const double *flux, const double *flux_err,
                  const ZSBLSConfig *config,
                  ZSBLSResult *out) {
    if (!time || !flux || n_time < 5 || !config || !out) return -1;
    (void)flux_err;

    memset(out, 0, sizeof(*out));

    double p_min = config->period_min_days;
    double p_max = config->period_max_days;
    if (p_min <= 0 || p_max <= p_min) return -1;

    /* ── Frequency grid (uniform in 1/P), matching detectors.py ────────────── */
    double baseline = time[n_time - 1] - time[0];
    if (baseline <= 0) return -1;
    double f_min = 1.0 / p_max;
    double f_max = 1.0 / p_min;
    double df = 1.0 / (baseline * config->frequency_factor);
    long n_freqs = (long)floor((f_max - f_min) / df);
    if (n_freqs < 2000) n_freqs = 2000;
    if (config->n_freq > 0) n_freqs = config->n_freq;

    /* ── Duration grid ─────────────────────────────────────────────────────── */
    double dur_grid[ZS_N_DUR_GRID];
    int n_dur = 0;
    double dur_limit_days = (p_min * 24.0 > 0.5) ? p_min : 0.5 / 24.0;
    dur_limit_days = fmax(p_min * 24.0, 0.5) / 24.0;
    if (config->n_dur > 0 && config->duration_min_hrs > 0) {
        /* caller-supplied linear grid */
        for (int d = 0; d < config->n_dur && d < ZS_N_DUR_GRID; d++) {
            double hrs = config->duration_min_hrs +
                (config->duration_max_hrs - config->duration_min_hrs) * d / (config->n_dur - 1);
            if (hrs / 24.0 >= dur_limit_days) break;
            dur_grid[d] = hrs / 24.0;
            n_dur++;
        }
    } else {
        /* Python reference fixed grid, trimmed by 12h < max(period_min*24, 0.5) */
        double lim_hrs = fmax(p_min * 24.0, 0.5);
        for (int d = 0; d < ZS_N_DUR_GRID; d++) {
            if (ZS_DUR_HRS_GRID[d] >= lim_hrs) break;
            dur_grid[d] = ZS_DUR_HRS_GRID[d] / 24.0;
            n_dur++;
        }
    }
    if (n_dur == 0) return -1;
    double min_duration = dur_grid[0];

    /* ── Allocate grids ────────────────────────────────────────────────────── */
    out->n_periods = (int)n_freqs;
    out->n_durations = n_dur;
    out->period_grid = malloc((size_t)n_freqs * sizeof(double));
    out->power_spectrum = malloc((size_t)n_freqs * sizeof(double));
    out->duration_spectrum = malloc((size_t)n_freqs * sizeof(double));
    out->t0_spectrum = malloc((size_t)n_freqs * sizeof(double));
    out->depth_spectrum = malloc((size_t)n_freqs * sizeof(double));
    out->power_grid = (config->store_grid)
        ? malloc((size_t)n_freqs * (size_t)n_dur * sizeof(double)) : NULL;
    if (!out->period_grid || !out->power_spectrum ||
        !out->duration_spectrum || !out->t0_spectrum || !out->depth_spectrum ||
        (config->store_grid && !out->power_grid)) {
        zs_bls_result_free(out);
        return -1;
    }
    for (long p = 0; p < n_freqs; p++) {
        out->power_spectrum[p] = 0.0;
        out->period_grid[p] = 0.0;
        out->duration_spectrum[p] = 0.0;
        out->t0_spectrum[p] = 0.0;
        out->depth_spectrum[p] = 0.0;
    }

    /* ── Pre-accumulate totals ─────────────────────────────────────────────── */
    double min_t = time[0];
    double sum_y = 0.0, sum_ivar = 0.0;
    for (size_t i = 0; i < n_time; i++) {
        if (time[i] < min_t) min_t = time[i];
        sum_y += flux[i];
        sum_ivar += 1.0;
    }

    /* ── Binning parameters (astropy: oversample=10) ───────────────────────── */
    int oversample = (config->oversample > 0) ? config->oversample : 10;
    double bin_duration = min_duration / (double)oversample;

    /* Per-thread scratch */
    int nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel
    { nthreads = omp_get_num_threads(); }
#endif
    size_t blocksize_max = (size_t)(ceil(p_max / bin_duration)) + oversample + 2;
    double *mean_y_0 = calloc((size_t)nthreads * blocksize_max, sizeof(double));
    double *mean_ivar_0 = calloc((size_t)nthreads * blocksize_max, sizeof(double));
    if (!mean_y_0 || !mean_ivar_0) {
        free(mean_y_0); free(mean_ivar_0);
        zs_bls_result_free(out);
        return -1;
    }

    /* ── Period scan ───────────────────────────────────────────────────────── */
    double best_power = -1.0;
    double best_period = 0.0, best_dur = 0.0, best_depth = 0.0, best_t0 = 0.0;

#pragma omp parallel if(!omp_in_parallel())
    {
        double local_best_power = -1.0;
        double local_period = 0.0, local_dur = 0.0, local_depth = 0.0, local_t0 = 0.0;

#ifdef _OPENMP
        int ithread = omp_get_thread_num();
#else
        int ithread = 0;
#endif
        double *mean_y = mean_y_0 + (size_t)ithread * blocksize_max;
        double *mean_ivar = mean_ivar_0 + (size_t)ithread * blocksize_max;

#pragma omp for schedule(guided, 32)
        for (long p = 0; p < n_freqs; p++) {
            /* periods low -> high: freq_grid = linspace(f_min, f_max); reverse */
            double freq = f_min + (f_max - f_min) * (double)p / (double)(n_freqs - 1);
            double period = 1.0 / freq;
            out->period_grid[p] = period;

            int n_bins = (int)(ceil(period / bin_duration)) + oversample;
            if ((size_t)n_bins + 2 > blocksize_max) n_bins = (int)blocksize_max - 2;
            if (n_bins < oversample + 2) continue;

            /* Bin the data */
            for (int n = 0; n < n_bins + 1; n++) { mean_y[n] = 0.0; mean_ivar[n] = 0.0; }
            for (size_t i = 0; i < n_time; i++) {
                double ph = wrap_into(time[i] - min_t, period);
                int ind = (int)(ph / bin_duration) + 1;
                if (ind > n_bins) ind = n_bins;
                mean_y[ind] += flux[i];
                mean_ivar[ind] += 1.0;
            }

            /* Wrap-around pad: copy first `oversample` bins to the end */
            for (int n = 1, idx = n_bins - oversample; n <= oversample; n++, idx++) {
                mean_y[idx] = mean_y[n];
                mean_ivar[idx] = mean_ivar[n];
            }

            /* Prefix sums */
            for (int n = 1; n <= n_bins; n++) {
                mean_y[n] += mean_y[n - 1];
                mean_ivar[n] += mean_ivar[n - 1];
            }

            double period_best_power = -1.0;
            double period_best_dur = 0.0, period_best_depth = 0.0, period_best_t0 = 0.0;

            for (int k = 0; k < n_dur; k++) {
                int dur = (int)lround(dur_grid[k] / bin_duration);
                if (dur < 1) dur = 1;
                int n_max = n_bins - dur;
                double k_best_ll = -1.0;
                double k_best_depth = 0.0, k_best_t0 = 0.0;
                for (int n = 0; n <= n_max; n++) {
                    double y_in = mean_y[n + dur] - mean_y[n];
                    double ivar_in = mean_ivar[n + dur] - mean_ivar[n];
                    double y_out = sum_y - y_in;
                    double ivar_out = sum_ivar - ivar_in;
                    if (ivar_in < DBL_EPSILON || ivar_out < DBL_EPSILON) continue;

                    y_in /= ivar_in;
                    y_out /= ivar_out;

                    /* likelihood objective (astropy default) */
                    double arg = y_out - y_in;
                    double log_like = 0.5 * ivar_in * arg * arg;

                    if (y_out >= y_in && log_like > k_best_ll) {
                        k_best_ll = log_like;
                        k_best_depth = y_out - y_in;
                        k_best_t0 = fmod(n * bin_duration +
                                         0.5 * dur * bin_duration + min_t, period);
                    }
                }
                if (config->store_grid)
                out->power_grid[(size_t)p * n_dur + k] =
                    (k_best_ll > 0) ? k_best_ll : 0.0;
                if (k_best_ll > period_best_power) {
                    period_best_power = k_best_ll;
                    period_best_dur = dur * bin_duration;
                    period_best_depth = k_best_depth;
                    period_best_t0 = k_best_t0;
                }
            }
            out->power_spectrum[p] = (period_best_power > 0) ? period_best_power : 0.0;
            out->duration_spectrum[p] = period_best_dur;
            out->t0_spectrum[p] = period_best_t0;
            out->depth_spectrum[p] = period_best_depth;

            if (period_best_power > local_best_power) {
                local_best_power = period_best_power;
                local_period = period;
                local_dur = period_best_dur;
                local_depth = period_best_depth;
                local_t0 = period_best_t0;
            }
        }

#pragma omp critical
        {
            if (local_best_power > best_power) {
                best_power = local_best_power;
                best_period = local_period;
                best_dur = local_dur;
                best_depth = local_depth;
                best_t0 = local_t0;
            }
        }
    }

    free(mean_y_0);
    free(mean_ivar_0);

    if (best_power < 0) {
        zs_bls_result_free(out);
        return -1;
    }

    out->best_period_days = best_period;
    out->best_power = best_power;
    out->best_duration_hrs = best_dur * 24.0;
    out->best_depth = best_depth;
    out->best_t0_days = best_t0;

    /* Matched-filter SNR (detectors.py _finalize_result) */
    double sigma_oot = 0.0;
    int n_in = 0, n_out = 0;
    {
        double half_dur_phase = (best_dur / best_period) / 2.0;
        double s_mean = 0.0;
        for (size_t i = 0; i < n_time; i++) {
            double ph = wrap_into(time[i] - best_t0, best_period) / best_period;
            if (ph > 0.5) ph -= 1.0;
            double aph = fabs(ph);
            if (aph <= half_dur_phase) n_in++;
            else if (aph < 0.4) { n_out++; s_mean += flux[i]; }
        }
        if (n_out > 0) {
            s_mean /= n_out;
            double var = 0.0;
            double half_dur_phase2 = (best_dur / best_period) / 2.0;
            for (size_t i = 0; i < n_time; i++) {
                double ph = wrap_into(time[i] - best_t0, best_period) / best_period;
                if (ph > 0.5) ph -= 1.0;
                double aph = fabs(ph);
                if (aph <= half_dur_phase2) { /* in */ }
                else if (aph < 0.4) { double d = flux[i] - s_mean; var += d * d; }
            }
            sigma_oot = sqrt(var / (n_out > 1 ? n_out - 1 : 1.0));
        }
    }
    if (n_in >= 3 && n_out >= 10 && sigma_oot > 0) {
        out->best_snr = fabs(best_depth) / (sigma_oot * sqrt(1.0 / n_in + 1.0 / n_out));
    } else {
        out->best_snr = 0.0;
    }
    out->has_detection = (out->best_snr > config->snr_threshold);

    /* ── FAP (detectors.py FAPValidator.from_power_spectrum + SNR form) ────── */
    double baseline_d = time[n_time - 1] - time[0];
    double delta_f = 1.0 / fmax(baseline_d, 1.0);
    int n_indep = (int)lround((1.0 / p_min - 1.0 / p_max) / delta_f);
    if (n_indep < 1) n_indep = 1;
    out->n_independent = n_indep;

    double fap_power = 1.0, fap_snr = 1.0;
    if (n_freqs < 10) {
        fap_power = 1.0;
    } else {
        /* noise floor / RMS of the power spectrum */
        double *ps_copy = malloc((size_t)n_freqs * sizeof(double));
        if (ps_copy) {
            memcpy(ps_copy, out->power_spectrum, (size_t)n_freqs * sizeof(double));
            qsort(ps_copy, (size_t)n_freqs, sizeof(double), cmp_double_asc);
            double noise_floor = ps_copy[n_freqs / 2];
            double mean = 0.0;
            for (long i = 0; i < n_freqs; i++) mean += ps_copy[i];
            mean /= n_freqs;
            double var = 0.0;
            for (long i = 0; i < n_freqs; i++) {
                double d = ps_copy[i] - mean;
                var += d * d;
            }
            double noise_rms = sqrt(var / n_freqs);

            out->noise_floor = noise_floor;
            out->noise_rms = noise_rms;

            if (noise_rms < 1e-30) {
                fap_power = 0.0;
            } else {
                long cut = (long)((double)n_freqs * 0.98);
                long min_cut = (long)n_freqs - 2;
                if (cut < min_cut) cut = min_cut;
                if (cut < 1) cut = 1;
                /* loc = percentile(noise_samples, 50) = median of first cut */
                double loc = ps_copy[cut / 2];
                /* MAD of noise_samples (first cut entries) */
                double *abs_dev = malloc((size_t)cut * sizeof(double));
                if (!abs_dev) {
                    free(ps_copy);
                    fap_power = 1.0;
                } else {
                    for (long i = 0; i < cut; i++) abs_dev[i] = fabs(ps_copy[i] - loc);
                    /* median of abs_dev */
                    double *dev_copy = malloc((size_t)cut * sizeof(double));
                    if (!dev_copy) {
                        free(abs_dev); free(ps_copy);
                        fap_power = 1.0;
                    } else {
                        memcpy(dev_copy, abs_dev, (size_t)cut * sizeof(double));
                        qsort(dev_copy, (size_t)cut, sizeof(double), cmp_double_asc);
                        double mad = 1.4826 * dev_copy[cut / 2];
                        double scale = fmax(mad, noise_rms * 1e-3);
                        out->fap_loc = loc;
                        out->fap_scale = scale;
                        double p_single = exp(-fmax(best_power - loc, 0.0) / scale);
                        if (p_single < 1e-300) p_single = 1e-300;
                        if (p_single > 1.0) p_single = 1.0;
                        fap_power = 1.0 - pow(1.0 - p_single, (double)n_indep);
                        if (fap_power > 1.0) fap_power = 1.0;
                        if (fap_power < 0.0) fap_power = 0.0;
                        free(dev_copy);
                    }
                    free(abs_dev);
                }
            }
            free(ps_copy);
        }
    }

    /* FAP — matched-filter SNR trial-corrected (SPOC-style) */
    if (out->best_snr > 0) {
        double p_noise = 0.5 * erfc(out->best_snr / sqrt(2.0));
        fap_snr = 1.0 - pow(1.0 - p_noise, (double)n_indep);
        if (fap_snr > 1.0) fap_snr = 1.0;
        if (fap_snr < 0.0) fap_snr = 0.0;
    } else {
        fap_snr = 1.0;
    }

    out->fap_power = fap_power;
    out->fap_snr = fap_snr;
    out->best_fap = fmin(fap_power, fap_snr);

    /* Periodicity score S_P */
    if (out->best_fap >= config->fap_threshold ||
        out->best_snr <= config->snr_threshold) {
        out->s_periodicity = 0.0;
    } else {
        double snr_ref = 50.0;
        out->s_periodicity = fmin(1.0, (out->best_snr - config->snr_threshold) /
                                        (snr_ref - config->snr_threshold));
    }

    return 0;
}

/* ── Run BLS at a specific period (prior ladder refinement) ───────────────── */

int zs_bls_run_at_period(const double *time, size_t n_time,
                         const double *flux, const double *flux_err,
                         double target_period_days, double duration_days,
                         ZSBLSCandidate *out) {
    if (!time || !flux || !out || target_period_days <= 0) return -1;
    (void)flux_err;

    memset(out, 0, sizeof(*out));

    double min_t = time[0];
    double sum_y = 0.0, sum_ivar = 0.0;
    for (size_t i = 0; i < n_time; i++) {
        if (time[i] < min_t) min_t = time[i];
        sum_y += flux[i];
        sum_ivar += 1.0;
    }

    int oversample = 10;
    double bin_duration = (duration_days / 2.0) / (double)oversample;
    int n_bins = (int)(ceil(target_period_days / bin_duration)) + oversample;

    double *mean_y = calloc((size_t)n_bins + 2, sizeof(double));
    double *mean_ivar = calloc((size_t)n_bins + 2, sizeof(double));
    if (!mean_y || !mean_ivar) { free(mean_y); free(mean_ivar); return -1; }

    for (size_t i = 0; i < n_time; i++) {
        double ph = wrap_into(time[i] - min_t, target_period_days);
        int ind = (int)(ph / bin_duration) + 1;
        if (ind > n_bins) ind = n_bins;
        mean_y[ind] += flux[i];
        mean_ivar[ind] += 1.0;
    }
    for (int n = 1, idx = n_bins - oversample; n <= oversample; n++, idx++) {
        mean_y[idx] = mean_y[n];
        mean_ivar[idx] = mean_ivar[n];
    }
    for (int n = 1; n <= n_bins; n++) {
        mean_y[n] += mean_y[n - 1];
        mean_ivar[n] += mean_ivar[n - 1];
    }

    int dur = (int)lround(duration_days / bin_duration);
    if (dur < 1) dur = 1;
    int n_max = n_bins - dur;
    double best_power = -1.0, best_depth = 0.0, best_t0 = 0.0;

    for (int n = 0; n <= n_max; n++) {
        double y_in = mean_y[n + dur] - mean_y[n];
        double ivar_in = mean_ivar[n + dur] - mean_ivar[n];
        double y_out = sum_y - y_in;
        double ivar_out = sum_ivar - ivar_in;
        if (ivar_in < DBL_EPSILON || ivar_out < DBL_EPSILON) continue;
        y_in /= ivar_in;
        y_out /= ivar_out;
        double arg = y_out - y_in;
        double log_like = 0.5 * ivar_in * arg * arg;
        if (y_out >= y_in && log_like > best_power) {
            best_power = log_like;
            best_depth = y_out - y_in;
            best_t0 = fmod(n * bin_duration + 0.5 * dur * bin_duration + min_t,
                           target_period_days);
        }
    }

    free(mean_y);
    free(mean_ivar);

    out->period_days = target_period_days;
    out->power = best_power > 0 ? best_power : 0.0;
    out->snr = 0.0;
    out->fap = exp(-2.0 * out->power);
    out->duration_hrs = duration_days * 24.0;
    out->t0_days = best_t0;
    out->depth = best_depth;
    return 0;
}

/* ── Ladder candidates from prior ─────────────────────────────────────────── */

int zs_bls_ladder_candidates(const ZSBLSResult *global,
                             double period_prior_days,
                             int ladder_k,
                             double ladder_min_rel_snr,
                             ZSBLSCandidate *ladder_out) {
    if (!global || !ladder_out || period_prior_days <= 0) return 0;
    (void)ladder_min_rel_snr;

    int n = 0;
    double factors[] = {1.0, 2.0, 0.5, 3.0, 1.0 / 3.0, 4.0, 0.25};
    for (size_t f = 0; f < sizeof(factors) / sizeof(factors[0]); f++) {
        double target = period_prior_days * factors[f];
        double tol = 0.05 * period_prior_days;
        double best_period = 0.0, best_power = -1.0;
        for (int i = 0; i < global->n_periods; i++) {
            double p = global->period_grid[i];
            if (fabs(p - target) <= tol) {
                double pw = 0.0;
                for (int d = 0; d < global->n_durations; d++) {
                    double v = global->power_grid[(size_t)i * global->n_durations + d];
                    if (v > pw) pw = v;
                }
                if (pw > best_power) { best_power = pw; best_period = p; }
            }
        }
        if (best_power > 0) {
            ladder_out[n].period_days = best_period;
            ladder_out[n].power = best_power;
            ladder_out[n].duration_hrs = 0.0;
            ladder_out[n].t0_days = global->best_t0_days;
            ladder_out[n].depth = 0.0;
            ladder_out[n].snr = sqrt(best_power);
            ladder_out[n].fap = exp(-2.0 * best_power);
            n++;
        }
    }
    return n;
}

/* ── Top-candidates ladder (local maxima, mirrors detectors.py) ──────────── */

static double zs_bls_matched_snr(const double *time, size_t n_time,
                                 const double *flux,
                                 double period, double t0, double dur_days,
                                 double depth) {
    if (!time || !flux || n_time < 5 || period <= 0 || dur_days <= 0) return 0.0;
    double half = dur_days / 2.0;
    /* in-transit vs out-of-transit (phase <0.4) like the main finalizer */
    int n_in = 0, n_out = 0;
    double sum_out = 0.0;
    for (size_t i = 0; i < n_time; i++) {
        double ph = wrap_into(time[i] - t0, period) / period;
        if (ph > 0.5) ph -= 1.0;
        double aph = fabs(ph);
        double thr = half / period;
        if (aph <= thr) n_in++;
        else if (aph < 0.4) { n_out++; sum_out += flux[i]; }
    }
    if (n_in < 3 || n_out < 10) return 0.0;
    double mean_out = sum_out / n_out;
    double var = 0.0;
    for (size_t i = 0; i < n_time; i++) {
        double ph = wrap_into(time[i] - t0, period) / period;
        if (ph > 0.5) ph -= 1.0;
        double aph = fabs(ph);
        double thr = half / period;
        if (aph <= thr) continue;
        else if (aph < 0.4) { double d = flux[i] - mean_out; var += d * d; }
    }
    double sigma = sqrt(var / (n_out > 1 ? n_out - 1 : 1));
    if (!(sigma > 0) || !isfinite(sigma)) return 0.0;
    return fabs(depth) / (sigma * sqrt(1.0 / n_in + 1.0 / n_out));
}

int zs_bls_top_candidates(const ZSBLSResult *global,
                          const double *time, size_t n_time,
                          const double *flux,
                          int k,
                          double min_relative_snr,
                          ZSBLSCandidate *out,
                          int max_out) {
    if (!global || !out || k <= 0 || max_out <= 0) return 0;
    int n_periods = global->n_periods;
    if (n_periods < 3) return 0;
    const double *ps = global->power_spectrum;
    const double *pg = global->period_grid;
    const double *dg = global->duration_spectrum;
    const double *tg = global->t0_spectrum;
    const double *depthg = global->depth_spectrum;
    if (!ps || !pg || !dg || !tg || !depthg) return 0;

    /* collect strict local maxima */
    int *peaks = (int *)malloc((size_t)n_periods * sizeof(int));
    if (!peaks) return 0;
    int n_peaks = 0;
    for (int i = 1; i < n_periods - 1; i++) {
        if (ps[i] > ps[i - 1] && ps[i] > ps[i + 1]) peaks[n_peaks++] = i;
    }
    if (n_peaks == 0) {
        int idx = 0;
        double mx = ps[0];
        for (int i = 1; i < n_periods; i++) if (ps[i] > mx) { mx = ps[i]; idx = i; }
        peaks[n_peaks++] = idx;
    }
    /* sort peaks by power desc (simple insertion, n_peaks ~ few k) */
    for (int i = 1; i < n_peaks; i++) {
        int key = peaks[i];
        double key_pow = ps[key];
        int j = i - 1;
        while (j >= 0 && ps[peaks[j]] < key_pow) {
            peaks[j + 1] = peaks[j];
            j--;
        }
        peaks[j + 1] = key;
    }
    /* global max power for relative threshold */
    double global_max = ps[peaks[0]];
    for (int i = 1; i < n_peaks; i++) if (ps[peaks[i]] > global_max) global_max = ps[peaks[i]];

    int n_chosen = 0;
    int limit = k < max_out ? k : max_out;
    for (int pi = 0; pi < n_peaks; pi++) {
        int idx = peaks[pi];
        double P = pg[idx];
        double tau = dg[idx];
        double pw = ps[idx];
        if (tau / fmax(P, 1e-9) > 0.15) continue;
        if (P < 0.3) continue;
        int alias = 0;
        for (int c = 0; c < n_chosen; c++) {
            double ratio = P / fmax(out[c].period_days, 1e-9);
            if (fabs(log(ratio)) < 0.10) { alias = 1; break; }
        }
        if (alias) continue;
        if (n_chosen >= limit) break;
        if (pw < min_relative_snr * global_max) continue;
        /* compute per-candidate SNR and FAP */
        double t0 = tg[idx];
        double depth = depthg[idx];
        double dur_hrs = tau * 24.0;
        double snr = 0.0;
        if (time && flux && n_time > 0) {
            snr = zs_bls_matched_snr(time, n_time, flux, P, t0, tau, depth);
        }
        double fap = 1.0;
        if (global->n_independent > 0 && global->fap_scale > 1e-12) {
            double p_single = exp(-fmax(pw - global->fap_loc, 0.0) / global->fap_scale);
            if (p_single < 1e-300) p_single = 1e-300;
            if (p_single > 1.0) p_single = 1.0;
            fap = 1.0 - pow(1.0 - p_single, (double)global->n_independent);
            if (fap > 1.0) fap = 1.0;
            if (fap < 0.0) fap = 0.0;
        } else {
            fap = exp(-2.0 * pw);
        }
        out[n_chosen].period_days = P;
        out[n_chosen].power = pw;
        out[n_chosen].snr = snr;
        out[n_chosen].fap = fap;
        out[n_chosen].duration_hrs = dur_hrs;
        out[n_chosen].t0_days = t0;
        out[n_chosen].depth = depth;
        n_chosen++;
        if (n_chosen >= limit) break;
    }
    free(peaks);
    return n_chosen;
}

/* ── Free functions ───────────────────────────────────────────────────────── */

void zs_bls_result_free(ZSBLSResult *res) {
    if (!res) return;
    free(res->period_grid);
    free(res->power_spectrum);
    free(res->duration_spectrum);
    free(res->t0_spectrum);
    free(res->depth_spectrum);
    free(res->power_grid);
    res->period_grid = NULL;
    res->power_spectrum = NULL;
    res->duration_spectrum = NULL;
    res->t0_spectrum = NULL;
    res->depth_spectrum = NULL;
    res->power_grid = NULL;
    res->n_periods = 0;
    res->n_durations = 0;
}

void zs_bls_candidate_free(ZSBLSCandidate *cand) {
    (void)cand;
}