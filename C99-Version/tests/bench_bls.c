/*
 * bench_bls.c  ·  Benchmark C99 BLS periodogram
 * Generates synthetic light curve with injected transit, times the search.
 * Usage: bench_bls [n_points] [n_freq] [n_dur]
 */
#include "zspace_bls.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static double now_sec(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
#endif
}

int main(int argc, char **argv) {
    size_t n_points = (argc > 1) ? (size_t)atol(argv[1]) : 3000;
    int n_freq = (argc > 2) ? atoi(argv[2]) : 20000;
    int n_dur = (argc > 3) ? atoi(argv[3]) : 15;

    unsigned seed = 42;
    double *time = malloc(n_points * sizeof(double));
    double *flux = malloc(n_points * sizeof(double));
    if (!time || !flux) { fprintf(stderr, "OOM\n"); return 1; }

    /* Simulated 60-day baseline, 30-min cadence */
    double baseline = 60.0;
    double period_true = 3.695;
    double depth_true = 0.008;
    double dur_hrs_true = 2.5;

    srand(seed);
    for (size_t i = 0; i < n_points; i++) {
        time[i] = baseline * (double)i / (double)n_points;
        double phase = fmod((time[i] - 1.0) / period_true, 1.0);
        if (phase < 0.0) phase += 1.0;
        double dur_frac = (dur_hrs_true / 24.0) / period_true;
        int in_transit = (phase < dur_frac) || (phase > 1.0 - dur_frac);
        double signal = in_transit ? -depth_true : 0.0;
        double noise = ((double)rand() / RAND_MAX - 0.5) * 2.0 * 0.0015;
        flux[i] = 1.0 + signal + noise;
    }

    ZSBLSConfig cfg = ZS_BLS_CONFIG_DEFAULT;
    cfg.n_freq = n_freq;
    cfg.n_dur = n_dur;

    /* Warmup */
    ZSBLSResult res;
    zs_bls_search(time, n_points, flux, NULL, &cfg, &res);
    zs_bls_result_free(&res);

    /* Timed run */
    double t0 = now_sec();
    int rc = zs_bls_search(time, n_points, flux, NULL, &cfg, &res);
    double t1 = now_sec();

    if (rc != 0) { fprintf(stderr, "search failed\n"); return 1; }

    printf("n_points=%zu n_freq=%d n_dur=%d threads=%d\n",
           n_points, n_freq, n_dur,
#ifdef _OPENMP
           omp_get_max_threads()
#else
           1
#endif
    );
    printf("time_ms=%.2f\n", (t1 - t0) * 1000.0);
    printf("best_period=%.6f days\n", res.best_period_days);
    printf("best_power=%.6f\n", res.best_power);
    printf("best_snr=%.2f\n", res.best_snr);
    printf("best_duration=%.3f hrs\n", res.best_duration_hrs);
    printf("best_depth=%.6f\n", res.best_depth);
    printf("best_t0=%.6f\n", res.best_t0_days);
    printf("true_period=%.6f\n", period_true);

    /* Debug SNR internals */
    {
        double half_dur_phase = (res.best_duration_hrs / 24.0 / res.best_period_days) / 2.0;
        int n_in = 0, n_out = 0;
        double s_mean = 0.0;
        for (size_t i = 0; i < n_points; i++) {
            double ph = fmod((time[i] - res.best_t0_days) / res.best_period_days, 1.0);
            if (ph > 0.5) ph -= 1.0;
            double aph = fabs(ph);
            if (aph <= half_dur_phase) n_in++;
            else if (aph < 0.4) { n_out++; s_mean += flux[i]; }
        }
        s_mean /= n_out;
        double var = 0.0;
        for (size_t i = 0; i < n_points; i++) {
            double ph = fmod((time[i] - res.best_t0_days) / res.best_period_days, 1.0);
            if (ph > 0.5) ph -= 1.0;
            double aph = fabs(ph);
            if (aph > half_dur_phase && aph < 0.4) {
                double d = flux[i] - s_mean; var += d * d;
            }
        }
        printf("dbg: t0=%.6f n_in=%d n_out=%d sigma_oot=%.8f\n",
               res.best_t0_days, n_in, n_out, sqrt(var / (n_out - 1)));
    }

    zs_bls_result_free(&res);
    free(time);
    free(flux);
    return 0;
}