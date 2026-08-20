#include "zspace_bls.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    size_t n = 5000;
    double *time = malloc(n * sizeof(double));
    double *flux = malloc(n * sizeof(double));
    
    double period = 10.0;
    double depth = 0.01;
    double duration_hrs = 3.0;
    double t0 = 1.0;
    double cadence = 30.0 / (24.0 * 60.0);
    
    for (size_t i = 0; i < n; i++) {
        time[i] = i * cadence;
        double phase = fmod(time[i] - t0, period) / period;
        if (phase < 0) phase += 1.0;
        double dur_frac = duration_hrs / 24.0 / period;
        double in_transit = (phase < 0.5 * dur_frac) || (phase > 1.0 - 0.5 * dur_frac);
        flux[i] = 1.0 - (in_transit ? depth : 0.0);
        flux[i] += (rand() / (double)RAND_MAX - 0.5) * 0.001;
    }
    
    /* Manual boxcar test at true period */
    double freq = 1.0 / period;
    double dur_frac = duration_hrs / 24.0 / period;
    
    double sum_in = 0, sum_out = 0;
    int n_in = 0, n_out = 0;
    for (size_t i = 0; i < n; i++) {
        double ph = fmod(time[i] - 1.0, period) / period;
        if (ph < 0) ph += 1.0;
        int in_tr = (ph < 0.5 * dur_frac) || (ph > 1.0 - 0.5 * dur_frac);
        if (in_tr) { sum_in += flux[i]; n_in++; }
        else { sum_out += flux[i]; n_out++; }
    }
    
    double mean_in = sum_in / n_in;
    double mean_out = sum_out / n_out;
    double d = mean_out - mean_in;
    double snr = fabs(d) / 0.001 * sqrt(n_in);
    double power = snr * snr;
    
    printf("Manual boxcar at true period:\n");
    printf("  n_in=%d, n_out=%d\n", n_in, n_out);
    printf("  mean_in=%.6f, mean_out=%.6f\n", mean_in, mean_out);
    printf("  depth=%.6f, snr=%.2f, power=%.2f\n", d, snr, power);
    
    /* Test BLS at true period */
    ZSBLSConfig cfg = ZS_BLS_CONFIG_DEFAULT;
    cfg.period_min_days = period * 0.9;
    cfg.period_max_days = period * 1.1;
    cfg.n_freq = 100;
    cfg.n_dur = 5;
    
    ZSBLSResult res;
    ZSBLSCandidate cand;
    int ret = zs_bls_run_at_period(time, n, flux, NULL, period, duration_hrs/24.0, &cand);
    
    if (ret == 0) {
        printf("\nBLS run_at_period:\n");
        printf("  Period: %.4f days\n", cand.period_days);
        printf("  Power: %.4f\n", cand.power);
        printf("  SNR: %.2f\n", cand.snr);
        printf("  Depth: %.4f\n", cand.depth);
    } else {
        printf("BLS run_at_period failed\n");
    }
    
    free(time);
    free(flux);
    return 0;
}