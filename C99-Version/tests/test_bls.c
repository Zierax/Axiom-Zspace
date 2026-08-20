#include "zspace_bls.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    /* Generate synthetic light curve: 10 day period, 1% depth, 3hr duration */
    size_t n = 5000;
    double *time = malloc(n * sizeof(double));
    double *flux = malloc(n * sizeof(double));
    
    double period = 10.0;
    double depth = 0.01;
    double duration_hrs = 3.0;
    double t0 = 1.0;
    double cadence = 30.0 / (24.0 * 60.0);  /* 30 min cadence in days */
    
    for (size_t i = 0; i < n; i++) {
        time[i] = i * cadence;
        double phase = fmod(time[i] - t0, period) / period;
        if (phase < 0) phase += 1.0;
        double in_transit = (phase < 0.5 * duration_hrs / 24.0 / period) ||
                            (phase > 1.0 - 0.5 * duration_hrs / 24.0 / period);
        flux[i] = 1.0 - (in_transit ? depth : 0.0);
        flux[i] += (rand() / (double)RAND_MAX - 0.5) * 0.001;  /* noise */
    }
    
    ZSBLSConfig cfg = ZS_BLS_CONFIG_DEFAULT;
    cfg.period_min_days = 5.0;
    cfg.period_max_days = 15.0;
    cfg.n_freq = 5000;
    cfg.n_dur = 10;
    
    ZSBLSResult res;
    int ret = zs_bls_search(time, n, flux, NULL, &cfg, &res);
    
    if (ret == 0) {
        printf("BLS Search Result:\n");
        printf("  Detection: %s\n", res.has_detection ? "YES" : "NO");
        printf("  Period: %.4f days (true: %.1f)\n", res.best_period_days, period);
        printf("  Power: %.4f\n", res.best_power);
        printf("  SNR: %.2f\n", res.best_snr);
        printf("  FAP: %.2e\n", res.best_fap);
        printf("  Duration: %.2f hrs (true: %.1f)\n", res.best_duration_hrs, duration_hrs);
        printf("  Depth: %.4f (true: %.3f)\n", res.best_depth, depth);
        
        double period_error = fabs(res.best_period_days - period) / period * 100;
        printf("  Period error: %.2f%%\n", period_error);
    } else {
        printf("BLS search failed\n");
    }
    
    zs_bls_result_free(&res);
    free(time);
    free(flux);
    return 0;
}