#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(void) {
    size_t n_points = 3000;
    double baseline = 60.0, period_true = 3.695, depth_true = 0.008, dur_hrs_true = 2.5;
    double *time = malloc(n_points * sizeof(double));
    double *flux = malloc(n_points * sizeof(double));
    srand(42);
    for (size_t i = 0; i < n_points; i++) {
        time[i] = baseline * (double)i / (double)n_points;
        double phase = fmod((time[i] - 1.0) / period_true, 1.0);
        double dur_frac = (dur_hrs_true / 24.0) / period_true;
        int in_t = (phase < dur_frac) || (phase > 1.0 - dur_frac);
        double signal = in_t ? -depth_true : 0.0;
        double noise = ((double)rand() / RAND_MAX - 0.5) * 2.0 * 0.0015;
        flux[i] = 1.0 + signal + noise;
    }
    double t0b = 0.996875, pb = 3.695262, db = 5.0 / 24.0;
    double half = (db / pb) / 2.0;
    int n_in = 0, n_out = 0;
    double s_mean = 0.0;
    for (size_t i = 0; i < n_points; i++) {
        double ph = fmod((time[i] - t0b) / pb, 1.0);
        if (ph > 0.5) ph -= 1.0;
        double aph = fabs(ph);
        if (aph <= half) n_in++;
        else if (aph < 0.4) { n_out++; s_mean += flux[i]; }
    }
    s_mean /= n_out;
    double var = 0.0;
    double dmax = 0, dmin = 0;
    for (size_t i = 0; i < n_points; i++) {
        double ph = fmod((time[i] - t0b) / pb, 1.0);
        if (ph > 0.5) ph -= 1.0;
        double aph = fabs(ph);
        if (aph > half && aph < 0.4) {
            double d = flux[i] - s_mean;
            var += d * d;
            if (d > dmax) dmax = d;
            if (d < dmin) dmin = d;
        }
    }
    printf("n_in=%d n_out=%d s_mean=%.8f sigma=%.8f dmin=%.5f dmax=%.5f\n",
           n_in, n_out, s_mean, sqrt(var / (n_out - 1)), dmin, dmax);
    /* print offending points */
    for (size_t i = 0; i < n_points; i++) {
        double ph = fmod((time[i] - t0b) / pb, 1.0);
        if (ph > 0.5) ph -= 1.0;
        double aph = fabs(ph);
        if (aph > half && aph < 0.4 && fabs(flux[i] - s_mean) > 0.004) {
            printf("offender: t=%.3f ph=%.5f flux=%.6f\n", time[i], ph, flux[i]);
        }
    }
    return 0;
}