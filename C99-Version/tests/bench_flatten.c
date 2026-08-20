#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "zspace_ingestion.h"

static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

int main(int argc, char **argv) {
    long n = (argc > 1) ? atol(argv[1]) : 91566;
    int window = (argc > 2) ? atoi(argv[2]) : 1567;
    ZSLightCurve lc;
    memset(&lc, 0, sizeof(lc));
    lc.n = (size_t)n;
    lc.time = malloc(lc.n * sizeof(double));
    lc.flux = malloc(lc.n * sizeof(double));
    for (size_t i = 0; i < lc.n; i++) { lc.time[i] = i * 0.00139; lc.flux[i] = 1.0 + 0.001 * ((double)i / lc.n); }
    double t0 = now_ms();
    int rc = zs_savgol_flatten(&lc, window, 3);
    double t1 = now_ms();
    printf("n=%ld window=%d rc=%d time_ms=%.1f omp_threads=%d\n",
           n, window, rc, t1 - t0, omp_get_max_threads());
    free(lc.time); free(lc.flux);
    return 0;
}