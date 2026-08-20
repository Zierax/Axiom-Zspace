#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zspace_ingestion.h"

int main(int argc, char **argv) {
    const size_t n = 91566;
    int window = (argc > 1) ? atoi(argv[1]) : 1621;
    int poly = 3;
    ZSLightCurve lc;
    memset(&lc, 0, sizeof(lc));
    lc.n = n;
    lc.time = malloc(n * sizeof(double));
    lc.flux = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) { lc.time[i] = (double)i * 0.00139; lc.flux[i] = 1.0; }
    int rc = zs_savgol_flatten(&lc, window, poly);
    printf("rc=%d flux[0]=%.10f flux[809]=%.10f flux[810]=%.10f flux[90755]=%.10f flux[91565]=%.10f\n",
           rc, lc.flux[0], lc.flux[809], lc.flux[810], lc.flux[90755], lc.flux[91565]);
    int bad = 0;
    for (size_t i = 0; i < n; i++)
        if (lc.flux[i] != 1.0) { bad++; if (bad < 6) printf("  bad[%zu]=%.10f\n", i, lc.flux[i]); }
    printf("bad=%d\n", bad);
    free(lc.time); free(lc.flux);
    return 0;
}