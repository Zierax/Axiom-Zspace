#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "zspace_ingestion.h"

/* minimal re-declaration to inspect internals */
typedef struct {
    double *V;
    double *vnorm;
    double rdiag[11];
    double xmid;
    double xscale;
    size_t W;
    int p;
} SGQR_t;

static void qr_dump(const double *V, size_t W, int p, const double *rdiag, const double *vnorm) {
    size_t n = (size_t)p + 1;
    printf("rdiag: %.6e %.6e %.6e\n", rdiag[0], rdiag[1], rdiag[2]);
    printf("vnorm: %.6e %.6e %.6e\n", vnorm[0], vnorm[1], vnorm[2]);
    for (size_t k = 0; k < n; k++)
        printf("  Rrow[%zu]: %.6e %.6e %.6e %.6e\n", k,
               V[k*n+0], V[k*n+1], V[k*n+2], V[k*n+3]);
}

int main(void) {
    const size_t n = 10000;
    int window = 1621, poly = 3;
    ZSLightCurve lc;
    memset(&lc, 0, sizeof(lc));
    lc.n = n;
    lc.time = malloc(n * sizeof(double));
    lc.flux = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) { lc.time[i] = (double)i; lc.flux[i] = 1.0; }
    int rc = zs_savgol_flatten(&lc, window, poly);
    printf("rc=%d flux[0]=%.6e flux[1]=%.6e flux[809]=%.6e flux[810]=%.6e flux[9000]=%.6e flux[9999]=%.6e\n",
           rc, lc.flux[0], lc.flux[1], lc.flux[809], lc.flux[810], lc.flux[9000], lc.flux[9999]);
    free(lc.time); free(lc.flux);
    return 0;
}