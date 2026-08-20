/*
 * zspace_ingestion.c  ·  Light curve ingestion implementation
 * Requires cfitsio (https://heasarc.gsfc.nasa.gov/fitsio/)
 */
#include "zspace_ingestion.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef HAVE_CFITSIO
#include <fitsio.h>
#else
/* Stub definitions when cfitsio not available */
typedef int fitsfile;
#define READONLY 0
#endif

/* ── Internal helpers ──────────────────────────────────────────────────────── */

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

static double median(double *arr, size_t n) {
    if (n == 0) return NAN;
    qsort(arr, n, sizeof(double), cmp_double);
    if (n % 2 == 0) return (arr[n/2 - 1] + arr[n/2]) * 0.5;
    return arr[n/2];
}

static double mad(double *arr, size_t n, double med) {
    if (n == 0) return NAN;
    double *devs = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) devs[i] = fabs(arr[i] - med);
    double m = median(devs, n);
    free(devs);
    return m * 1.4826;  /* scale to std for normal */
}

/* Solve (p+1)x(p+1) linear system A x = b with partial pivoting. */
static int solve_small(double *A, double *b, int n) {
    for (int col = 0; col < n; col++) {
        int piv = col;
        double best = fabs(A[col * n + col]);
        for (int r = col + 1; r < n; r++) {
            double v = fabs(A[r * n + col]);
            if (v > best) { best = v; piv = r; }
        }
        if (best == 0.0) return -1;
        if (piv != col) {
            for (int c = 0; c < n; c++) {
                double tmp = A[col * n + c]; A[col * n + c] = A[piv * n + c]; A[piv * n + c] = tmp;
            }
            double tmp = b[col]; b[col] = b[piv]; b[piv] = tmp;
        }
        double d = A[col * n + col];
        for (int c = 0; c < n; c++) A[col * n + c] /= d;
        b[col] /= d;
        for (int r = 0; r < n; r++) {
            if (r == col) continue;
            double f = A[r * n + col];
            if (f == 0.0) continue;
            for (int c = 0; c < n; c++) A[r * n + c] -= f * A[col * n + c];
            b[r] -= f * b[col];
        }
    }
    return 0;
}

/* SG filter coefficients via least squares on the Vandermonde matrix
 * (matches scipy.signal._savitzky_golay.sgolay_coeffs: first row of lstsq(V, I)). */
static int sg_coeffs(int window_len, int polyorder, double **coeffs_out) {
    int p = polyorder, h = window_len / 2;
    if (p < 0 || p >= window_len) return -1;
    if (p > 10) return -1;
    double *V = malloc((size_t)window_len * (p + 1) * sizeof(double));
    double *A = malloc((size_t)(p + 1) * (p + 1) * sizeof(double));
    double *b = malloc((size_t)(p + 1) * sizeof(double));
    double *z = malloc((size_t)(p + 1) * sizeof(double));
    if (!V || !A || !b || !z) { free(V); free(A); free(b); free(z); return -1; }
    for (int i = 0; i < window_len; i++) {
        double x = (double)(i - h);
        double xp = 1.0;
        for (int j = 0; j <= p; j++) {
            V[i * (p + 1) + j] = xp;
            xp *= x;
        }
    }
    for (int j = 0; j <= p; j++)
        for (int k = 0; k <= p; k++) {
            double s = 0.0;
            for (int i = 0; i < window_len; i++) s += V[i * (p + 1) + j] * V[i * (p + 1) + k];
            A[j * (p + 1) + k] = s;
        }
    for (int j = 0; j <= p; j++) b[j] = (j == 0) ? 1.0 : 0.0;
    if (solve_small(A, b, p + 1) != 0) {
        free(V); free(A); free(b); free(z); return -1;
    }
    double *c = malloc((size_t)window_len * sizeof(double));
    if (!c) { free(V); free(A); free(b); free(z); return -1; }
    for (int i = 0; i < window_len; i++) {
        double s = 0.0;
        for (int j = 0; j <= p; j++) s += V[i * (p + 1) + j] * b[j];
        c[i] = s;
    }
    free(V); free(A); free(b); free(z);
    *coeffs_out = c;
    return 0;
}

/* ── Light curve allocation/free ───────────────────────────────────────────── */

void zs_lightcurve_free(ZSLightCurve *lc) {
    if (!lc) return;
    free(lc->time);
    free(lc->flux);
    free(lc->flux_err);
    free(lc->quality);
    lc->time = lc->flux = lc->flux_err = NULL;
    lc->quality = NULL;
    lc->n = 0;
}

static int lc_ensure_capacity(ZSLightCurve *lc, size_t needed) {
    if (lc->n + needed <= (lc->time ? (size_t)(-1) : 0)) return 0; /* already allocated */
    size_t new_cap = lc->n ? lc->n * 2 : 4096;
    while (new_cap < needed) new_cap *= 2;
    double *nt = realloc(lc->time, new_cap * sizeof(double));
    double *nf = realloc(lc->flux, new_cap * sizeof(double));
    double *ne = realloc(lc->flux_err, new_cap * sizeof(double));
    uint8_t *nq = realloc(lc->quality, new_cap * sizeof(uint8_t));
    if (!nt || !nf || !ne || !nq) {
        free(nt); free(nf); free(ne); free(nq);
        return -1;
    }
    lc->time = nt; lc->flux = nf; lc->flux_err = ne; lc->quality = nq;
    return 0;
}

/* ── FITS reading ──────────────────────────────────────────────────────────── */

int zs_read_fits(const char *path, ZSLightCurve *out, char *errbuf, size_t errbuf_sz) {
#ifdef HAVE_CFITSIO
    fitsfile *fptr = NULL;
    int status = 0;
    int hdutype;
    long nrows = 0;
    int col_time, col_flux, col_flux_err, col_qual;
    double *time_arr = NULL, *flux_arr = NULL, *flux_err_arr = NULL;
    uint8_t *qual_arr = NULL;

    if (!out) { if (errbuf) snprintf(errbuf, errbuf_sz, "null output"); return -1; }

    if (fits_open_file(&fptr, path, READONLY, &status)) {
        if (errbuf) snprintf(errbuf, errbuf_sz, "fits_open_file failed: %d", status);
        return -1;
    }

    /* Move to first BINTABLE HDU (usually HDU 1 for light curves) */
    if (fits_movabs_hdu(fptr, 2, &hdutype, &status)) {
        if (fits_movabs_hdu(fptr, 1, &hdutype, &status)) {
            if (errbuf) snprintf(errbuf, errbuf_sz, "no BINTABLE HDU found");
            fits_close_file(fptr, &status);
            return -1;
        }
    }

    if (fits_get_num_rows(fptr, &nrows, &status) || nrows == 0) {
        if (errbuf) snprintf(errbuf, errbuf_sz, "no rows in FITS table");
        fits_close_file(fptr, &status);
        return -1;
    }

    /* Find column indices (TESS/Kepler standard names) */
    fits_get_colnum(fptr, CASEINSEN, "TIME", &col_time, &status);
    fits_get_colnum(fptr, CASEINSEN, "SAP_FLUX", &col_flux, &status);
    if (status) { fits_get_colnum(fptr, CASEINSEN, "PDCSAP_FLUX", &col_flux, &status); }
    if (status) { fits_get_colnum(fptr, CASEINSEN, "FLUX", &col_flux, &status); }
    fits_get_colnum(fptr, CASEINSEN, "SAP_FLUX_ERR", &col_flux_err, &status);
    if (status) { fits_get_colnum(fptr, CASEINSEN, "PDCSAP_FLUX_ERR", &col_flux_err, &status); }
    if (status) { fits_get_colnum(fptr, CASEINSEN, "FLUX_ERR", &col_flux_err, &status); }
    fits_get_colnum(fptr, CASEINSEN, "QUALITY", &col_qual, &status);

    if (status || col_time <= 0 || col_flux <= 0) {
        if (errbuf) snprintf(errbuf, errbuf_sz, "required columns not found");
        fits_close_file(fptr, &status);
        return -1;
    }

    time_arr = calloc(nrows, sizeof(double));
    flux_arr = calloc(nrows, sizeof(double));
    flux_err_arr = calloc(nrows, sizeof(double));
    qual_arr = calloc(nrows, sizeof(uint8_t));
    if (!time_arr || !flux_arr || !flux_err_arr || !qual_arr) {
        free(time_arr); free(flux_arr); free(flux_err_arr); free(qual_arr);
        if (errbuf) snprintf(errbuf, errbuf_sz, "memory allocation failed");
        fits_close_file(fptr, &status);
        return -1;
    }

    /* Read columns */
    long firstrow = 1, firstelem = 1;
    fits_read_col(fptr, TDOUBLE, col_time, firstrow, firstelem, nrows, NULL, time_arr, NULL, &status);
    fits_read_col(fptr, TDOUBLE, col_flux, firstrow, firstelem, nrows, NULL, flux_arr, NULL, &status);
    if (col_flux_err > 0) {
        fits_read_col(fptr, TDOUBLE, col_flux_err, firstrow, firstelem, nrows, NULL, flux_err_arr, NULL, &status);
    }
    if (col_qual > 0) {
        long *qual_long = malloc(nrows * sizeof(long));
        fits_read_col(fptr, TLONG, col_qual, firstrow, firstelem, nrows, NULL, qual_long, NULL, &status);
        for (long i = 0; i < nrows; i++) qual_arr[i] = (uint8_t)qual_long[i];
        free(qual_long);
    }

    fits_close_file(fptr, &status);

    if (status) {
        free(time_arr); free(flux_arr); free(flux_err_arr); free(qual_arr);
        if (errbuf) snprintf(errbuf, errbuf_sz, "column read failed: %d", status);
        return -1;
    }

    /* Populate output */
    out->time = time_arr;
    out->flux = flux_arr;
    out->flux_err = flux_err_arr;
    out->quality = qual_arr;
    out->n = nrows;

    /* Compute metadata */
    out->t_start = time_arr[0];
    out->t_stop = time_arr[nrows - 1];
    if (nrows > 1) {
        double sum_dt = 0;
        for (long i = 1; i < nrows; i++) sum_dt += time_arr[i] - time_arr[i-1];
        out->cadence_days = sum_dt / (nrows - 1);
    } else {
        out->cadence_days = 0;
    }

    return 0;
#else
    (void)path; (void)out;
    if (errbuf) snprintf(errbuf, errbuf_sz, "cfitsio not available - compile with cfitsio");
    return -1;
#endif
}

/* ── Quality masking ───────────────────────────────────────────────────────── */

size_t zs_mask_quality(ZSLightCurve *lc, uint32_t quality_mask) {
    if (!lc || lc->n == 0) return 0;
    size_t write_idx = 0;
    for (size_t i = 0; i < lc->n; i++) {
        if ((lc->quality[i] & quality_mask) == 0) {
            if (write_idx != i) {
                lc->time[write_idx] = lc->time[i];
                lc->flux[write_idx] = lc->flux[i];
                lc->flux_err[write_idx] = lc->flux_err[i];
                lc->quality[write_idx] = lc->quality[i];
            }
            write_idx++;
        }
    }
    lc->n = write_idx;
    return write_idx;
}

/* ── Sigma clipping ────────────────────────────────────────────────────────── */

size_t zs_sigma_clip(ZSLightCurve *lc, double sigma, int max_iter) {
    if (!lc || lc->n == 0) return 0;
    double *flux = lc->flux;
    size_t n = lc->n;
    for (int iter = 0; iter < max_iter; iter++) {
        double med = median(flux, n);
        double mad_val = mad(flux, n, med);
        if (mad_val <= 0) break;
        double lo = med - sigma * mad_val;
        double hi = med + sigma * mad_val;
        size_t keep = 0;
        for (size_t i = 0; i < n; i++) {
            if (flux[i] >= lo && flux[i] <= hi) {
                if (keep != i) {
                    flux[keep] = flux[i];
                    lc->time[keep] = lc->time[i];
                    lc->flux_err[keep] = lc->flux_err[i];
                    lc->quality[keep] = lc->quality[i];
                }
                keep++;
            }
        }
        if (keep == n) break;
        n = keep;
    }
    lc->n = n;
    return n;
}

/* ── Savitzky-Golay flattening ─────────────────────────────────────────────── */

/* Savitzky-Golay flattening.
 *
 * Middle points: fixed-window convolution with the precomputed SG
 * coefficients (mode='interp' as in scipy.signal.savgol_filter).
 * Edge points: polynomial least-squares fit of the first/last `window_len`
 * samples on LOCAL x = 0..W-1 (scipy `_fit_edge` semantics), evaluated at
 * local indices 0..halflen-1 (left) and halflen..W-1 (right).
 *
 * Numerical stability: the local x coordinates are centred and scaled to
 * [-1, 1] before building the Vandermonde; a raw x in [0, 1620] makes the
 * QR factorisation lose all precision (Vandermonde cond ~ 1e15) while
 * numpy/scipy use SVD.  The centred fit is the same least-squares fit
 * (affine change of basis -> identical fitted values), just well
 * conditioned.
 *
 * Matches scipy.signal.savgol_filter(window_length, polyorder,
 * mode='interp') point-by-point. */

typedef struct {
    double *V;      /* factored V, W x (p+1), row-major, x in [-1, 1] */
    double *vnorm;  /* per-column ||v||^2 (implicit v[k]=1) */
    double rdiag[11]; /* true R diagonal = +/-||v|| per column (p <= 10) */
    double xmid;    /* (W-1)/2, for mapping local index -> centred x */
    double xscale;  /* (W-1)/2, so xc = (i - xmid)/xscale in [-1, 1] */
    size_t W;
    int p;
} SGQR;

static double sg_qr_x(const SGQR *q, long local_idx) {
    return ((double)local_idx - q->xmid) / q->xscale;
}

static int sg_qr_build(size_t W, int p, SGQR *q) {
    size_t m = W, n = (size_t)p + 1;
    q->V = malloc(m * n * sizeof(double));
    q->vnorm = malloc((size_t)p * sizeof(double));
    if (!q->V || !q->vnorm) { free(q->V); free(q->vnorm); q->V = NULL; q->vnorm = NULL; return -1; }
    q->W = W; q->p = p;
    q->xmid = 0.5 * (double)(W - 1);
    q->xscale = (W > 1) ? 0.5 * (double)(W - 1) : 1.0;
    for (size_t i = 0; i < m; i++) {
        double x = ((double)i - q->xmid) / q->xscale;
        double xp = 1.0;
        for (size_t j = 0; j < n; j++) { q->V[i * n + j] = xp; xp *= x; }
    }
    for (size_t k = 0; k < (size_t)p; k++) {
        double r_norm = 0.0;
        for (size_t i = k; i < m; i++) r_norm += q->V[i * n + k] * q->V[i * n + k];
        r_norm = sqrt(r_norm);
        if (r_norm == 0.0) return -1;
        double alpha = (q->V[k * n + k] >= 0) ? -r_norm : r_norm;
        q->rdiag[k] = alpha;
        q->V[k * n + k] -= alpha;
        for (size_t i = k + 1; i < m; i++) q->V[i * n + k] /= q->V[k * n + k];
        double vnorm2 = 1.0;
        for (size_t i = k + 1; i < m; i++) vnorm2 += q->V[i * n + k] * q->V[i * n + k];
        q->vnorm[k] = vnorm2;
        for (size_t j = k + 1; j < n; j++) {
            double dot = q->V[k * n + j];
            for (size_t i = k + 1; i < m; i++) dot += q->V[i * n + k] * q->V[i * n + j];
            double tau = 2.0 * dot / vnorm2;
            q->V[k * n + j] -= tau;
            for (size_t i = k + 1; i < m; i++) q->V[i * n + j] -= tau * q->V[i * n + k];
        }
    }
    /* The last column (k == p) never gets a reflector; its diagonal IS the
     * R diagonal for the back-substitution. */
    q->rdiag[p] = q->V[p * n + p];
    return 0;
}

/* Evaluate the edge fit at LOCAL index `local_idx` (0..W-1 within the
 * fitted window starting at data index n0). */
static double sg_qr_eval(const SGQR *q, const double *y, size_t n0, long local_idx) {
    size_t m = q->W, n = (size_t)q->p + 1;
    const double *V = q->V;
    double Qy[m];
    for (size_t i = 0; i < m; i++) Qy[i] = y[n0 + i];
    for (size_t k = 0; k < (size_t)q->p; k++) {
        double dotq = Qy[k];
        for (size_t i = k + 1; i < m; i++) dotq += V[i * n + k] * Qy[i];
        double tauq = 2.0 * dotq / q->vnorm[k];
        Qy[k] -= tauq;
        for (size_t i = k + 1; i < m; i++) Qy[i] -= tauq * V[i * n + k];
    }
    double c[n];
    for (long k = (long)n - 1; k >= 0; k--) {
        double s = Qy[k];
        for (size_t j = (size_t)k + 1; j < n; j++) s -= V[k * n + j] * c[j];
        c[k] = s / q->rdiag[k];
    }
    double x = sg_qr_x(q, local_idx);
    double xp = 1.0, val = 0.0;
    for (size_t j = 0; j < n; j++) { val += c[j] * xp; xp *= x; }
    return val;
}

static void sg_qr_free(SGQR *q) {
    free(q->V); free(q->vnorm);
    q->V = NULL; q->vnorm = NULL;
}

int zs_savgol_flatten(ZSLightCurve *lc, int window_len, int polyorder) {
    if (!lc || lc->n == 0) return -1;
    if (window_len % 2 == 0 || window_len < 3) return -1;
    if (polyorder < 0 || polyorder >= window_len) return -1;

    double *coeffs = NULL;
    if (sg_coeffs(window_len, polyorder, &coeffs) != 0) return -1;

    double *orig = malloc(lc->n * sizeof(double));
    if (!orig) { free(coeffs); return -1; }
    memcpy(orig, lc->flux, lc->n * sizeof(double));

    int half = window_len / 2;
    size_t edge = (size_t)half;
    if (edge > lc->n / 2) edge = lc->n / 2;

    SGQR q;
    q.V = q.vnorm = NULL;
    if (edge > 0 && sg_qr_build((size_t)window_len, polyorder, &q) != 0) {
        free(coeffs); free(orig);
        return -1;
    }

    /* Edge points (scipy _fit_edge): left fits data[0..W] evaluated at
     * local 0..halflen-1; right fits data[n-W..n] evaluated at local
     * halflen..W-1, both with x = arange(W). */
    if (edge > 0) {
        #pragma omp parallel for schedule(static) if (!omp_in_parallel() && window_len >= 129)
        for (long i = 0; i < (long)edge; i++)
            lc->flux[i] = sg_qr_eval(&q, orig, 0, i);
        #pragma omp parallel for schedule(static) if (!omp_in_parallel() && window_len >= 129)
        for (long j = 0; j < (long)edge; j++)
            lc->flux[lc->n - edge + (size_t)j] =
                sg_qr_eval(&q, orig, lc->n - (size_t)window_len, (long)half + 1 + j);
        sg_qr_free(&q);
    }

    /* Middle: fixed-window convolution. */
    #pragma omp parallel for schedule(static) if (!omp_in_parallel() && window_len >= 129)
    for (long i = (long)edge; i <= (long)(lc->n - 1 - edge); i++) {
        size_t start = (i >= (long)half) ? (size_t)i - (size_t)half : 0;
        double sum = 0;
        for (size_t k = start; k < start + (size_t)window_len; k++)
            sum += orig[k] * coeffs[k - start];
        lc->flux[i] = sum;
    }

    free(orig);
    free(coeffs);
    return 0;
}

/* ── Full ingestion pipeline ───────────────────────────────────────────────── */

int zs_ingest_fits(const char *path,
                   ZSLightCurve *out,
                   double sigma_clip,
                   int savgol_window,
                   int savgol_polyorder,
                   char *errbuf, size_t errbuf_sz) {
    if (zs_read_fits(path, out, errbuf, errbuf_sz) != 0) return -1;
    if (zs_mask_quality(out, 0) == 0) {
        if (errbuf) snprintf(errbuf, errbuf_sz, "all cadences masked");
        zs_lightcurve_free(out);
        return -1;
    }
    if (sigma_clip > 0) {
        zs_sigma_clip(out, sigma_clip, 5);
        if (out->n < 10) {
            if (errbuf) snprintf(errbuf, errbuf_sz, "too few points after sigma clip");
            zs_lightcurve_free(out);
            return -1;
        }
    }
    if (savgol_window > 0 && savgol_polyorder >= 0) {
        if (zs_savgol_flatten(out, savgol_window, savgol_polyorder) != 0) {
            if (errbuf) snprintf(errbuf, errbuf_sz, "SG flatten failed");
            zs_lightcurve_free(out);
            return -1;
        }
    }
    return 0;
}

/* ── Synthetic generation ──────────────────────────────────────────────────── */

int zs_generate_synthetic(double period_days, double snr,
                          double transit_depth, double duration_hrs,
                          double t0_days,
                          double baseline_days, double cadence_min,
                          ZSLightCurve *out) {
    if (!out) return -1;

    double cadence_days = cadence_min / (24.0 * 60.0);
    size_t n = (size_t)(baseline_days / cadence_days) + 1;
    if (n < 100) n = 100;

    out->time = malloc(n * sizeof(double));
    out->flux = malloc(n * sizeof(double));
    out->flux_err = malloc(n * sizeof(double));
    out->quality = malloc(n * sizeof(uint8_t));
    if (!out->time || !out->flux || !out->flux_err || !out->quality) {
        zs_lightcurve_free(out);
        return -1;
    }

    double duration_days = duration_hrs / 24.0;
    double depth = transit_depth;
    double noise = depth / snr;

    for (size_t i = 0; i < n; i++) {
        double t = i * cadence_days;
        out->time[i] = t;
        double phase = fmod(t - t0_days, period_days) / period_days;
        if (phase < 0) phase += 1.0;
        double in_transit = (phase < 0.5 * duration_days / period_days) ||
                            (phase > 1.0 - 0.5 * duration_days / period_days);
        double flux = 1.0 - (in_transit ? depth : 0.0);
        out->flux[i] = flux + ((double)rand() / RAND_MAX - 0.5) * 2.0 * noise;
        out->flux_err[i] = noise;
        out->quality[i] = 0;
    }

    out->n = n;
    out->t_start = 0;
    out->t_stop = baseline_days;
    out->cadence_days = cadence_days;
    return 0;
}