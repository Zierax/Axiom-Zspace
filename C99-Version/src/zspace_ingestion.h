/*
 * zspace_ingestion.h  ·  Light curve ingestion (FITS I/O, quality masking,
 * sigma-clipping, Savitzky-Golay flattening)
 * Mirrors zspace_engine/ingestion.py functionality in C99.
 */
#ifndef ZSPACE_INGESTION_H
#define ZSPACE_INGESTION_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Light curve product ───────────────────────────────────────────────────── */
typedef struct {
    double *time;
    double *flux;
    double *flux_err;
    uint8_t *quality;
    size_t n;
    double t_start;
    double t_stop;
    double cadence_days;
} ZSLightCurve;

/* Free a ZSLightCurve allocated by ingestion functions */
void zs_lightcurve_free(ZSLightCurve *lc);

/* ── FITS reading (requires cfitsio) ───────────────────────────────────────── */
/* Read a TESS/Kepler FITS file into a ZSLightCurve.
 * Returns 0 on success, -1 on error (sets error message in optional buffer). */
int zs_read_fits(const char *path, ZSLightCurve *out, char *errbuf, size_t errbuf_sz);

/* ── Quality masking ───────────────────────────────────────────────────────── */
/* Retain only cadences where quality == 0 (or quality_mask & quality == 0).
 * Modifies the light curve in place, compacting arrays.
 * Returns new n, or 0 if all masked. */
size_t zs_mask_quality(ZSLightCurve *lc, uint32_t quality_mask);

/* ── Sigma clipping ────────────────────────────────────────────────────────── */
/* Iterative sigma clipping on flux. Clips values outside [median - sigma*MAD, median + sigma*MAD].
 * Modifies in place. Returns new n. */
size_t zs_sigma_clip(ZSLightCurve *lc, double sigma, int max_iter);

/* ── Savitzky-Golay flattening ─────────────────────────────────────────────── */
/* SG filter for detrending long-term systematics.
 * window_len must be odd. polyorder < window_len.
 * Modifies flux in place (replaces with flattened flux).
 * Returns 0 on success, -1 on error. */
int zs_savgol_flatten(ZSLightCurve *lc, int window_len, int polyorder);

/* ── Full ingestion pipeline ───────────────────────────────────────────────── */
/* Run the standard pipeline: read FITS -> quality mask -> sigma clip -> SG flatten.
 * Returns 0 on success, -1 on error. */
int zs_ingest_fits(const char *path,
                   ZSLightCurve *out,
                   double sigma_clip,
                   int savgol_window,
                   int savgol_polyorder,
                   char *errbuf, size_t errbuf_sz);

/* ── Synthetic light curve generation (for testing) ────────────────────────── */
/* Generate a synthetic light curve with an injected transit.
 * Mirrors benchmarks_controlled/synthetic.py generate_true_planet.
 * Returns 0 on success. */
int zs_generate_synthetic(double period_days, double snr,
                          double transit_depth, double duration_hrs,
                          double t0_days,
                          double baseline_days, double cadence_min,
                          ZSLightCurve *out);

#ifdef __cplusplus
}
#endif

#endif /* ZSPACE_INGESTION_H */