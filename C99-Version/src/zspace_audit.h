/*
 * zspace_audit.h  ·  Transit vitality audits (C99)
 * =================================================
 * Direct port of zspace_engine/auditors.py hot paths:
 *   - extract_individual_transit_depths   (template matched-filter)
 *   - even_odd_test                       (Welch t-test)
 *   - depth_consistency_score             (chi2-reduced CV)
 *   - secondary_eclipse_test              (phase-0.5 window)
 *   - ingress_egress_test                 (trapezoid fit)
 */
#ifndef ZSPACE_AUDIT_H
#define ZSPACE_AUDIT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double *depths;      /* per-transit depths      */
    double *depth_errs;  /* per-transit uncertainties */
    int    *ns;          /* transit numbers (0=first) */
    int     n;
} ZSExtractResult;

typedef struct {
    int    n_even, n_odd;
    double depth_even, depth_odd;
    double depth_even_err, depth_odd_err;   /* standard errors of the means */
    double delta_sigma;
    double t_stat, p_value;
    int    is_eb_flag;
} ZSEvenOdd;

typedef struct {
    int    n;
    double mean_depth, std_depth, cv;
    double sigma_med, chi2_red, s_depth;
} ZSDepthCons;

typedef struct {
    double primary_depth, secondary_depth;
    double secondary_ratio, secondary_snr;
    int    n_primary, n_secondary, n_oot;
    int    ok;
} ZSSecondary;

typedef struct {
    double depth_fit, ingress_fraction, flat_fraction;
    double ingress_hrs, flat_hrs;
    int    is_v_shape;
    int    fp_risk;        /* 0 LOW, 1 MEDIUM, 2 HIGH, -1 UNKNOWN */
    int    fit_ok;
} ZSIngressEgress;

/* extract_individual_transit_depths — template matched filter per epoch.
   Returns 0 on success (out->n may be 0 when insufficient coverage). */
int  zs_extract_depths(const double *time, size_t n_time, const double *flux,
                       double period, double t0, double duration,
                       ZSExtractResult *out);
void zs_extract_free(ZSExtractResult *r);

void zs_even_odd(const double *depths, const int *ns, int n, ZSEvenOdd *out);

void zs_depth_consistency(const double *depths, const double *depth_errs, int n,
                          ZSDepthCons *out);

void zs_secondary_eclipse(const double *time, size_t n_time, const double *flux,
                          double period, double t0, double duration,
                          ZSSecondary *out);

/* bin_phase/bin_flux are the 200-bin full-orbit fold (NaN cells skipped). */
void zs_ingress_egress(const double *bin_phase, const double *bin_flux,
                       size_t n_bins, double period, double duration,
                       double transit_depth, ZSIngressEgress *out);

#ifdef __cplusplus
}
#endif

#endif /* ZSPACE_AUDIT_H */