"""
BLS Kernels for Purce translation
=================================
One line = one numpy call. No arithmetic inside call arguments.
No nested calls. No BinOp nodes inside arguments.
Only numpy ops supported by purce v0.1.0:
  add, subtract, multiply, divide, power, sqrt, square, abs,
  maximum, negative, exp, log, sin, cos, tan,
  linspace, sum, max, min, argmax, argmin, mean, std
  (no fmod, no logical_or, no where=, no out=, no nan_to_num,
   no isfinite, no logspace, no arange, no argsort)

The folding/scanning hot loop lives in zspace_bls.c (pure C);
these kernels cover grid generation + scalar statistics.
"""
import numpy as np

# ═══════════════════════════════════════════════════════════════════════════════
# FREQUENCY GRID
# ═══════════════════════════════════════════════════════════════════════════════

def bls_freq_grid(f_min, f_max, n_freq):
    """Linear frequency grid: f_min + i * df"""
    return np.linspace(f_min, f_max, n_freq)

def bls_duration_grid(dur_min, dur_max, n_dur):
    """Linear duration grid"""
    return np.linspace(dur_min, dur_max, n_dur)

def bls_period_from_freq(freq):
    """period = 1 / freq"""
    return np.divide(1.0, freq)

# ═══════════════════════════════════════════════════════════════════════════════
# SCALAR STATISTICS (called once per phase window)
# ═══════════════════════════════════════════════════════════════════════════════

def bls_sum(arr):
    """Sum of array"""
    return np.sum(arr)

def bls_max(arr):
    """Maximum value"""
    return np.max(arr)

def bls_argmax(arr):
    """Index of maximum (returns int index)"""
    return np.argmax(arr)

def bls_mean(arr):
    """Mean of array"""
    return np.mean(arr)

def bls_snr_from_stats(depth, std_out, n_in, n_out):
    """SNR = depth / std_out * sqrt(n_in * n_out / n_total)
    n_total is passed as n_out (caller handles ratio)"""
    return np.multiply(np.divide(depth, std_out, where=std_out > 0), np.sqrt(n_in * n_out))

# ═══════════════════════════════════════════════════════════════════════════════
# FAP MODEL (exponential tail)
# ═══════════════════════════════════════════════════════════════════════════════

def bls_fap_exponential(power, a, b):
    """FAP = exp(-a * power + b)"""
    return np.exp(np.add(np.multiply(a, power), b))

def bls_fap_threshold(fap_target, a, b):
    """Inverse: power threshold for given FAP"""
    return np.divide(np.add(b, np.log(fap_target)), a)

# ═══════════════════════════════════════════════════════════════════════════════
# HARMONIC REJECTION (elementwise on period arrays)
# ═══════════════════════════════════════════════════════════════════════════════

def bls_ratio(p1, p2):
    """p1 / p2 elementwise"""
    return np.divide(p1, p2)

def bls_round(arr):
    """Round to nearest integer"""
    return np.round(arr)

def bls_abs(arr):
    """Absolute value"""
    return np.abs(arr)

def bls_less_than(a, b):
    """a < b elementwise (1.0/0.0)"""
    return np.less(a, b).astype(np.float64)

# ═══════════════════════════════════════════════════════════════════════════════
# TRANSIT DURATION FROM PHYSICS (elementwise)
# ═══════════════════════════════════════════════════════════════════════════════

def bls_one_plus_rp(rp_rs):
    """1 + rp"""
    return np.add(1.0, rp_rs)

def bls_sq(x):
    """x^2"""
    return np.multiply(x, x)

def bls_one_minus_sq_b(sq_b):
    """1 - b^2"""
    return np.subtract(1.0, sq_b)

def bls_sqrt(x):
    """sqrt(x)"""
    return np.sqrt(x)

def bls_dur_term(term, a_rs):
    """term / a_rs"""
    return np.divide(term, a_rs)

def bls_dur_from_phys(period_days, dur_term):
    """T = P * dur_term / pi"""
    return np.multiply(period_days, np.divide(dur_term, np.pi))

def bls_var_approx(mean, mean_sq):
    """variance approx: E[x^2] - E[x]^2"""
    return np.subtract(mean_sq, np.multiply(mean, mean))