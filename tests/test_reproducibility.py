"""
Reproducibility guarantees for the v1 benchmark suite.

The controlled benchmark stands on these properties:

* same seed + same parameters -> byte-identical synthetic light curves;
* different seeds -> different noise realizations (no accidental aliasing);
* the deterministic sample builders in run_controlled.py pick identical
  targets when given the same seed (so benchmark runs are comparable).

These tests run offline: no lightkurve, no NEA, no engine invocation.
"""

import numpy as np
import pytest

from benchmarks_controlled import synthetic as S
from benchmarks_controlled import run_controlled as R


def test_same_seed_produces_byte_identical_lightcurve():
    t, q = S.make_time_axis(n_sectors=3, seed=11)
    f1 = S.make_noise(t, 150.0, seed=11)
    t2, _ = S.make_time_axis(n_sectors=3, seed=11)
    f2 = S.make_noise(t2, 150.0, seed=11)
    assert np.array_equal(f1, f2)


def test_different_seeds_produce_different_noise():
    t, _ = S.make_time_axis(n_sectors=3, seed=12)
    f1 = S.make_noise(t, 150.0, seed=12)
    f2 = S.make_noise(t, 150.0, seed=13)
    assert not np.array_equal(f1, f2)


def test_same_seed_produces_identical_true_planet():
    p1 = S.generate_true_planet(0, period_days=3.7, target_snr=20.0, seed=42)
    p2 = S.generate_true_planet(0, period_days=3.7, target_snr=20.0, seed=42)
    assert np.array_equal(p1.time, p2.time)
    assert np.array_equal(p1.flux, p2.flux)
    assert p1.subkind == p2.subkind
    assert p1.injected_depth == p2.injected_depth


def test_same_seed_produces_identical_false_eb():
    e1 = S.generate_false_eb(0, period_days=9.1, depth=0.01, seed=42)
    e2 = S.generate_false_eb(0, period_days=9.1, depth=0.01, seed=42)
    assert np.array_equal(e1.time, e2.time)
    assert np.array_equal(e1.flux, e2.flux)


def test_sample_builders_are_deterministic_for_seed():
    a1, b1 = R.build_true_set(6, seed=20260814), R.build_false_set(6, seed=20260814)
    a2, b2 = R.build_true_set(6, seed=20260814), R.build_false_set(6, seed=20260814)
    assert [t.tic_id for t in a1] == [t.tic_id for t in a2]
    assert [t.tic_id for t in b1] == [t.tic_id for t in b2]


def test_sample_builders_change_with_seed():
    a1 = R.build_true_set(6, seed=20260814)
    a2 = R.build_true_set(6, seed=20260815)
    assert [t.injected_depth for t in a1] != [t.injected_depth for t in a2]


def test_noise_only_lightcurve_has_no_injected_dip():
    t, _ = S.make_time_axis(n_sectors=3, seed=99)
    f = S.make_noise(t, 150.0, seed=99)
    med = np.median(f)
    # pure noise must not contain a >1% coherent dip at any single point
    assert (f < med * 0.99).sum() == 0 or not np.any(
        np.diff(np.where(f < med * 0.99)[0]) == 1)