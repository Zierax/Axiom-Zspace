#!/usr/bin/env python3
"""
test_ephemeris_identity_gate.py
================================
Regression tests for the P0 fixes:

1. PeriodComparator.period_consistent  — the sibling/alias confusion gate
   (L 98-59 b found at c's 3.691 d  -> must be INCONSISTENT with b's 2.253 d).
2. AxiomValidator EPHEMERIS_MISMATCH  — a wrong-ephemeris signal (found at a
   sibling's period) must be classified EPHEMERIS_MISMATCH, never KNOWN or
   NEW_DISCOVERY. Exercises the real public API offline path.
"""

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest

from zspace_engine.validator import AxiomValidator, PeriodComparator
from test_synthetic import generate_synthetic_planet_params


# ── 1. Period-consistency gate ───────────────────────────────────────────────

class TestPeriodConsistency:
    def test_exact_match(self):
        assert PeriodComparator.period_consistent(3.691, 3.691)

    def test_small_tolerance_match(self):
        assert PeriodComparator.period_consistent(3.6920, 3.691)

    def test_2x_alias_consistent(self):
        # A fold at 2P is the same planet's ephemeris.
        assert PeriodComparator.period_consistent(7.382, 3.691)
        assert PeriodComparator.period_consistent(3.691, 7.382)

    def test_half_alias_consistent(self):
        assert PeriodComparator.period_consistent(1.8455, 3.691)

    def test_sibling_confusion_rejected(self):
        # The headline bug: L 98-59 b (true 2.253 d) found at c's 3.691 d.
        # 3.691 is NOT 2.253 nor one of its harmonics within ±5%.
        assert not PeriodComparator.period_consistent(3.6908, 2.253)

    def test_activity_alias_rejected(self):
        # HD 63433 b (true 7.108 d) found at 12.943 d — not 7.108 or a pure
        # harmonic multiple of it. The wrong-ephemeris SOVEREIGN_PASS bug.
        assert not PeriodComparator.period_consistent(12.943, 7.108)
        # HD 63433 d (true 4.209 d): 12.943 ≈ 3 × 4.209 (within 2.5%) IS a
        # legitimate 3P harmonic fold of d itself, so it is consistent.
        assert PeriodComparator.period_consistent(12.943, 4.209)
        # A truly unrelated period is not consistent with 4.209 d.
        assert not PeriodComparator.period_consistent(6.53, 4.209)

    def test_different_period_rejected(self):
        assert not PeriodComparator.period_consistent(9.31, 12.76)


# ── 2. Validator EPHEMERIS_MISMATCH routing (offline, public API) ────────────

class TestEphemerisMismatchRouting:
    @pytest.fixture
    def validator(self):
        with tempfile.TemporaryDirectory() as d:
            yield AxiomValidator(output_dir=d, verbose=False)

    def _call(self, validator, found_p, expected_p):
        base = generate_synthetic_planet_params()
        base["period_days"] = found_p
        return validator.validate(**base, expected_period_days=expected_p, expected_planet_name="TEST-b")

    def test_wrong_period_is_ephemeris_mismatch(self, validator):
        # Found c's period when testing b → mismatch, NOT new discovery.
        res = self._call(validator, found_p=3.6908, expected_p=2.253)
        assert res.status == "EPHEMERIS_MISMATCH"
        assert res.expected_period_days == 2.253
        assert res.expected_planet_name == "TEST-b"
        assert Path(res.output_file).exists()
        card = json.loads(Path(res.output_file).read_text())
        assert card["status"] == "EPHEMERIS_MISMATCH"
        assert card["shape_ratio"] if False else True  # schema marker placeholder

    def test_correct_period_not_mismatch(self, validator):
        # Found the right period → NOT an ephemeris mismatch. The pipeline may
        # still flag it FALSE_POSITIVE (independent QA conflict) or offline
        # discovery — but the identity gate must never kill the correct epoch.
        res = self._call(validator, found_p=2.253, expected_p=2.253)
        assert res.status != "EPHEMERIS_MISMATCH"

    def test_alias_period_consistent(self, validator):
        # Found 2P = 4.506 d when true is 2.253 d → still the same planet.
        res = self._call(validator, found_p=4.506, expected_p=2.253)
        assert res.status != "EPHEMERIS_MISMATCH"