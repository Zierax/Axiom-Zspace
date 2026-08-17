"""
Test for Task 4.1: Verify FP-4 test behavior at the calibrated boundary.

The v1 calibration (THRESHOLDS_REPORT.md) set fp4_shape_min = 0.4, i.e.
FP-4 FAILS for shape_ratio <= 0.4 (deep V-shaped transits, characteristic
of EBs) and PASSES for shape_ratio > 0.4 (U-shaped, planet-like).

Requirements: 2.1, 2.5
"""

import pytest
from zspace_engine.validator import ProofEngine


def test_fp4_fails_when_shape_ratio_equals_floor():
    """
    Test that FP-4 fails at the calibrated boundary (shape_ratio = 0.4).

    Validates: Requirement 2.1
    """
    # Create ProofEngine with shape_ratio = 0.4 (boundary = floor)
    engine = ProofEngine(
        period_days=3.5,
        transit_depth=0.01,
        transit_duration_hrs=2.0,
        stellar_mass_solar=1.0,
        stellar_radius_solar=1.0,
        stellar_teff_k=5778.0,
        stellar_logg=4.44,
        planet_radius_earth=10.0,
        bls_snr=10.0,
        bls_fap=1e-6,
        even_odd_delta_sigma=1.5,
        shape_ratio=0.4,  # Boundary case (strictly-greater comparison)
    )
    
    # Call false_positive_ruling
    result = engine.false_positive_ruling(
        secondary_snr=1.0,
        centroid_sigma=1.0,
        density_ratio=1.0,
        is_grazing=False,
        tic_id="TIC 123456789",
    )
    
    # Find FP-4 test
    fp4_test = next(t for t in result["tests"] if t["test"] == "FP-4 Shape Ratio (U vs V)")
    
    # Verify FP-4 fails at the calibrated floor
    assert fp4_test["verdict"] == "FAIL", f"FP-4 should FAIL at shape_ratio = 0.4, but got {fp4_test['verdict']}"
    assert fp4_test["value"] == 0.4
    assert fp4_test["threshold"] == 0.4


def test_fp4_fails_when_shape_ratio_below_floor():
    """
    Test that FP-4 fails when shape_ratio < 0.4 (deep V-shaped transit).

    Validates: Requirement 2.1
    """
    # Create ProofEngine with shape_ratio = 0.2 (V-shaped)
    engine = ProofEngine(
        period_days=3.5,
        transit_depth=0.01,
        transit_duration_hrs=2.0,
        stellar_mass_solar=1.0,
        stellar_radius_solar=1.0,
        stellar_teff_k=5778.0,
        stellar_logg=4.44,
        planet_radius_earth=10.0,
        bls_snr=10.0,
        bls_fap=1e-6,
        even_odd_delta_sigma=1.5,
        shape_ratio=0.2,  # V-shaped (EB-like)
    )
    
    # Call false_positive_ruling
    result = engine.false_positive_ruling(
        secondary_snr=1.0,
        centroid_sigma=1.0,
        density_ratio=1.0,
        is_grazing=False,
        tic_id="TIC 123456789",
    )
    
    # Find FP-4 test
    fp4_test = next(t for t in result["tests"] if t["test"] == "FP-4 Shape Ratio (U vs V)")
    
    # Verify FP-4 fails when shape_ratio < 0.4
    assert fp4_test["verdict"] == "FAIL", f"FP-4 should FAIL when shape_ratio < 0.4, but got {fp4_test['verdict']}"
    assert fp4_test["value"] == 0.2


def test_fp4_passes_when_shape_ratio_above_floor():
    """
    Test that FP-4 passes when shape_ratio > 0.4 (U-shaped transit).

    Validates: Requirement 2.5
    """
    # Create ProofEngine with shape_ratio = 2.0 (U-shaped)
    engine = ProofEngine(
        period_days=3.5,
        transit_depth=0.01,
        transit_duration_hrs=2.0,
        stellar_mass_solar=1.0,
        stellar_radius_solar=1.0,
        stellar_teff_k=5778.0,
        stellar_logg=4.44,
        planet_radius_earth=10.0,
        bls_snr=10.0,
        bls_fap=1e-6,
        even_odd_delta_sigma=1.5,
        shape_ratio=2.0,  # U-shaped (planet-like)
    )
    
    # Call false_positive_ruling
    result = engine.false_positive_ruling(
        secondary_snr=1.0,
        centroid_sigma=1.0,
        density_ratio=1.0,
        is_grazing=False,
        tic_id="TIC 123456789",
    )
    
    # Find FP-4 test
    fp4_test = next(t for t in result["tests"] if t["test"] == "FP-4 Shape Ratio (U vs V)")
    
    # Verify FP-4 passes when shape_ratio > 0.4
    assert fp4_test["verdict"] == "PASS", f"FP-4 should PASS when shape_ratio > 0.4, but got {fp4_test['verdict']}"
    assert fp4_test["value"] == 2.0


def test_fp4_multiple_shape_ratios():
    """
    Test FP-4 at and around the calibrated 0.4 floor.

    Validates: Requirements 2.1, 2.5
    """
    test_cases = [
        (0.2, "FAIL"),   # deep V-shaped
        (0.3, "FAIL"),   # V-shaped
        (0.39, "FAIL"),  # just below floor
        (0.4, "FAIL"),   # boundary (strictly-greater comparison)
        (0.41, "PASS"),  # just above floor
        (1.1, "PASS"),   # U-shaped
        (2.0, "PASS"),   # U-shaped
        (4.5, "PASS"),   # U-shaped
    ]
    
    for shape_ratio, expected_verdict in test_cases:
        engine = ProofEngine(
            period_days=3.5,
            transit_depth=0.01,
            transit_duration_hrs=2.0,
            stellar_mass_solar=1.0,
            stellar_radius_solar=1.0,
            stellar_teff_k=5778.0,
            stellar_logg=4.44,
            planet_radius_earth=10.0,
            bls_snr=10.0,
            bls_fap=1e-6,
            even_odd_delta_sigma=1.5,
            shape_ratio=shape_ratio,
        )
        
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=1.0,
            is_grazing=False,
            tic_id="TIC 123456789",
        )
        
        fp4_test = next(t for t in result["tests"] if t["test"] == "FP-4 Shape Ratio (U vs V)")
        
        assert fp4_test["verdict"] == expected_verdict, \
            f"FP-4 with shape_ratio={shape_ratio} should be {expected_verdict}, but got {fp4_test['verdict']}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])