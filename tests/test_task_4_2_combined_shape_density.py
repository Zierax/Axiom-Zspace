"""
Test for Task 4.2: Combined shape+density failure detection

This test verifies that when both FP-4 (shape ratio) and FP-7 (density ratio) fail,
the logical_closure contains the message:
"SHAPE+DENSITY CONFLICT: V-shaped transit with density mismatch indicates EB"

Also verifies that when FP-7 is CRITICAL_FAIL (density outside range),
the overall verdict is FALSE_POSITIVE (not SOVEREIGN_FAIL).

Requirements: 2.3, 2.4
"""

import pytest
from unittest.mock import patch
from zspace_engine.validator import ProofEngine


def test_combined_shape_density_failure_adds_conflict_message():
    """
    Test that when both FP-4 and FP-7 fail, the logical_closure contains
    the SHAPE+DENSITY CONFLICT message.
    
    Validates: Requirement 2.4
    """
    # Mock check_external_catalogs to avoid network calls
    with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
        mock_catalog.return_value = {
            "is_multiple": False,
            "catalog_source": "None",
            "risk_level": "LOW",
            "classification": "Single",
            "query_latency_ms": 0.0,
        }
        
        # Create ProofEngine with:
        # - shape_ratio = 0.2 (< 0.4 calibrated floor, so FP-4 fails)
        # - density_ratio will be 6.0 (> 5.0 calibrated band, so FP-7 fails)
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
            shape_ratio=0.2,  # V-shaped (FP-4 will fail)
        )
        
        # Call false_positive_ruling with density_ratio = 6.0 (FP-7 will fail)
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=6.0,  # Outside calibrated [0.2, 5.0] band
            is_grazing=False,
            tic_id="TIC 123456789",
        )
        
        # Verify both FP-4 and FP-7 failed
        fp4_test = next(t for t in result["tests"] if t["test"] == "FP-4 Shape Ratio (U vs V)")
        fp7_test = next(t for t in result["tests"] if t["test"] == "FP-7 Density Ratio")
        
        assert fp4_test["verdict"] == "FAIL", "FP-4 should fail with shape_ratio=0.2"
        assert fp7_test["verdict"] == "FAIL", "FP-7 should fail with density_ratio=6.0"
        
        # Verify the SHAPE+DENSITY CONFLICT message is in logical_closure
        conflict_message = "SHAPE+DENSITY CONFLICT: V-shaped transit with density mismatch indicates EB"
        assert conflict_message in result["logical_closure"], \
            f"Expected conflict message in logical_closure, but got: {result['logical_closure']}"


def test_fp7_critical_fail_results_in_false_positive():
    """
    Test that when FP-7 is CRITICAL_FAIL (density outside range),
    the overall verdict is FALSE_POSITIVE (not SOVEREIGN_FAIL).
    
    Validates: Requirement 2.3
    """
    # Mock check_external_catalogs to avoid network calls
    with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
        mock_catalog.return_value = {
            "is_multiple": False,
            "catalog_source": "None",
            "risk_level": "LOW",
            "classification": "Single",
            "query_latency_ms": 0.0,
        }
        
        # Create ProofEngine with:
        # - shape_ratio = 0.8 (above 0.4 calibrated floor, so FP-4 passes)
        # - All other tests will pass
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
            shape_ratio=0.8,
        )
        
        # Call false_positive_ruling with density_ratio = 6.0 (FP-7 critical fail)
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=6.0,  # Outside calibrated [0.2, 5.0] band (critical fail)
            is_grazing=False,
            tic_id="TIC 123456789",
        )
        
        # Verify FP-7 is critical and failed
        fp7_test = next(t for t in result["tests"] if t["test"] == "FP-7 Density Ratio")
        assert fp7_test["verdict"] == "FAIL", "FP-7 should fail with density_ratio=6.0"
        assert fp7_test["weight"] == "critical", "FP-7 should have critical weight"
        
        # Verify overall verdict is FALSE_POSITIVE (not SOVEREIGN_FAIL)
        assert result["overall_verdict"] == "FALSE_POSITIVE", \
            f"Expected FALSE_POSITIVE when FP-7 critical fails, but got {result['overall_verdict']}"


def test_no_conflict_message_when_only_fp4_fails():
    """
    Test that the SHAPE+DENSITY CONFLICT message is NOT added when only FP-4 fails.
    """
    # Mock check_external_catalogs to avoid network calls
    with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
        mock_catalog.return_value = {
            "is_multiple": False,
            "catalog_source": "None",
            "risk_level": "LOW",
            "classification": "Single",
            "query_latency_ms": 0.0,
        }
        
        # Create ProofEngine with:
        # - shape_ratio = 0.2 (< 0.4 calibrated floor, so FP-4 fails)
        # - density_ratio will be 1.0 (within calibrated band, so FP-7 passes)
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
            shape_ratio=0.2,  # V-shaped (FP-4 will fail)
        )
        
        # Call false_positive_ruling with density_ratio = 1.0 (FP-7 will pass)
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=1.0,  # Within calibrated [0.2, 5.0] band
            is_grazing=False,
            tic_id="TIC 123456789",
        )
        
        # Verify FP-4 failed but FP-7 passed
        fp4_test = next(t for t in result["tests"] if t["test"] == "FP-4 Shape Ratio (U vs V)")
        fp7_test = next(t for t in result["tests"] if t["test"] == "FP-7 Density Ratio")
        
        assert fp4_test["verdict"] == "FAIL", "FP-4 should fail with shape_ratio=0.2"
        assert fp7_test["verdict"] == "PASS", "FP-7 should pass with density_ratio=1.0"
        
        # Verify the SHAPE+DENSITY CONFLICT message is NOT in logical_closure
        conflict_message = "SHAPE+DENSITY CONFLICT: V-shaped transit with density mismatch indicates EB"
        assert conflict_message not in result["logical_closure"], \
            f"Conflict message should not appear when only FP-4 fails"


def test_no_conflict_message_when_only_fp7_fails():
    """
    Test that the SHAPE+DENSITY CONFLICT message is NOT added when only FP-7 fails.
    """
    # Mock check_external_catalogs to avoid network calls
    with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
        mock_catalog.return_value = {
            "is_multiple": False,
            "catalog_source": "None",
            "risk_level": "LOW",
            "classification": "Single",
            "query_latency_ms": 0.0,
        }
        
        # Create ProofEngine with:
        # - shape_ratio = 2.0 (above floor, so FP-4 passes)
        # - density_ratio will be 0.1 (< 0.2 calibrated band, so FP-7 fails)
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
            shape_ratio=2.0,  # U-shaped (FP-4 will pass)
        )
        
        # Call false_positive_ruling with density_ratio = 0.1 (FP-7 will fail)
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=0.1,  # Outside calibrated [0.2, 5.0] band
            is_grazing=False,
            tic_id="TIC 123456789",
        )
        
        # Verify FP-4 passed but FP-7 failed
        fp4_test = next(t for t in result["tests"] if t["test"] == "FP-4 Shape Ratio (U vs V)")
        fp7_test = next(t for t in result["tests"] if t["test"] == "FP-7 Density Ratio")
        
        assert fp4_test["verdict"] == "PASS", "FP-4 should pass with shape_ratio=2.0"
        assert fp7_test["verdict"] == "FAIL", "FP-7 should fail with density_ratio=0.1"
        
        # Verify the SHAPE+DENSITY CONFLICT message is NOT in logical_closure
        conflict_message = "SHAPE+DENSITY CONFLICT: V-shaped transit with density mismatch indicates EB"
        assert conflict_message not in result["logical_closure"], \
            f"Conflict message should not appear when only FP-7 fails"


def test_both_tests_pass_no_conflict_message():
    """
    Test that the SHAPE+DENSITY CONFLICT message is NOT added when both tests pass.
    """
    # Mock check_external_catalogs to avoid network calls
    with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
        mock_catalog.return_value = {
            "is_multiple": False,
            "catalog_source": "None",
            "risk_level": "LOW",
            "classification": "Single",
            "query_latency_ms": 0.0,
        }
        
        # Create ProofEngine with:
        # - shape_ratio = 2.0 (> 1.0, so FP-4 passes)
        # - density_ratio will be 1.0 (within [0.5, 2.0], so FP-7 passes)
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
            shape_ratio=2.0,  # U-shaped (FP-4 will pass)
        )
        
        # Call false_positive_ruling with density_ratio = 1.0 (FP-7 will pass)
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=1.0,  # Within [0.5, 2.0] range
            is_grazing=False,
            tic_id="TIC 123456789",
        )
        
        # Verify both tests passed
        fp4_test = next(t for t in result["tests"] if t["test"] == "FP-4 Shape Ratio (U vs V)")
        fp7_test = next(t for t in result["tests"] if t["test"] == "FP-7 Density Ratio")
        
        assert fp4_test["verdict"] == "PASS", "FP-4 should pass with shape_ratio=2.0"
        assert fp7_test["verdict"] == "PASS", "FP-7 should pass with density_ratio=1.0"
        
        # Verify the SHAPE+DENSITY CONFLICT message is NOT in logical_closure
        conflict_message = "SHAPE+DENSITY CONFLICT: V-shaped transit with density mismatch indicates EB"
        assert conflict_message not in result["logical_closure"], \
            f"Conflict message should not appear when both tests pass"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
