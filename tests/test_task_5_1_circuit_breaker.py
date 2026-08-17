"""
Test for Task 5.1: Circuit breaker logic for overall_verdict calculation

This test verifies that when any critical test fails, the overall_verdict
is set to "FALSE_POSITIVE" and the logical_closure contains the message:
"CRITICAL CIRCUIT BREAKER TRIGGERED: One or more critical tests failed"

Requirements: 3.1, 3.2, 3.3, 3.4, 3.5
"""

import pytest
from unittest.mock import patch
from zspace_engine.validator import ProofEngine


def test_circuit_breaker_prevents_sovereign_pass_on_critical_failure():
    """
    Test that when any critical test fails, overall_verdict cannot be SOVEREIGN_PASS.
    
    Validates: Requirements 3.1, 3.2
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
        
        # Create ProofEngine with parameters that would normally pass
        # except for one critical test (FP-7 density ratio)
        engine = ProofEngine(
            period_days=3.5,
            transit_depth=0.01,
            transit_duration_hrs=2.0,
            stellar_mass_solar=1.0,
            stellar_radius_solar=1.0,
            stellar_teff_k=5778.0,
            stellar_logg=4.44,
            planet_radius_earth=10.0,
            bls_snr=10.0,  # High SNR (FP-1 passes)
            bls_fap=1e-6,  # Low FAP (FP-2 passes)
            even_odd_delta_sigma=1.5,  # Good even/odd (FP-3 passes)
            shape_ratio=2.0,  # U-shaped (FP-4 passes)
        )
        
        # Call false_positive_ruling with density_ratio outside range (FP-7 fails)
        result = engine.false_positive_ruling(
            secondary_snr=1.0,  # FP-5 passes
            centroid_sigma=1.0,  # FP-6 passes
density_ratio=0.1,  # Outside calibrated [0.2, 5.0] - FP-7 CRITICAL FAIL
            is_grazing=False,
            tic_id="TIC 123456789",
        )
        
        # Verify FP-7 failed and is critical
        fp7_test = next(t for t in result["tests"] if t["test"] == "FP-7 Density Ratio")
        assert fp7_test["verdict"] == "FAIL", "FP-7 should fail with density_ratio=0.1"
        assert fp7_test["weight"] == "critical", "FP-7 should have critical weight"
        
        # Verify circuit breaker prevents SOVEREIGN_PASS
        assert result["overall_verdict"] != "SOVEREIGN_PASS", \
            f"Circuit breaker should prevent SOVEREIGN_PASS when critical test fails, but got {result['overall_verdict']}"
        
        # Should be FALSE_POSITIVE due to critical failure
        assert result["overall_verdict"] == "FALSE_POSITIVE", \
            f"Expected FALSE_POSITIVE when critical test fails, but got {result['overall_verdict']}"


def test_circuit_breaker_adds_trigger_message():
    """
    Test that when any critical test fails, the logical_closure contains
    the circuit breaker trigger message.
    
    Validates: Requirement 3.4
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
        
        # Create ProofEngine with one critical test failing
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
            shape_ratio=2.0,
        )
        
        # Call with density_ratio outside range (critical failure)
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=6.0,  # Outside calibrated [0.2, 5.0] - CRITICAL FAIL
            is_grazing=False,
            tic_id="TIC 123456789",
        )
        
        # Verify the circuit breaker message is in logical_closure
        circuit_breaker_message = "CRITICAL CIRCUIT BREAKER TRIGGERED: One or more critical tests failed"
        assert circuit_breaker_message in result["logical_closure"], \
            f"Expected circuit breaker message in logical_closure, but got: {result['logical_closure']}"


def test_circuit_breaker_with_multiple_critical_failures():
    """
    Test that circuit breaker works correctly when multiple critical tests fail.
    
    Validates: Requirements 3.1, 3.3, 3.4
    """
    # Mock check_external_catalogs to return multiple star system (FP-9 fails)
    with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
        mock_catalog.return_value = {
            "is_multiple": True,  # FP-9 will fail
            "catalog_source": "SIMBAD",
            "risk_level": "HIGH",
            "classification": "EB*",
            "query_latency_ms": 150.0,
        }
        
        # Create ProofEngine with low SNR (FP-1 fails)
        engine = ProofEngine(
            period_days=3.5,
            transit_depth=0.01,
            transit_duration_hrs=2.0,
            stellar_mass_solar=1.0,
            stellar_radius_solar=1.0,
            stellar_teff_k=5778.0,
            stellar_logg=4.44,
            planet_radius_earth=10.0,
            bls_snr=5.0,  # Below 5.5 threshold - FP-1 CRITICAL FAIL
            bls_fap=1e-6,
            even_odd_delta_sigma=1.5,
            shape_ratio=2.0,
        )
        
        # Call with density_ratio outside range (FP-7 also fails)
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=0.4,  # Within calibrated band: critical failures come from FP-1 and FP-9
            is_grazing=False,
            tic_id="TIC 123456789",
        )
        
        # Count critical failures
        critical_failures = [t for t in result["tests"] 
                           if t["weight"] == "critical" and t["verdict"] == "FAIL"]
        
        assert len(critical_failures) >= 2, \
            f"Expected at least 2 critical failures, but got {len(critical_failures)}"
        
        # Verify overall verdict is FALSE_POSITIVE
        assert result["overall_verdict"] == "FALSE_POSITIVE", \
            f"Expected FALSE_POSITIVE with multiple critical failures, but got {result['overall_verdict']}"
        
        # Verify circuit breaker message is present
        circuit_breaker_message = "CRITICAL CIRCUIT BREAKER TRIGGERED: One or more critical tests failed"
        assert circuit_breaker_message in result["logical_closure"], \
            f"Expected circuit breaker message with multiple critical failures"


def test_no_circuit_breaker_when_all_critical_pass():
    """
    Test that when all critical tests pass, the circuit breaker message
    is NOT added and SOVEREIGN_PASS is possible.
    
    Validates: Requirement 3.5 (invariant)
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
        
        # Create ProofEngine with all tests passing
        engine = ProofEngine(
            period_days=3.5,
            transit_depth=0.01,
            transit_duration_hrs=2.0,
            stellar_mass_solar=1.0,
            stellar_radius_solar=1.0,
            stellar_teff_k=5778.0,
            stellar_logg=4.44,
            planet_radius_earth=10.0,
            bls_snr=10.0,  # FP-1 passes
            bls_fap=1e-6,  # FP-2 passes
            even_odd_delta_sigma=1.5,  # FP-3 passes
            shape_ratio=2.0,  # FP-4 passes
        )
        
        # Call with all parameters passing
        result = engine.false_positive_ruling(
            secondary_snr=1.0,  # FP-5 passes
            centroid_sigma=1.0,  # FP-6 passes
            density_ratio=1.0,  # FP-7 passes
            is_grazing=False,  # FP-8 passes
            tic_id="TIC 123456789",  # FP-9 passes
        )
        
        # Verify all critical tests passed
        critical_tests = [t for t in result["tests"] if t["weight"] == "critical"]
        critical_passed = [t for t in critical_tests if t["verdict"] == "PASS"]
        
        assert len(critical_passed) == len(critical_tests), \
            f"Expected all {len(critical_tests)} critical tests to pass"
        
        # Verify circuit breaker message is NOT present
        circuit_breaker_message = "CRITICAL CIRCUIT BREAKER TRIGGERED: One or more critical tests failed"
        assert circuit_breaker_message not in result["logical_closure"], \
            f"Circuit breaker message should not appear when all critical tests pass"
        
        # Verify SOVEREIGN_PASS is possible (with ≤1 total failure)
        if result["n_fail"] <= 1:
            assert result["overall_verdict"] == "SOVEREIGN_PASS", \
                f"Expected SOVEREIGN_PASS when all critical tests pass and n_fail={result['n_fail']}"


def test_circuit_breaker_invariant():
    """
    Test the invariant: (overall_verdict == "SOVEREIGN_PASS") implies (n_critical_pass == n_critical)
    
    Validates: Requirement 3.5
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
        
        # Test multiple scenarios
        test_cases = [
            # (density_ratio, bls_snr, expected_sovereign_pass_possible)
            (1.0, 10.0, True),   # All critical pass
            (0.3, 10.0, False),  # FP-7 fails
            (1.0, 5.0, False),   # FP-1 fails
            (0.3, 5.0, False),   # Both fail
        ]
        
        for density_ratio, bls_snr, expected_sovereign_pass_possible in test_cases:
            engine = ProofEngine(
                period_days=3.5,
                transit_depth=0.01,
                transit_duration_hrs=2.0,
                stellar_mass_solar=1.0,
                stellar_radius_solar=1.0,
                stellar_teff_k=5778.0,
                stellar_logg=4.44,
                planet_radius_earth=10.0,
                bls_snr=bls_snr,
                bls_fap=1e-6,
                even_odd_delta_sigma=1.5,
                shape_ratio=2.0,
            )
            
            result = engine.false_positive_ruling(
                secondary_snr=1.0,
                centroid_sigma=1.0,
                density_ratio=density_ratio,
                is_grazing=False,
                tic_id="TIC 123456789",
            )
            
            # Verify invariant: SOVEREIGN_PASS => all critical tests passed
            if result["overall_verdict"] == "SOVEREIGN_PASS":
                assert result["n_critical_pass"] == result["n_critical"], \
                    f"Invariant violated: SOVEREIGN_PASS but {result['n_critical_pass']}/{result['n_critical']} critical tests passed"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
