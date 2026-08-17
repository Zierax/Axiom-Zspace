"""
Integration Test for Task 7.1: Conflict Detection

This test demonstrates the complete conflict detection workflow with realistic
scenarios showing how SNR vs density/shape conflicts are detected and logged.
"""

from unittest.mock import patch
from zspace_engine.validator import ProofEngine


def test_high_snr_eb_with_conflicts():
    """
    Integration test: High SNR signal that is actually an EB.
    
    Scenario: A candidate has strong detection (SNR=25.0) but shows clear
    EB indicators (V-shaped transit, poor density match). The system should:
    1. Detect both SNR_DENSITY_CONFLICT and SNR_SHAPE_CONFLICT
    2. Set overall_verdict to FALSE_POSITIVE (due to critical test failures)
    3. Include conflict metadata with timestamps and latency
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
        
        # Create ProofEngine with EB-like parameters but high SNR
        engine = ProofEngine(
            period_days=1.5,
            transit_depth=0.02,
            transit_duration_hrs=1.5,
            stellar_mass_solar=1.2,
            stellar_radius_solar=1.1,
            stellar_teff_k=6000.0,
            stellar_logg=4.3,
            planet_radius_earth=12.0,
            bls_snr=25.0,  # Very high SNR - strong detection
            bls_fap=1e-8,
            even_odd_delta_sigma=0.8,
            shape_ratio=0.3,  # Deep V-shaped transit: fails FP-4 (floor 0.4)
        )
        
        # Call false_positive_ruling with poor density match
        result = engine.false_positive_ruling(
            secondary_snr=0.5,
            centroid_sigma=0.5,
            density_ratio=6.0,  # Way outside calibrated [0.2, 5.0] band - EB indicator
            is_grazing=False,
            tic_id="TIC 307210830",  # Known EB system
        )
        
        # Verify overall verdict is FALSE_POSITIVE
        assert result["overall_verdict"] == "FALSE_POSITIVE", \
            "High SNR EB should be rejected as FALSE_POSITIVE"
        
        # Verify both conflicts are detected
        conflicts = result["conflicts"]
        assert len(conflicts) == 2, "Should detect both SNR_DENSITY and SNR_SHAPE conflicts"
        
        conflict_types = {c["conflict_type"] for c in conflicts}
        assert "SNR_DENSITY_CONFLICT" in conflict_types
        assert "SNR_SHAPE_CONFLICT" in conflict_types
        
        # Verify conflict metadata
        for conflict in conflicts:
            assert conflict["tic_id"] == "TIC 307210830"
            assert conflict["snr"] == 25.0
            assert conflict["overall_verdict"] == "FALSE_POSITIVE"
            assert conflict["resolution_latency_ms"] >= 0
            
        # Verify FP-4 and FP-7 both failed
        fp4_test = next(t for t in result["tests"] if t["test"] == "FP-4 Shape Ratio (U vs V)")
        fp7_test = next(t for t in result["tests"] if t["test"] == "FP-7 Density Ratio")
        
        assert fp4_test["verdict"] == "FAIL", "FP-4 should fail for V-shaped transit"
        assert fp7_test["verdict"] == "FAIL", "FP-7 should fail for poor density match"
        
        # Verify SHAPE+DENSITY CONFLICT message is in logical_closure
        conflict_message = "SHAPE+DENSITY CONFLICT: V-shaped transit with density mismatch indicates EB"
        assert conflict_message in result["logical_closure"]
        
        print("\n=== High SNR EB Detection Test ===")
        print(f"TIC ID: {conflicts[0]['tic_id']}")
        print(f"SNR: {conflicts[0]['snr']}")
        print(f"Conflicts detected: {len(conflicts)}")
        for c in conflicts:
            print(f"  - {c['conflict_type']}: latency={c['resolution_latency_ms']:.2f}ms")
        print(f"Overall verdict: {result['overall_verdict']}")


def test_clean_planet_no_conflicts():
    """
    Integration test: Clean planetary signal with no conflicts.
    
    Scenario: A candidate has high SNR and passes all physical constraints.
    No conflicts should be detected.
    """
    # Mock check_external_catalogs to avoid network calls
    with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
        mock_catalog.return_value = {
            "is_multiple": False,
            "catalog_source": "SIMBAD",
            "risk_level": "LOW",
            "classification": "Single",
            "query_latency_ms": 150.0,
        }
        
        # Create ProofEngine with clean planetary parameters
        engine = ProofEngine(
            period_days=3.5,
            transit_depth=0.008,
            transit_duration_hrs=2.5,
            stellar_mass_solar=1.0,
            stellar_radius_solar=1.0,
            stellar_teff_k=5778.0,
            stellar_logg=4.44,
            planet_radius_earth=9.5,
            bls_snr=18.0,  # High SNR
            bls_fap=1e-7,
            even_odd_delta_sigma=0.5,
            shape_ratio=3.2,  # U-shaped transit (planet indicator)
        )
        
        # Call false_positive_ruling with good density match
        result = engine.false_positive_ruling(
            secondary_snr=0.5,
            centroid_sigma=0.5,
            density_ratio=1.1,  # Within [0.5, 2.0] - good match
            is_grazing=False,
            tic_id="TIC 260128333",
        )
        
        # Verify overall verdict is SOVEREIGN_PASS
        assert result["overall_verdict"] == "SOVEREIGN_PASS", \
            "Clean planet should pass as SOVEREIGN_PASS"
        
        # Verify no conflicts detected
        conflicts = result["conflicts"]
        assert len(conflicts) == 0, "Clean planet should have no conflicts"
        
        # Verify all critical tests passed
        assert result["n_critical_pass"] == result["n_critical"], \
            "All critical tests should pass for clean planet"
        
        print("\n=== Clean Planet Detection Test ===")
        print(f"TIC ID: TIC 260128333")
        print(f"SNR: {engine.snr}")
        print(f"Shape ratio: {engine.shape_r} (U-shape)")
        print(f"Density ratio: 1.1 (good match)")
        print(f"Conflicts detected: {len(conflicts)}")
        print(f"Overall verdict: {result['overall_verdict']}")


def test_moderate_snr_eb_no_conflict_logging():
    """
    Integration test: EB with moderate SNR (< 10.0).
    
    Scenario: A candidate has EB indicators but SNR < 10.0, so no conflicts
    are logged (conflicts only logged for high SNR cases).
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
        
        # Create ProofEngine with EB parameters but moderate SNR
        engine = ProofEngine(
            period_days=2.0,
            transit_depth=0.015,
            transit_duration_hrs=1.8,
            stellar_mass_solar=1.1,
            stellar_radius_solar=1.05,
            stellar_teff_k=5900.0,
            stellar_logg=4.4,
            planet_radius_earth=11.0,
            bls_snr=8.5,  # Moderate SNR < 10.0
            bls_fap=1e-5,
            even_odd_delta_sigma=0.9,
            shape_ratio=0.3,  # Deep V-shaped transit: fails FP-4 (floor 0.4)
        )
        
        # Call false_positive_ruling with poor density match
        result = engine.false_positive_ruling(
            secondary_snr=0.5,
            centroid_sigma=0.5,
            density_ratio=6.0,  # Outside calibrated [0.2, 5.0] band
            is_grazing=False,
            tic_id="TIC 999888777",
        )
        
        # Verify overall verdict is FALSE_POSITIVE (due to test failures)
        assert result["overall_verdict"] == "FALSE_POSITIVE", \
            "EB should be rejected as FALSE_POSITIVE"
        
        # Verify NO conflicts detected (SNR < 10.0)
        conflicts = result["conflicts"]
        assert len(conflicts) == 0, \
            "No conflicts should be logged when SNR < 10.0"
        
        # Verify FP-4 and FP-7 still failed (tests still work)
        fp4_test = next(t for t in result["tests"] if t["test"] == "FP-4 Shape Ratio (U vs V)")
        fp7_test = next(t for t in result["tests"] if t["test"] == "FP-7 Density Ratio")
        
        assert fp4_test["verdict"] == "FAIL"
        assert fp7_test["verdict"] == "FAIL"
        
        print("\n=== Moderate SNR EB Test ===")
        print(f"TIC ID: TIC 999888777")
        print(f"SNR: {engine.snr} (< 10.0, no conflict logging)")
        print(f"Shape ratio: {engine.shape_r} (V-shape)")
        print(f"Density ratio: 2.8 (poor match)")
        print(f"Conflicts detected: {len(conflicts)} (expected 0)")
        print(f"Overall verdict: {result['overall_verdict']}")
