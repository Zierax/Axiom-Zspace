"""
Demo test for Task 3.2: FP-9 Catalog Multiplicity Integration

This test demonstrates the complete integration of the FP-9 test
within the false_positive_ruling() method.
"""

from unittest.mock import patch
from zspace_engine.validator import ProofEngine


def test_fp9_integration_demo():
    """
    Demonstrates the FP-9 Catalog Multiplicity test integration.
    
    This test shows:
    1. check_external_catalogs() is called with the TIC ID
    2. FP-9 test is created with critical weight
    3. Verdict is set based on is_multiple flag
    4. CATALOG_OFFLINE flag is handled correctly
    """
    # Create a ProofEngine instance
    engine = ProofEngine(
        period_days=3.5,
        transit_depth=0.01,
        transit_duration_hrs=2.0,
        stellar_mass_solar=1.0,
        stellar_radius_solar=1.0,
        stellar_teff_k=5778.0,
        stellar_logg=4.44,
        planet_radius_earth=10.0,
        bls_snr=15.0,
        bls_fap=1e-6,
        even_odd_delta_sigma=1.5,
        shape_ratio=2.0,
    )
    
    # Test Case 1: Known multiple star system (should FAIL)
    print("\n=== Test Case 1: Known Multiple Star System ===")
    with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
        mock_catalog.return_value = {
            "is_multiple": True,
            "catalog_source": "SIMBAD",
            "risk_level": "HIGH",
            "classification": "EB*",
            "query_latency_ms": 150.0,
        }
        
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=1.0,
            is_grazing=False,
            tic_id="TIC 307210830",
        )
        
        fp9_test = next(t for t in result["tests"] if t["test"] == "FP-9 Catalog Multiplicity")
        print(f"FP-9 Test: {fp9_test['test']}")
        print(f"Weight: {fp9_test['weight']}")
        print(f"Verdict: {fp9_test['verdict']}")
        print(f"Description: {fp9_test['description']}")
        print(f"Comparison: {fp9_test['comparison']}")
        
        assert fp9_test["verdict"] == "FAIL"
        assert fp9_test["weight"] == "critical"
        assert result["overall_verdict"] == "FALSE_POSITIVE"  # Critical test failed
    
    # Test Case 2: Single star (should PASS)
    print("\n=== Test Case 2: Single Star ===")
    with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
        mock_catalog.return_value = {
            "is_multiple": False,
            "catalog_source": "SIMBAD",
            "risk_level": "LOW",
            "classification": "Star",
            "query_latency_ms": 120.0,
        }
        
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=1.0,
            is_grazing=False,
            tic_id="TIC 123456",
        )
        
        fp9_test = next(t for t in result["tests"] if t["test"] == "FP-9 Catalog Multiplicity")
        print(f"FP-9 Test: {fp9_test['test']}")
        print(f"Weight: {fp9_test['weight']}")
        print(f"Verdict: {fp9_test['verdict']}")
        print(f"Description: {fp9_test['description']}")
        print(f"Comparison: {fp9_test['comparison']}")
        
        assert fp9_test["verdict"] == "PASS"
        assert fp9_test["weight"] == "critical"
    
    # Test Case 3: Catalog offline (should PASS with flag)
    print("\n=== Test Case 3: Catalog Offline ===")
    with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
        mock_catalog.return_value = {
            "is_multiple": False,
            "catalog_source": "Offline",
            "risk_level": "UNKNOWN",
            "classification": "Network unavailable",
            "query_latency_ms": 0.0,
        }
        
        result = engine.false_positive_ruling(
            secondary_snr=1.0,
            centroid_sigma=1.0,
            density_ratio=1.0,
            is_grazing=False,
            tic_id="TIC 789012",
        )
        
        fp9_test = next(t for t in result["tests"] if t["test"] == "FP-9 Catalog Multiplicity")
        print(f"FP-9 Test: {fp9_test['test']}")
        print(f"Weight: {fp9_test['weight']}")
        print(f"Verdict: {fp9_test['verdict']}")
        print(f"Description: {fp9_test['description']}")
        print(f"Catalog Offline Flag: {fp9_test.get('catalog_offline')}")
        
        assert fp9_test["verdict"] == "PASS"
        assert fp9_test["weight"] == "critical"
        assert fp9_test.get("catalog_offline") is True
        assert "CATALOG_OFFLINE" in fp9_test["description"]
    
    print("\n=== All Test Cases Passed ===")


if __name__ == "__main__":
    test_fp9_integration_demo()
