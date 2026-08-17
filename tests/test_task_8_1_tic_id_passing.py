"""
Test Task 8.1: Verify tic_id is passed from validate() to false_positive_ruling()

This test verifies that the tic_id parameter flows correctly from the
AxiomValidator.validate() method through to the false_positive_ruling() method.
"""

import pytest
from unittest.mock import patch, MagicMock
from zspace_engine.validator import AxiomValidator, ProofEngine


def test_tic_id_passed_to_false_positive_ruling():
    """
    Verify that tic_id from validate() is passed to false_positive_ruling().
    
    This test mocks the false_positive_ruling method to capture the tic_id
    parameter and verifies it matches the input to validate().
    """
    # Create validator instance
    validator = AxiomValidator(output_dir=".")
    
    # Test TIC ID
    test_tic_id = "307210830"
    
    # Track if false_positive_ruling was called with correct tic_id
    captured_tic_id = None
    
    # Mock the false_positive_ruling method to capture tic_id
    original_fp_ruling = ProofEngine.false_positive_ruling
    
    def mock_fp_ruling(self, secondary_snr, centroid_sigma, density_ratio,
                       is_grazing, tic_id, n_transits=None, **kwargs):
        nonlocal captured_tic_id
        captured_tic_id = tic_id
        # Call the original method
        return original_fp_ruling(self, secondary_snr, centroid_sigma, density_ratio, is_grazing, tic_id)
    
    with patch.object(ProofEngine, 'false_positive_ruling', mock_fp_ruling):
        # Mock network queries to avoid actual API calls
        with patch.object(validator._querier, 'query', return_value=([], [], None)):
            # Call validate with test_tic_id
            result = validator.validate(
                tic_id="307210830",
                period_days=3.6986,
                transit_depth=0.00836,
                transit_duration_hrs=2.0,
                t0_btjd=1201.0,
                stellar_mass_solar=1.0,
                stellar_radius_solar=1.0,
                stellar_teff_k=5778.0,
                stellar_logg=4.44,
                planet_radius_earth=9.97,
                cvs_score=0.83,
                cvs_verdict="PLANET CANDIDATE",
                cvs_proof_chain=[],
                bls_snr=12.5,
                bls_fap=1e-6,
                even_odd_delta_sigma=0.975,
                shape_ratio=4.711,
                secondary_snr=0.5,
                centroid_sigma=0.5,
            )
    
    # Verify that false_positive_ruling was called with the correct tic_id
    assert captured_tic_id == test_tic_id, \
        f"Expected tic_id '{test_tic_id}' but got '{captured_tic_id}'"
    
    print(f"✓ tic_id '{test_tic_id}' successfully passed to false_positive_ruling()")


def test_tic_id_used_in_catalog_check():
    """
    Verify that the tic_id passed to false_positive_ruling is used in catalog check.
    
    This test mocks check_external_catalogs to verify it receives the correct tic_id.
    """
    from zspace_engine import validator as validator_module
    
    # Test TIC ID
    test_tic_id = "260128333"
    
    # Track if check_external_catalogs was called with correct tic_id
    captured_catalog_tic_id = None
    
    def mock_catalog_check(tic_id, timeout=10):
        nonlocal captured_catalog_tic_id
        captured_catalog_tic_id = tic_id
        return {
            "is_multiple": False,
            "catalog_source": "SIMBAD",
            "risk_level": "LOW",
            "classification": "Single",
            "query_latency_ms": 100.0,
        }
    
    # Create ProofEngine instance
    engine = ProofEngine(
        period_days=5.0,
        transit_depth=0.01,
        transit_duration_hrs=2.5,
        stellar_mass_solar=1.0,
        stellar_radius_solar=1.0,
        stellar_teff_k=5778.0,
        stellar_logg=4.44,
        planet_radius_earth=10.0,
        bls_snr=12.0,
        bls_fap=1e-5,
        even_odd_delta_sigma=1.5,
        shape_ratio=3.0,
    )
    
    # Mock check_external_catalogs
    with patch.object(validator_module, 'check_external_catalogs', mock_catalog_check):
        # Call false_positive_ruling with test_tic_id
        result = engine.false_positive_ruling(
            secondary_snr=0.5,
            centroid_sigma=0.5,
            density_ratio=1.2,
            is_grazing=False,
            tic_id=test_tic_id,
        )
    
    # Verify that check_external_catalogs was called with the correct tic_id
    assert captured_catalog_tic_id == test_tic_id, \
        f"Expected catalog check with tic_id '{test_tic_id}' but got '{captured_catalog_tic_id}'"
    
    print(f"✓ tic_id '{test_tic_id}' successfully used in check_external_catalogs()")


if __name__ == "__main__":
    test_tic_id_passed_to_false_positive_ruling()
    test_tic_id_used_in_catalog_check()
    print("\n✓ All Task 8.1 tests passed!")
