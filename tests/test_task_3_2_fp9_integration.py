"""
Test for Task 3.2: FP-9 Catalog Multiplicity test integration.

Verifies that check_external_catalogs() is called within false_positive_ruling()
and that the FP-9 test entry is created correctly.
"""

import pytest
from unittest.mock import patch, MagicMock
from zspace_engine.validator import ProofEngine


class TestFP9Integration:
    """Test suite for FP-9 Catalog Multiplicity test integration."""
    
    def test_fp9_test_present_in_tests_list(self):
        """Verify that FP-9 test is present in the tests list."""
        # Create a ProofEngine instance with minimal parameters
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
        
        # Mock check_external_catalogs to avoid network calls
        with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
            mock_catalog.return_value = {
                "is_multiple": False,
                "catalog_source": "SIMBAD",
                "risk_level": "LOW",
                "classification": "Star",
                "query_latency_ms": 100.0,
            }
            
            # Call false_positive_ruling
            result = engine.false_positive_ruling(
                secondary_snr=1.0,
                centroid_sigma=1.0,
                density_ratio=1.0,
                is_grazing=False,
                tic_id="TIC 123456",
            )
            
            # Verify FP-9 test is present
            test_names = [t["test"] for t in result["tests"]]
            assert "FP-9 Catalog Multiplicity" in test_names
            
            # Verify check_external_catalogs was called
            mock_catalog.assert_called_once_with("TIC 123456")
    
    def test_fp9_test_has_critical_weight(self):
        """Verify that FP-9 test has critical weight."""
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
        
        with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
            mock_catalog.return_value = {
                "is_multiple": False,
                "catalog_source": "SIMBAD",
                "risk_level": "LOW",
                "classification": "Star",
                "query_latency_ms": 100.0,
            }
            
            result = engine.false_positive_ruling(
                secondary_snr=1.0,
                centroid_sigma=1.0,
                density_ratio=1.0,
                is_grazing=False,
                tic_id="TIC 123456",
            )
            
            # Find FP-9 test
            fp9_test = next(t for t in result["tests"] if t["test"] == "FP-9 Catalog Multiplicity")
            
            # Verify weight is critical
            assert fp9_test["weight"] == "critical"
    
    def test_fp9_fails_when_is_multiple_true(self):
        """Verify that FP-9 test fails when is_multiple is True."""
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
        
        with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
            mock_catalog.return_value = {
                "is_multiple": True,
                "catalog_source": "SIMBAD",
                "risk_level": "HIGH",
                "classification": "EB*",
                "query_latency_ms": 100.0,
            }
            
            result = engine.false_positive_ruling(
                secondary_snr=1.0,
                centroid_sigma=1.0,
                density_ratio=1.0,
                is_grazing=False,
                tic_id="TIC 123456",
            )
            
            # Find FP-9 test
            fp9_test = next(t for t in result["tests"] if t["test"] == "FP-9 Catalog Multiplicity")
            
            # Verify verdict is FAIL
            assert fp9_test["verdict"] == "FAIL"
    
    def test_fp9_passes_when_is_multiple_false(self):
        """Verify that FP-9 test passes when is_multiple is False."""
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
        
        with patch('zspace_engine.validator.check_external_catalogs') as mock_catalog:
            mock_catalog.return_value = {
                "is_multiple": False,
                "catalog_source": "SIMBAD",
                "risk_level": "LOW",
                "classification": "Star",
                "query_latency_ms": 100.0,
            }
            
            result = engine.false_positive_ruling(
                secondary_snr=1.0,
                centroid_sigma=1.0,
                density_ratio=1.0,
                is_grazing=False,
                tic_id="TIC 123456",
            )
            
            # Find FP-9 test
            fp9_test = next(t for t in result["tests"] if t["test"] == "FP-9 Catalog Multiplicity")
            
            # Verify verdict is PASS
            assert fp9_test["verdict"] == "PASS"
    
    def test_fp9_handles_catalog_offline(self):
        """Verify that FP-9 test handles CATALOG_OFFLINE flag correctly."""
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
                tic_id="TIC 123456",
            )
            
            # Find FP-9 test
            fp9_test = next(t for t in result["tests"] if t["test"] == "FP-9 Catalog Multiplicity")
            
            # Verify verdict is PASS (benefit of the doubt when offline)
            assert fp9_test["verdict"] == "PASS"
            
            # Verify catalog_offline flag is set
            assert fp9_test.get("catalog_offline") is True
            
            # Verify description mentions CATALOG_OFFLINE
            assert "CATALOG_OFFLINE" in fp9_test["description"]