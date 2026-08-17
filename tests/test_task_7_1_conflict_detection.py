"""
Test Task 7.1: Conflict Detection Logic for SNR vs Density/Shape Mismatches

This test verifies that the false_positive_ruling() method correctly detects
conflicts when:
1. SNR > 10.0 AND density_ratio outside [0.5, 2.0] → SNR_DENSITY_CONFLICT
2. SNR > 10.0 AND shape_ratio <= 1.0 → SNR_SHAPE_CONFLICT

Validates: Requirements 5.1, 5.2, 5.3, 5.4
"""

from unittest.mock import patch
from zspace_engine.validator import ProofEngine


def test_snr_density_conflict_detection():
    """
    Test that SNR > 10.0 AND density_ratio outside [0.5, 2.0] triggers
    SNR_DENSITY_CONFLICT.
    
    Validates: Requirement 5.1
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
        
        # Create ProofEngine with high SNR (> 10.0)
        engine = ProofEngine(
            period_days=3.5,
            transit_depth=0.01,
            transit_duration_hrs=2.0,
            stellar_mass_solar=1.0,
            stellar_radius_solar=1.0,
            stellar_teff_k=5778.0,
            stellar_logg=4.44,
            planet_radius_earth=10.0,
            bls_snr=15.0,  # High SNR > 10.0
            bls_fap=1e-6,
            even_odd_delta_sigma=0.5,
            shape_ratio=2.0,  # U-shape (not triggering shape conflict)
        )
        
        # Call false_positive_ruling with density_ratio outside [0.5, 2.0]
        result = engine.false_positive_ruling(
            secondary_snr=0.5,
            centroid_sigma=0.5,
            density_ratio=0.3,  # Outside [0.5, 2.0]
            is_grazing=False,
            tic_id="TIC 123456789",
        )
        
        # Verify conflicts list exists
        assert "conflicts" in result, "Result should contain 'conflicts' field"
        
        # Verify SNR_DENSITY_CONFLICT is detected
        conflicts = result["conflicts"]
        assert len(conflicts) >= 1, "Should detect at least one conflict"
        
        snr_density_conflicts = [c for c in conflicts if c["conflict_type"] == "SNR_DENSITY_CONFLICT"]
        assert len(snr_density_conflicts) == 1, "Should detect exactly one SNR_DENSITY_CONFLICT"
        
        conflict = snr_density_conflicts[0]
        assert conflict["tic_id"] == "TIC 123456789"
        assert conflict["snr"] == 15.0
        assert conflict["density_ratio"] == 0.3
        assert conflict["shape_ratio"] is None
        assert "timestamp" in conflict
        assert "resolution_latency_ms" in conflict
        assert "overall_verdict" in conflict


def test_snr_shape_conflict_detection():
    """
    Test that SNR > 10.0 AND shape_ratio <= 1.0 triggers SNR_SHAPE_CONFLICT.
    
    Validates: Requirement 5.2
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
        
        # Create ProofEngine with high SNR (> 10.0) and V-shape (shape_ratio <= 1.0)
        engine = ProofEngine(
            period_days=3.5,
            transit_depth=0.01,
            transit_duration_hrs=2.0,
            stellar_mass_solar=1.0,
            stellar_radius_solar=1.0,
            stellar_teff_k=5778.0,
            stellar_logg=4.44,
            planet_radius_earth=10.0,
            bls_snr=12.0,  # High SNR > 10.0
            bls_fap=1e-6,
            even_odd_delta_sigma=0.5,
            shape_ratio=0.8,  # V-shape <= 1.0
        )
        
        # Call false_positive_ruling with density_ratio in valid range
        result = engine.false_positive_ruling(
            secondary_snr=0.5,
            centroid_sigma=0.5,
            density_ratio=1.0,  # Within [0.5, 2.0]
            is_grazing=False,
            tic_id="TIC 987654321",
        )
        
        # Verify conflicts list exists
        assert "conflicts" in result, "Result should contain 'conflicts' field"
        
        # Verify SNR_SHAPE_CONFLICT is detected
        conflicts = result["conflicts"]
        assert len(conflicts) >= 1, "Should detect at least one conflict"
        
        snr_shape_conflicts = [c for c in conflicts if c["conflict_type"] == "SNR_SHAPE_CONFLICT"]
        assert len(snr_shape_conflicts) == 1, "Should detect exactly one SNR_SHAPE_CONFLICT"
        
        conflict = snr_shape_conflicts[0]
        assert conflict["tic_id"] == "TIC 987654321"
        assert conflict["snr"] == 12.0
        assert conflict["density_ratio"] is None
        assert conflict["shape_ratio"] == 0.8
        assert "timestamp" in conflict
        assert "resolution_latency_ms" in conflict
        assert "overall_verdict" in conflict


def test_both_conflicts_detected():
    """
    Test that both SNR_DENSITY_CONFLICT and SNR_SHAPE_CONFLICT are detected
    when both conditions are met.
    
    Validates: Requirements 5.1, 5.2
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
        
        # Create ProofEngine with high SNR, V-shape, and bad density
        engine = ProofEngine(
            period_days=3.5,
            transit_depth=0.01,
            transit_duration_hrs=2.0,
            stellar_mass_solar=1.0,
            stellar_radius_solar=1.0,
            stellar_teff_k=5778.0,
            stellar_logg=4.44,
            planet_radius_earth=10.0,
            bls_snr=20.0,  # High SNR > 10.0
            bls_fap=1e-6,
            even_odd_delta_sigma=0.5,
            shape_ratio=0.9,  # V-shape <= 1.0
        )
        
        # Call false_positive_ruling with density_ratio outside range
        result = engine.false_positive_ruling(
            secondary_snr=0.5,
            centroid_sigma=0.5,
            density_ratio=2.5,  # Outside [0.5, 2.0]
            is_grazing=False,
            tic_id="TIC 111222333",
        )
        
        # Verify conflicts list exists
        assert "conflicts" in result, "Result should contain 'conflicts' field"
        
        # Verify both conflicts are detected
        conflicts = result["conflicts"]
        assert len(conflicts) == 2, "Should detect exactly two conflicts"
        
        conflict_types = {c["conflict_type"] for c in conflicts}
        assert "SNR_DENSITY_CONFLICT" in conflict_types
        assert "SNR_SHAPE_CONFLICT" in conflict_types


def test_no_conflict_when_snr_low():
    """
    Test that no conflicts are detected when SNR <= 10.0, even if density
    and shape ratios are outside acceptable ranges.
    
    Validates: Requirements 5.1, 5.2
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
        
        # Create ProofEngine with low SNR (<= 10.0)
        engine = ProofEngine(
            period_days=3.5,
            transit_depth=0.01,
            transit_duration_hrs=2.0,
            stellar_mass_solar=1.0,
            stellar_radius_solar=1.0,
            stellar_teff_k=5778.0,
            stellar_logg=4.44,
            planet_radius_earth=10.0,
            bls_snr=8.0,  # Low SNR <= 10.0
            bls_fap=1e-6,
            even_odd_delta_sigma=0.5,
            shape_ratio=0.8,  # V-shape <= 1.0
        )
        
        # Call false_positive_ruling with density_ratio outside range
        result = engine.false_positive_ruling(
            secondary_snr=0.5,
            centroid_sigma=0.5,
            density_ratio=0.3,  # Outside [0.5, 2.0]
            is_grazing=False,
            tic_id="TIC 444555666",
        )
        
        # Verify conflicts list exists but is empty
        assert "conflicts" in result, "Result should contain 'conflicts' field"
        assert len(result["conflicts"]) == 0, "Should not detect any conflicts when SNR <= 10.0"


def test_conflict_metadata_completeness():
    """
    Test that conflict metadata contains all required keys.
    
    Validates: Requirements 5.3, 5.4
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
        
        # Create ProofEngine with high SNR and bad density
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
            even_odd_delta_sigma=0.5,
            shape_ratio=2.0,
        )
        
        # Call false_positive_ruling
        result = engine.false_positive_ruling(
            secondary_snr=0.5,
            centroid_sigma=0.5,
            density_ratio=0.4,
            is_grazing=False,
            tic_id="TIC 777888999",
        )
        
        # Verify conflict metadata completeness
        conflicts = result["conflicts"]
        assert len(conflicts) >= 1, "Should detect at least one conflict"
        
        for conflict in conflicts:
            # Required keys from requirements 5.3, 5.4
            required_keys = [
                "conflict_type",
                "timestamp",
                "tic_id",
                "snr",
                "density_ratio",
                "shape_ratio",
                "resolution_latency_ms",
                "overall_verdict",
            ]
            
            for key in required_keys:
                assert key in conflict, f"Conflict metadata should contain '{key}'"
            
            # Verify timestamp is ISO 8601 format
            assert "T" in conflict["timestamp"], "Timestamp should be in ISO 8601 format"
            assert conflict["timestamp"].endswith("Z") or "+" in conflict["timestamp"], \
                "Timestamp should include timezone information"
            
            # Verify resolution_latency_ms is a number
            assert isinstance(conflict["resolution_latency_ms"], (int, float)), \
                "resolution_latency_ms should be a number"
            assert conflict["resolution_latency_ms"] >= 0, \
                "resolution_latency_ms should be non-negative"


def test_no_conflict_when_conditions_not_met():
    """
    Test that no conflicts are detected when SNR is high but density and shape
    ratios are within acceptable ranges.
    
    Validates: Requirements 5.1, 5.2
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
        
        # Create ProofEngine with high SNR but good density and shape
        engine = ProofEngine(
            period_days=3.5,
            transit_depth=0.01,
            transit_duration_hrs=2.0,
            stellar_mass_solar=1.0,
            stellar_radius_solar=1.0,
            stellar_teff_k=5778.0,
            stellar_logg=4.44,
            planet_radius_earth=10.0,
            bls_snr=15.0,  # High SNR > 10.0
            bls_fap=1e-6,
            even_odd_delta_sigma=0.5,
            shape_ratio=2.5,  # U-shape > 1.0
        )
        
        # Call false_positive_ruling with density_ratio in valid range
        result = engine.false_positive_ruling(
            secondary_snr=0.5,
            centroid_sigma=0.5,
            density_ratio=1.2,  # Within [0.5, 2.0]
            is_grazing=False,
            tic_id="TIC 123123123",
        )
        
        # Verify conflicts list exists but is empty
        assert "conflicts" in result, "Result should contain 'conflicts' field"
        assert len(result["conflicts"]) == 0, \
            "Should not detect any conflicts when all conditions are within acceptable ranges"
