"""
Unit Test for Task 7.2: Conflict Logging with Resource Logger

This test verifies that conflict detection events are logged at INFO level
using the resource logger from logging_config, and that logging failures
do not cause validation to fail.
"""

import logging
from unittest.mock import patch, MagicMock
from zspace_engine.validator import ProofEngine


def test_conflict_logging_snr_density():
    """
    Test that SNR_DENSITY_CONFLICT is logged at INFO level.
    
    Verifies:
    - Logger is obtained from logging_config.get_logger(__name__)
    - Log message format: "Conflict detected: SNR_DENSITY_CONFLICT for {tic_id}"
    - Logging is at INFO level
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
        
        # Mock the logger to capture log calls
        with patch('zspace_engine.validator.get_logger') as mock_get_logger:
            mock_logger = MagicMock()
            mock_get_logger.return_value = mock_logger
            
            # Create ProofEngine with high SNR and poor density match
            engine = ProofEngine(
                period_days=1.5,
                transit_depth=0.02,
                transit_duration_hrs=1.5,
                stellar_mass_solar=1.2,
                stellar_radius_solar=1.1,
                stellar_teff_k=6000.0,
                stellar_logg=4.3,
                planet_radius_earth=12.0,
                bls_snr=15.0,  # High SNR > 10.0
                bls_fap=1e-8,
                even_odd_delta_sigma=0.8,
                shape_ratio=2.5,  # U-shaped (no shape conflict)
            )
            
            # Call false_positive_ruling with poor density match
            result = engine.false_positive_ruling(
                secondary_snr=0.5,
                centroid_sigma=0.5,
                density_ratio=3.5,  # Outside [0.5, 2.0] - triggers conflict
                is_grazing=False,
                tic_id="TIC 123456789",
            )
            
            # Verify logger was obtained with correct module name
            mock_get_logger.assert_called_once_with('zspace_engine.validator')
            
            # Verify INFO level log was called with correct message
            mock_logger.info.assert_called_once_with(
                "Conflict detected: SNR_DENSITY_CONFLICT for TIC 123456789"
            )
            
            # Verify conflict was added to result
            assert len(result["conflicts"]) == 1
            assert result["conflicts"][0]["conflict_type"] == "SNR_DENSITY_CONFLICT"
            
            print("\n=== SNR_DENSITY_CONFLICT Logging Test ===")
            print(f"Logger obtained: {mock_get_logger.called}")
            print(f"Log message: {mock_logger.info.call_args[0][0]}")
            print(f"Conflict detected: {result['conflicts'][0]['conflict_type']}")


def test_conflict_logging_snr_shape():
    """
    Test that SNR_SHAPE_CONFLICT is logged at INFO level.
    
    Verifies:
    - Logger is obtained from logging_config.get_logger(__name__)
    - Log message format: "Conflict detected: SNR_SHAPE_CONFLICT for {tic_id}"
    - Logging is at INFO level
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
        
        # Mock the logger to capture log calls
        with patch('zspace_engine.validator.get_logger') as mock_get_logger:
            mock_logger = MagicMock()
            mock_get_logger.return_value = mock_logger
            
            # Create ProofEngine with high SNR and V-shaped transit
            engine = ProofEngine(
                period_days=1.5,
                transit_depth=0.02,
                transit_duration_hrs=1.5,
                stellar_mass_solar=1.2,
                stellar_radius_solar=1.1,
                stellar_teff_k=6000.0,
                stellar_logg=4.3,
                planet_radius_earth=12.0,
                bls_snr=20.0,  # High SNR > 10.0
                bls_fap=1e-8,
                even_odd_delta_sigma=0.8,
                shape_ratio=0.8,  # V-shaped transit (< 1.0)
            )
            
            # Call false_positive_ruling with good density match
            result = engine.false_positive_ruling(
                secondary_snr=0.5,
                centroid_sigma=0.5,
                density_ratio=1.2,  # Within [0.5, 2.0] - no density conflict
                is_grazing=False,
                tic_id="TIC 987654321",
            )
            
            # Verify logger was obtained with correct module name
            mock_get_logger.assert_called_once_with('zspace_engine.validator')
            
            # Verify INFO level log was called with correct message
            mock_logger.info.assert_called_once_with(
                "Conflict detected: SNR_SHAPE_CONFLICT for TIC 987654321"
            )
            
            # Verify conflict was added to result
            assert len(result["conflicts"]) == 1
            assert result["conflicts"][0]["conflict_type"] == "SNR_SHAPE_CONFLICT"
            
            print("\n=== SNR_SHAPE_CONFLICT Logging Test ===")
            print(f"Logger obtained: {mock_get_logger.called}")
            print(f"Log message: {mock_logger.info.call_args[0][0]}")
            print(f"Conflict detected: {result['conflicts'][0]['conflict_type']}")


def test_conflict_logging_both_conflicts():
    """
    Test that both conflicts are logged when both conditions are met.
    
    Verifies:
    - Both SNR_DENSITY_CONFLICT and SNR_SHAPE_CONFLICT are logged
    - Each has its own log message
    - Both conflicts are added to result
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
        
        # Mock the logger to capture log calls
        with patch('zspace_engine.validator.get_logger') as mock_get_logger:
            mock_logger = MagicMock()
            mock_get_logger.return_value = mock_logger
            
            # Create ProofEngine with high SNR, V-shaped transit, and poor density
            engine = ProofEngine(
                period_days=1.5,
                transit_depth=0.02,
                transit_duration_hrs=1.5,
                stellar_mass_solar=1.2,
                stellar_radius_solar=1.1,
                stellar_teff_k=6000.0,
                stellar_logg=4.3,
                planet_radius_earth=12.0,
                bls_snr=25.0,  # High SNR > 10.0
                bls_fap=1e-8,
                even_odd_delta_sigma=0.8,
                shape_ratio=0.7,  # V-shaped transit
            )
            
            # Call false_positive_ruling with poor density match
            result = engine.false_positive_ruling(
                secondary_snr=0.5,
                centroid_sigma=0.5,
                density_ratio=3.5,  # Outside [0.5, 2.0]
                is_grazing=False,
                tic_id="TIC 111222333",
            )
            
            # Verify logger was obtained
            mock_get_logger.assert_called_once_with('zspace_engine.validator')
            
            # Verify both log messages were called
            assert mock_logger.info.call_count == 2
            log_messages = [call[0][0] for call in mock_logger.info.call_args_list]
            
            assert "Conflict detected: SNR_DENSITY_CONFLICT for TIC 111222333" in log_messages
            assert "Conflict detected: SNR_SHAPE_CONFLICT for TIC 111222333" in log_messages
            
            # Verify both conflicts were added to result
            assert len(result["conflicts"]) == 2
            conflict_types = {c["conflict_type"] for c in result["conflicts"]}
            assert "SNR_DENSITY_CONFLICT" in conflict_types
            assert "SNR_SHAPE_CONFLICT" in conflict_types
            
            print("\n=== Both Conflicts Logging Test ===")
            print(f"Logger obtained: {mock_get_logger.called}")
            print(f"Number of log calls: {mock_logger.info.call_count}")
            print(f"Log messages:")
            for msg in log_messages:
                print(f"  - {msg}")
            print(f"Conflicts detected: {len(result['conflicts'])}")


def test_logging_failure_does_not_break_validation():
    """
    Test that logging failures do not cause validation to fail.
    
    Verifies:
    - If logger.info() raises an exception, validation continues
    - Conflict is still added to result
    - Overall verdict is still computed correctly
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
        
        # Mock the logger to raise an exception
        with patch('zspace_engine.validator.get_logger') as mock_get_logger:
            mock_logger = MagicMock()
            mock_logger.info.side_effect = Exception("Logging system failure")
            mock_get_logger.return_value = mock_logger
            
            # Create ProofEngine with high SNR and poor density match
            engine = ProofEngine(
                period_days=1.5,
                transit_depth=0.02,
                transit_duration_hrs=1.5,
                stellar_mass_solar=1.2,
                stellar_radius_solar=1.1,
                stellar_teff_k=6000.0,
                stellar_logg=4.3,
                planet_radius_earth=12.0,
                bls_snr=15.0,  # High SNR > 10.0
                bls_fap=1e-8,
                even_odd_delta_sigma=0.8,
                shape_ratio=2.5,
            )
            
            # Call false_positive_ruling - should NOT raise exception
            try:
                result = engine.false_positive_ruling(
                    secondary_snr=0.5,
                    centroid_sigma=0.5,
                    density_ratio=6.0,  # Outside calibrated [0.2, 5.0] band
                    is_grazing=False,
                    tic_id="TIC 444555666",
                )
                validation_succeeded = True
            except Exception as e:
                validation_succeeded = False
                print(f"Validation failed with exception: {e}")
            
            # Verify validation succeeded despite logging failure
            assert validation_succeeded, "Validation should succeed even if logging fails"
            
            # Verify conflict was still added to result
            assert len(result["conflicts"]) == 1
            assert result["conflicts"][0]["conflict_type"] == "SNR_DENSITY_CONFLICT"
            
            # Verify overall verdict was still computed
            assert result["overall_verdict"] == "FALSE_POSITIVE"
            
            print("\n=== Logging Failure Resilience Test ===")
            print(f"Validation succeeded: {validation_succeeded}")
            print(f"Conflict detected: {result['conflicts'][0]['conflict_type']}")
            print(f"Overall verdict: {result['overall_verdict']}")
            print("✓ Validation continues even when logging fails")


def test_no_logging_when_snr_below_threshold():
    """
    Test that no logging occurs when SNR < 10.0.
    
    Verifies:
    - Logger is not called when SNR < 10.0
    - No conflicts are added to result
    - Validation still works correctly
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
        
        # Mock the logger to capture log calls
        with patch('zspace_engine.validator.get_logger') as mock_get_logger:
            mock_logger = MagicMock()
            mock_get_logger.return_value = mock_logger
            
            # Create ProofEngine with low SNR and poor density match
            engine = ProofEngine(
                period_days=1.5,
                transit_depth=0.02,
                transit_duration_hrs=1.5,
                stellar_mass_solar=1.2,
                stellar_radius_solar=1.1,
                stellar_teff_k=6000.0,
                stellar_logg=4.3,
                planet_radius_earth=12.0,
                bls_snr=8.0,  # Low SNR < 10.0
                bls_fap=1e-5,
                even_odd_delta_sigma=0.8,
                shape_ratio=0.7,  # V-shaped transit
            )
            
            # Call false_positive_ruling with poor density match
            result = engine.false_positive_ruling(
                secondary_snr=0.5,
                centroid_sigma=0.5,
                density_ratio=6.0,  # Outside calibrated [0.2, 5.0] band
                is_grazing=False,
                tic_id="TIC 777888999",
            )
            
            # Verify logger was NOT called (SNR < 10.0)
            mock_logger.info.assert_not_called()
            
            # Verify no conflicts were added to result
            assert len(result["conflicts"]) == 0
            
            # Verify validation still works
            assert result["overall_verdict"] == "FALSE_POSITIVE"
            
            print("\n=== No Logging Below SNR Threshold Test ===")
            print(f"SNR: {engine.snr} (< 10.0)")
            print(f"Logger called: {mock_logger.info.called}")
            print(f"Conflicts detected: {len(result['conflicts'])}")
            print(f"Overall verdict: {result['overall_verdict']}")
            print("✓ No logging when SNR < 10.0")
