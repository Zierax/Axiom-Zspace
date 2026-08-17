"""
Unit tests for check_external_catalogs() function.

Tests TIC ID validation, catalog cross-matching, and error handling.
"""

import pytest
from unittest.mock import patch, MagicMock
from zspace_engine.validator import check_external_catalogs


class TestCheckExternalCatalogs:
    """Test suite for check_external_catalogs function."""
    
    def test_invalid_tic_id_format(self):
        """Test that invalid TIC ID format returns appropriate error response."""
        result = check_external_catalogs("INVALID_ID")
        
        assert result["is_multiple"] is False
        assert result["catalog_source"] == "None"
        assert result["risk_level"] == "UNKNOWN"
        assert result["classification"] == "Invalid TIC ID"
        assert result["query_latency_ms"] >= 0
    
    def test_valid_tic_id_formats(self):
        """Test that various valid TIC ID formats are accepted."""
        # Clear cache before test
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        valid_formats = [
            "TIC 307210830",
            "TIC307210830",
            "307210830",
            "tic 307210830",
        ]
        
        for tic_id in valid_formats:
            # Mock SIMBAD to avoid actual network calls
            with patch('astroquery.simbad.Simbad') as mock_simbad:
                mock_instance = MagicMock()
                mock_simbad.return_value = mock_instance
                mock_instance.query_object.return_value = None
                
                # Mock Gaia to avoid actual network calls
                with patch('astroquery.gaia.Gaia') as mock_gaia:
                    mock_gaia.launch_job_async.side_effect = Exception("Network error")
                    
                    result = check_external_catalogs(tic_id, timeout=1)
                    
                    # Should return offline status when both queries fail
                    assert result["catalog_source"] == "Offline"
                    assert result["risk_level"] == "UNKNOWN"
    
    @patch('astroquery.simbad.Simbad')
    def test_simbad_multiplicity_detection(self, mock_simbad_class):
        """Test that SIMBAD multiplicity keywords are detected correctly."""
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        # Mock SIMBAD response with eclipsing binary classification
        mock_instance = MagicMock()
        mock_simbad_class.return_value = mock_instance
        
        mock_table = MagicMock()
        mock_table.__len__.return_value = 1
        mock_table.colnames = ['OTYPE', 'MAIN_ID']
        mock_table.__getitem__.side_effect = lambda key: {
            'OTYPE': ['EB*'],
            'MAIN_ID': ['TIC 307210830']
        }[key]
        
        mock_instance.query_object.return_value = mock_table
        
        result = check_external_catalogs("TIC 307210830", timeout=10)
        
        assert result["is_multiple"] is True
        assert result["catalog_source"] == "SIMBAD"
        assert result["risk_level"] == "HIGH"
        assert "EB*" in result["classification"]
    
    @patch('astroquery.simbad.Simbad')
    def test_simbad_single_star(self, mock_simbad_class):
        """Test that single stars are correctly identified."""
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        # Mock SIMBAD response with single star classification
        mock_instance = MagicMock()
        mock_simbad_class.return_value = mock_instance
        
        mock_table = MagicMock()
        mock_table.__len__.return_value = 1
        mock_table.colnames = ['OTYPE', 'MAIN_ID']
        mock_table.__getitem__.side_effect = lambda key: {
            'OTYPE': ['Star'],
            'MAIN_ID': ['TIC 123456']
        }[key]
        
        mock_instance.query_object.return_value = mock_table
        
        result = check_external_catalogs("TIC 123456", timeout=10)
        
        assert result["is_multiple"] is False
        assert result["catalog_source"] == "SIMBAD"
        assert result["risk_level"] == "LOW"
    
    @patch('astroquery.gaia.Gaia')
    @patch('astroquery.simbad.Simbad')
    def test_gaia_fallback(self, mock_simbad_class, mock_gaia_class):
        """Test that Gaia is queried when SIMBAD fails."""
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        # Mock SIMBAD failure
        mock_simbad_instance = MagicMock()
        mock_simbad_class.return_value = mock_simbad_instance
        mock_simbad_instance.query_object.side_effect = Exception("SIMBAD timeout")
        
        # Mock Gaia success with non_single_star flag
        mock_job = MagicMock()
        mock_table = MagicMock()
        mock_table.__len__.return_value = 1
        mock_table.colnames = ['non_single_star', 'phot_variable_flag']
        mock_table.__getitem__.side_effect = lambda key: {
            'non_single_star': [1],
            'phot_variable_flag': ['VARIABLE']
        }[key]
        
        mock_job.get_results.return_value = mock_table
        mock_gaia_class.launch_job_async.return_value = mock_job
        
        result = check_external_catalogs("TIC 789012", timeout=10)
        
        assert result["is_multiple"] is True
        assert result["catalog_source"] == "Gaia"
        assert result["risk_level"] == "HIGH"
        assert "non_single_star=1" in result["classification"]
    
    def test_caching(self):
        """Test that results are cached to avoid redundant queries."""
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        with patch('astroquery.simbad.Simbad') as mock_simbad_class:
            mock_instance = MagicMock()
            mock_simbad_class.return_value = mock_instance
            
            mock_table = MagicMock()
            mock_table.__len__.return_value = 1
            mock_table.colnames = ['OTYPE', 'MAIN_ID']
            mock_table.__getitem__.side_effect = lambda key: {
                'OTYPE': ['Star'],
                'MAIN_ID': ['TIC 999999']
            }[key]
            
            mock_instance.query_object.return_value = mock_table
            
            # First call - should query SIMBAD
            result1 = check_external_catalogs("TIC 999999", timeout=10)
            assert mock_instance.query_object.call_count == 1
            
            # Second call - should use cache
            result2 = check_external_catalogs("TIC 999999", timeout=10)
            assert mock_instance.query_object.call_count == 1  # No additional call
            
            # Results should be identical (except latency)
            assert result1["is_multiple"] == result2["is_multiple"]
            assert result1["catalog_source"] == result2["catalog_source"]
            assert result1["risk_level"] == result2["risk_level"]
            assert result2["query_latency_ms"] == 0.0  # Cache hit
    
    @patch('astroquery.gaia.Gaia')
    @patch('astroquery.simbad.Simbad')
    def test_network_failure_graceful_degradation(self, mock_simbad_class, mock_gaia_class):
        """Test graceful handling when both SIMBAD and Gaia are unavailable."""
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        # Mock both services failing
        mock_simbad_instance = MagicMock()
        mock_simbad_class.return_value = mock_simbad_instance
        mock_simbad_instance.query_object.side_effect = Exception("Network error")
        
        mock_gaia_class.launch_job_async.side_effect = Exception("Network error")
        
        result = check_external_catalogs("TIC 111111", timeout=1)
        
        assert result["is_multiple"] is False
        assert result["catalog_source"] == "Offline"
        assert result["risk_level"] == "UNKNOWN"
        assert result["classification"] == "Network unavailable"
        assert result["query_latency_ms"] >= 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
