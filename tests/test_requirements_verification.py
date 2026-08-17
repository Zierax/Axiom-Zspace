"""
Requirements verification tests for check_external_catalogs() function.

Validates that the implementation meets all requirements from the design document.
"""

import pytest
from unittest.mock import patch, MagicMock
from zspace_engine.validator import check_external_catalogs


class TestRequirementsVerification:
    """Verify implementation meets all design requirements."""
    
    def test_requirement_4_1_function_exists(self):
        """
        Requirement 4.1: THE Validator SHALL implement a function 
        check_external_catalogs(tic_id) that accepts a TIC identifier.
        """
        # Verify function exists and is callable
        assert callable(check_external_catalogs)
        
        # Verify function accepts tic_id parameter
        import inspect
        sig = inspect.signature(check_external_catalogs)
        assert 'tic_id' in sig.parameters
    
    def test_requirement_4_2_return_structure(self):
        """
        Requirement 4.2: WHEN check_external_catalogs is called, 
        THE Catalog_Matcher SHALL return a dictionary containing keys 
        "is_multiple", "catalog_source", and "risk_level".
        """
        result = check_external_catalogs("INVALID")
        
        # Verify all required keys are present
        assert "is_multiple" in result
        assert "catalog_source" in result
        assert "risk_level" in result
        
        # Additional keys from design spec
        assert "classification" in result
        assert "query_latency_ms" in result
    
    @patch('astroquery.simbad.Simbad')
    def test_requirement_4_3_multiplicity_detection(self, mock_simbad_class):
        """
        Requirement 4.3: WHEN THE Catalog_Matcher identifies a TIC ID 
        associated with stellar classification containing "Double" or "Multiple", 
        THE Catalog_Matcher SHALL set "is_multiple" to True.
        """
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        # Test each multiplicity keyword
        multiplicity_keywords = ['Double', 'Multiple', 'EB*', 'V*', 'Binary']
        
        for keyword in multiplicity_keywords:
            # Clear cache for each test
            check_external_catalogs._cache.clear()
            
            mock_instance = MagicMock()
            mock_simbad_class.return_value = mock_instance
            
            mock_table = MagicMock()
            mock_table.__len__.return_value = 1
            mock_table.colnames = ['OTYPE', 'MAIN_ID']
            mock_table.__getitem__.side_effect = lambda key: {
                'OTYPE': [keyword],
                'MAIN_ID': ['TIC 123456']
            }[key]
            
            mock_instance.query_object.return_value = mock_table
            
            result = check_external_catalogs("TIC 123456", timeout=10)
            
            assert result["is_multiple"] is True, f"Failed for keyword: {keyword}"
    
    @patch('astroquery.simbad.Simbad')
    def test_requirement_4_4_high_risk_level(self, mock_simbad_class):
        """
        Requirement 4.4: WHEN THE Catalog_Matcher sets "is_multiple" to True, 
        THE Catalog_Matcher SHALL set "risk_level" to "HIGH".
        """
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        mock_instance = MagicMock()
        mock_simbad_class.return_value = mock_instance
        
        mock_table = MagicMock()
        mock_table.__len__.return_value = 1
        mock_table.colnames = ['OTYPE', 'MAIN_ID']
        mock_table.__getitem__.side_effect = lambda key: {
            'OTYPE': ['EB*'],
            'MAIN_ID': ['TIC 123456']
        }[key]
        
        mock_instance.query_object.return_value = mock_table
        
        result = check_external_catalogs("TIC 123456", timeout=10)
        
        assert result["is_multiple"] is True
        assert result["risk_level"] == "HIGH"
    
    @patch('astroquery.gaia.Gaia')
    @patch('astroquery.simbad.Simbad')
    def test_requirement_4_7_network_unavailable(self, mock_simbad_class, mock_gaia_class):
        """
        Requirement 4.7: WHERE network access is unavailable, 
        THE Catalog_Matcher SHALL return "risk_level" as "UNKNOWN" 
        and THE FP_Ruling_Engine SHALL set FP-9 verdict to "PASS" 
        with a flag "CATALOG_OFFLINE".
        """
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        # Mock both services failing
        mock_simbad_instance = MagicMock()
        mock_simbad_class.return_value = mock_simbad_instance
        mock_simbad_instance.query_object.side_effect = Exception("Network error")
        
        mock_gaia_class.launch_job_async.side_effect = Exception("Network error")
        
        result = check_external_catalogs("TIC 123456", timeout=1)
        
        assert result["risk_level"] == "UNKNOWN"
        assert result["catalog_source"] == "Offline"
    
    def test_tic_id_validation(self):
        """
        Design requirement: Implement TIC ID validation (pattern: TIC \\d+).
        """
        # Valid formats should be accepted
        valid_ids = ["TIC 123456", "TIC123456", "123456", "tic 123456"]
        
        for tic_id in valid_ids:
            result = check_external_catalogs(tic_id)
            # Should not return "Invalid TIC ID" for valid formats
            # (may return "Offline" or other status depending on network)
            assert result["classification"] != "Invalid TIC ID", f"Failed for: {tic_id}"
        
        # Invalid formats should be rejected
        invalid_ids = ["ABC123", "TIC-123", ""]
        
        for tic_id in invalid_ids:
            result = check_external_catalogs(tic_id)
            assert result["classification"] == "Invalid TIC ID", f"Failed for: {tic_id}"
            assert result["risk_level"] == "UNKNOWN"
    
    @patch('astroquery.simbad.Simbad')
    def test_graceful_error_handling(self, mock_simbad_class):
        """
        Design requirement: Implement graceful error handling for network failures.
        """
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        # Mock SIMBAD raising various exceptions
        mock_instance = MagicMock()
        mock_simbad_class.return_value = mock_instance
        mock_instance.query_object.side_effect = Exception("Timeout")
        
        # Should not raise exception, should return graceful degradation
        result = check_external_catalogs("TIC 123456", timeout=1)
        
        assert isinstance(result, dict)
        assert "is_multiple" in result
        assert result["catalog_source"] in ["Offline", "Gaia"]  # Falls back to Gaia or offline
    
    def test_in_memory_cache(self):
        """
        Design requirement: Add simple in-memory cache for query results.
        """
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
            
            # First call
            result1 = check_external_catalogs("TIC 999999", timeout=10)
            call_count_1 = mock_instance.query_object.call_count
            
            # Second call - should use cache
            result2 = check_external_catalogs("TIC 999999", timeout=10)
            call_count_2 = mock_instance.query_object.call_count
            
            # Verify cache was used (no additional call)
            assert call_count_2 == call_count_1
            
            # Verify cache hit has zero latency
            assert result2["query_latency_ms"] == 0.0
    
    def test_query_latency_measurement(self):
        """
        Design requirement: Return dict with key query_latency_ms.
        """
        result = check_external_catalogs("INVALID")
        
        assert "query_latency_ms" in result
        assert isinstance(result["query_latency_ms"], float)
        assert result["query_latency_ms"] >= 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
