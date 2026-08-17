"""
Verification test for Task 2.2: Stellar multiplicity detection logic.

This test verifies that the check_external_catalogs function correctly:
1. Parses SIMBAD OTYPE field for "Double", "Multiple", "EB*", "V*"
2. Sets is_multiple=True when multiplicity is detected
3. Sets risk_level="HIGH" when multiplicity is detected
"""

import pytest
from unittest.mock import patch, MagicMock
from zspace_engine.validator import check_external_catalogs


class TestTask22MultiplicitDetection:
    """Test suite verifying Task 2.2 requirements."""
    
    @pytest.mark.parametrize("otype_value,expected_multiple", [
        ("EB*", True),           # Eclipsing Binary
        ("V*", True),            # Variable Star
        ("Double", True),        # Double star
        ("Multiple", True),      # Multiple star system
        ("Star", False),         # Single star
        ("PM*", False),          # Proper motion star (single)
    ])
    @patch('astroquery.simbad.Simbad')
    def test_otype_multiplicity_keywords(self, mock_simbad_class, otype_value, expected_multiple):
        """
        Test that SIMBAD OTYPE field is correctly parsed for multiplicity keywords.
        
        Requirements 4.3, 4.4:
        - Parse OTYPE for "Double", "Multiple", "EB*", "V*"
        - Set is_multiple=True and risk_level="HIGH" when detected
        """
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        # Mock SIMBAD response
        mock_instance = MagicMock()
        mock_simbad_class.return_value = mock_instance
        
        mock_table = MagicMock()
        mock_table.__len__.return_value = 1
        mock_table.colnames = ['OTYPE', 'MAIN_ID']
        mock_table.__getitem__.side_effect = lambda key: {
            'OTYPE': [otype_value],
            'MAIN_ID': ['TIC 123456']
        }[key]
        
        mock_instance.query_object.return_value = mock_table
        
        # Execute
        result = check_external_catalogs("TIC 123456", timeout=10)
        
        # Verify
        assert result["is_multiple"] == expected_multiple, \
            f"OTYPE '{otype_value}' should set is_multiple={expected_multiple}"
        
        if expected_multiple:
            assert result["risk_level"] == "HIGH", \
                f"OTYPE '{otype_value}' should set risk_level='HIGH'"
            assert result["catalog_source"] == "SIMBAD"
        else:
            assert result["risk_level"] == "LOW", \
                f"OTYPE '{otype_value}' should set risk_level='LOW'"
    
    @patch('astroquery.simbad.Simbad')
    def test_case_insensitive_matching(self, mock_simbad_class):
        """Test that multiplicity keyword matching is case-insensitive."""
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        # Mock SIMBAD response with lowercase
        mock_instance = MagicMock()
        mock_simbad_class.return_value = mock_instance
        
        mock_table = MagicMock()
        mock_table.__len__.return_value = 1
        mock_table.colnames = ['OTYPE', 'MAIN_ID']
        mock_table.__getitem__.side_effect = lambda key: {
            'OTYPE': ['eb*'],  # lowercase
            'MAIN_ID': ['TIC 789012']
        }[key]
        
        mock_instance.query_object.return_value = mock_table
        
        result = check_external_catalogs("TIC 789012", timeout=10)
        
        assert result["is_multiple"] is True
        assert result["risk_level"] == "HIGH"
    
    @patch('astroquery.simbad.Simbad')
    def test_partial_keyword_matching(self, mock_simbad_class):
        """Test that partial keyword matches work (e.g., 'EB*WUMa' contains 'EB*')."""
        # Clear cache
        if hasattr(check_external_catalogs, '_cache'):
            check_external_catalogs._cache.clear()
        
        # Mock SIMBAD response with extended classification
        mock_instance = MagicMock()
        mock_simbad_class.return_value = mock_instance
        
        mock_table = MagicMock()
        mock_table.__len__.return_value = 1
        mock_table.colnames = ['OTYPE', 'MAIN_ID']
        mock_table.__getitem__.side_effect = lambda key: {
            'OTYPE': ['EB*WUMa'],  # Specific EB subtype
            'MAIN_ID': ['TIC 456789']
        }[key]
        
        mock_instance.query_object.return_value = mock_table
        
        result = check_external_catalogs("TIC 456789", timeout=10)
        
        assert result["is_multiple"] is True
        assert result["risk_level"] == "HIGH"
        assert "EB*" in result["classification"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
