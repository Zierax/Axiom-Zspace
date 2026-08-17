"""
Integration test for check_external_catalogs() function.

Tests the function can be imported and called from the validator module.
"""

import pytest
from zspace_engine.validator import check_external_catalogs


def test_function_import():
    """Test that check_external_catalogs can be imported from validator module."""
    assert callable(check_external_catalogs)


def test_function_signature():
    """Test that function has correct signature and returns expected structure."""
    # Test with invalid ID to avoid network calls
    result = check_external_catalogs("INVALID")
    
    # Verify return structure
    assert isinstance(result, dict)
    assert "is_multiple" in result
    assert "catalog_source" in result
    assert "risk_level" in result
    assert "classification" in result
    assert "query_latency_ms" in result
    
    # Verify types
    assert isinstance(result["is_multiple"], bool)
    assert isinstance(result["catalog_source"], str)
    assert isinstance(result["risk_level"], str)
    assert isinstance(result["classification"], str)
    assert isinstance(result["query_latency_ms"], float)


def test_timeout_parameter():
    """Test that timeout parameter is accepted."""
    result = check_external_catalogs("INVALID", timeout=5)
    assert result is not None


def test_return_values():
    """Test that return values match specification."""
    result = check_external_catalogs("INVALID")
    
    # catalog_source should be one of the specified values
    assert result["catalog_source"] in ["SIMBAD", "Gaia", "None", "Offline"]
    
    # risk_level should be one of the specified values
    assert result["risk_level"] in ["HIGH", "LOW", "UNKNOWN"]
    
    # is_multiple should be boolean
    assert result["is_multiple"] in [True, False]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
