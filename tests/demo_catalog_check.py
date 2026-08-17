"""
Demo script to test check_external_catalogs() function manually.

This script demonstrates the function with various TIC IDs.
"""

from zspace_engine.validator import check_external_catalogs


def demo_catalog_check():
    """Demonstrate catalog checking with various TIC IDs."""
    
    print("=" * 70)
    print("Demo: check_external_catalogs() Function")
    print("=" * 70)
    
    # Test cases
    test_cases = [
        ("TIC 307210830", "Known eclipsing binary (if in SIMBAD)"),
        ("TIC 123456789", "Random TIC ID"),
        ("307210830", "TIC ID without prefix"),
        ("INVALID_ID", "Invalid format"),
    ]
    
    for tic_id, description in test_cases:
        print(f"\nTest: {description}")
        print(f"TIC ID: {tic_id}")
        print("-" * 70)
        
        try:
            result = check_external_catalogs(tic_id, timeout=10)
            
            print(f"  is_multiple:      {result['is_multiple']}")
            print(f"  catalog_source:   {result['catalog_source']}")
            print(f"  risk_level:       {result['risk_level']}")
            print(f"  classification:   {result['classification']}")
            print(f"  query_latency_ms: {result['query_latency_ms']:.2f} ms")
            
        except Exception as e:
            print(f"  ERROR: {e}")
    
    print("\n" + "=" * 70)
    print("Demo complete!")
    print("=" * 70)


if __name__ == "__main__":
    demo_catalog_check()
