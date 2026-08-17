#!/usr/bin/env python3
"""Test fetching known planets from NASA Exoplanet Archive."""

try:
    import astroquery
    from astropy import log as astropy_log
    astropy_log.setLevel('WARNING')
except ImportError:
    pass

from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

print("Querying NASA Exoplanet Archive for TESS planets...")

try:
    table = NasaExoplanetArchive.query_criteria(
        table="pscomppars",
        select="pl_name,hostname,tic_id,pl_orbper,pl_rade,pl_trandep",
        where="tran_flag=1 AND tic_id IS NOT NULL",
        order="pl_orbper",
    )
    
    print(f"Found {len(table)} planets")
    print("\nFirst 5 planets:")
    print(table[:5])
    
    # Test accessing values
    for i, row in enumerate(table[:5]):
        print(f"\n{i+1}. {row['pl_name']}")
        print(f"   TIC ID: {row['tic_id']}")
        print(f"   Period: {row['pl_orbper']}")
        print(f"   Depth: {row['pl_trandep']}")
        
        # Check if masked
        if hasattr(row['tic_id'], 'mask'):
            print(f"   TIC masked: {row['tic_id'].mask}")
        if hasattr(row['pl_orbper'], 'mask'):
            print(f"   Period masked: {row['pl_orbper'].mask}")
        if hasattr(row['pl_trandep'], 'mask'):
            print(f"   Depth masked: {row['pl_trandep'].mask}")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
