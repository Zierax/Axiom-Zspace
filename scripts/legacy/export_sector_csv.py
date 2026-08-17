#!/usr/bin/env python3
"""
export_sector_csv.py — Bulk-export sector light curves to CSV for C engine
==========================================================================
Uses lightkurve (with caching) to download all sector data once,
then exports simple time,flux CSV files for the C engine to process
at full speed with zero network overhead.

Usage:
    python export_sector_csv.py 5            # Export sector 5
    python export_sector_csv.py 5 --max 100  # Export first 100 only
"""

import sys
import os
import json
import time
import logging
import argparse
import numpy as np
from pathlib import Path

# Removed lightkurve/astroquery logger overrides because they crash astropy logging in Py3.12

def load_tic_list(sector: int) -> list:
    """Load cached TIC list for sector."""
    cache_file = Path(f".cache/sector_{sector}_tics.json")
    if cache_file.exists():
        with open(cache_file) as f:
            data = json.load(f)
        
        # Handle dict format: {"sector": N, "tic_ids": [...]}
        if isinstance(data, dict):
            if 'tic_ids' in data:
                return [str(x) for x in data['tic_ids']]
            elif 'data' in data:
                return [str(x) for x in data['data']]
        
        # Handle plain list format: ["123", "456", ...]
        if isinstance(data, list):
            result = []
            for item in data:
                if isinstance(item, dict):
                    result.append(str(item.get('tic_id', item.get('ID', ''))))
                else:
                    result.append(str(item))
            return result
    
    # Download if not cached
    print(f"Downloading TIC list for sector {sector} from MAST...")
    try:
        import lightkurve as lk
        search = lk.search_lightcurve(f"TESS Sector {sector:03d}", mission="TESS")
        if search and len(search) > 0:
            tics = list(set(str(r.target_name).replace('TIC ', '') for r in search))
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_file, 'w') as f:
                json.dump({"sector": sector, "tic_ids": tics}, f)
            return tics
    except Exception as e:
        print(f"ERROR: Could not get TIC list: {e}")
    
    return []


def export_lightcurve_csv(tic_id: str, sector: int, output_dir: Path) -> dict:
    """
    Download light curve via lightkurve and export as CSV.
    Returns metadata dict or None on failure.
    """
    import lightkurve as lk
    
    csv_path = output_dir / f"TIC_{tic_id}.csv"
    meta_path = output_dir / f"TIC_{tic_id}.meta"
    
    # Skip if already exported
    if csv_path.exists() and csv_path.stat().st_size > 100:
        # Read existing metadata
        if meta_path.exists():
            with open(meta_path) as f:
                return json.load(f)
        return {"tic_id": tic_id, "status": "cached", "n_points": -1}
    
    try:
        # Search and download
        search = lk.search_lightcurve(f"TIC {tic_id}", mission="TESS", author="SPOC")
        if search is None or len(search) == 0:
            search = lk.search_lightcurve(f"TIC {tic_id}", mission="TESS")
        
        if search is None or len(search) == 0:
            return None
        
        # Download first available
        lc = search[0].download()
        if lc is None:
            return None
        
        # Extract arrays
        t = np.asarray(lc.time.value, dtype=np.float64)
        
        # Prefer PDCSAP_FLUX
        if hasattr(lc, 'pdcsap_flux') and lc.pdcsap_flux is not None:
            f = np.asarray(lc.pdcsap_flux.value, dtype=np.float64)
        else:
            f = np.asarray(lc.flux.value, dtype=np.float64)
        
        # Quality mask
        if hasattr(lc, 'quality') and lc.quality is not None:
            q = np.asarray(lc.quality.value, dtype=int)
            mask = (q == 0) & np.isfinite(f) & np.isfinite(t) & (f > 0)
        else:
            mask = np.isfinite(f) & np.isfinite(t) & (f > 0)
        
        t = t[mask]
        f = f[mask]
        
        if len(t) < 100:
            return None
        
        # Get stellar params from FITS header
        meta_info = {
            "tic_id": tic_id,
            "status": "exported",
            "n_points": len(t),
            "sector": getattr(lc.meta, 'SECTOR', sector) or sector,
            "stellar_mass": float(getattr(lc.meta, 'STELLARM', 0) or 0),
            "stellar_radius": float(getattr(lc.meta, 'STELLARR', 0) or 0),
            "stellar_teff": float(getattr(lc.meta, 'TEFF', 0) or 0),
            "stellar_logg": float(getattr(lc.meta, 'LOGG', 0) or 0),
        }
        
        # Write CSV (compact, no header — just time,flux)
        with open(csv_path, 'w') as fout:
            for i in range(len(t)):
                fout.write(f"{t[i]:.8f},{f[i]:.6f}\n")
        
        # Write metadata
        with open(meta_path, 'w') as fout:
            json.dump(meta_info, fout)
        
        return meta_info
        
    except Exception as e:
        print(f"\nError for TIC {tic_id}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Export sector light curves to CSV for C engine")
    parser.add_argument("sector", type=int, help="TESS sector number")
    parser.add_argument("--max", type=int, default=0, help="Max targets (0=all)")
    parser.add_argument("--parallel", type=int, default=4, help="Parallel downloads")
    args = parser.parse_args()
    
    sector = args.sector
    output_dir = Path(f".cache/sector_{sector}_csv")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"══════════════════════════════════════════════════════")
    print(f"  BULK EXPORT: Sector {sector} → CSV for C Engine")
    print(f"══════════════════════════════════════════════════════")
    
    # Load TIC list
    tic_list = load_tic_list(sector)
    if not tic_list:
        print("ERROR: No TIC IDs found")
        sys.exit(1)
    
    if args.max > 0:
        tic_list = tic_list[:args.max]
    
    total = len(tic_list)
    print(f"  Targets: {total}")
    print(f"  Output:  {output_dir}/")
    print()
    
    # Check how many are already cached
    cached = sum(1 for t in tic_list if (output_dir / f"TIC_{t}.csv").exists())
    if cached > 0:
        print(f"  Already cached: {cached}/{total} ({cached*100/total:.0f}%)")
    
    exported = 0
    failed = 0
    skipped = 0
    t_start = time.time()
    
    for idx, tic_id in enumerate(tic_list, 1):
        elapsed = time.time() - t_start
        rate = idx / max(elapsed, 0.01) * 60
        eta = (total - idx) / max(idx / max(elapsed, 0.01), 0.001) / 60
        
        pct = idx / total * 100
        bar_w = 30
        filled = int(bar_w * idx / total)
        bar = '█' * filled + '░' * (bar_w - filled)
        
        sys.stderr.write(
            f"\r  [{bar}] {pct:5.1f}%  ({idx}/{total})  "
            f"TIC {tic_id}  | OK:{exported} F:{failed} S:{skipped}  "
            f"| {rate:.1f}/min  ETA:{eta:.0f}m   "
        )
        sys.stderr.flush()
        
        result = export_lightcurve_csv(tic_id, sector, output_dir)
        
        if result is None:
            failed += 1
        elif result.get("status") == "cached":
            skipped += 1
            exported += 1
        else:
            exported += 1
    
    sys.stderr.write("\n\n")
    
    total_time = (time.time() - t_start) / 60
    
    # Write manifest file for C engine
    manifest = {
        "sector": sector,
        "csv_dir": str(output_dir),
        "total_exported": exported,
        "failed": failed,
        "elapsed_minutes": round(total_time, 1),
        "tic_ids": [t for t in tic_list if (output_dir / f"TIC_{t}.csv").exists()]
    }
    manifest_path = output_dir / "manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"  ╔══════════════════════════════════════════════════╗")
    print(f"  ║  EXPORT COMPLETE                                ║")
    print(f"  ║  Exported: {exported:>5}  |  Failed: {failed:>5}             ║")
    print(f"  ║  Time: {total_time:.1f} min                               ║")
    print(f"  ╚══════════════════════════════════════════════════╝")
    print(f"\n  Now run the C engine:")
    print(f"    ./scan_sector.sh {sector}")


if __name__ == "__main__":
    main()
