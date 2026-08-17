#!/usr/bin/env python3
"""
regenerate_sector_aggregates.py
Regenerates missing/zeroed discoveries.json and summary.json for all sectors.
Scans individual discovery files and builds aggregate manifests.
"""

import json
import os
from datetime import datetime, timezone
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent / "axiom_output"
SECTORS = [7, 36, 41, 42, 55, 67]

def collect_discoveries(sector_dir: Path):
    discoveries_dir = sector_dir / "discoveries"
    known_dir = sector_dir / "known"
    discoveries = []
    known = []
    if discoveries_dir.exists():
        for f in sorted(discoveries_dir.glob("Discovery_ZS-T-*.json")):
            try:
                with open(f) as fh:
                    data = json.load(fh)
                tic_id = data.get("tic_id", f.stem.split("-")[2] if "ZS-T-" in f.stem else "0")
                period = data.get("input_parameters", {}).get("period_days", 0)
                period = period if isinstance(period, (int, float)) else 0
                cvs = data.get("cvs_score", 0)
                cvs = cvs if isinstance(cvs, (int, float)) else 0
                verdict = data.get("cvs_verdict", "UNKNOWN")
                zid = data.get("zspace_id", f.stem)
                discoveries.append({
                    "file": f.name,
                    "zspace_id": zid,
                    "tic_id": tic_id,
                    "period_days": period,
                    "cvs": cvs,
                    "verdict": verdict,
                })
            except Exception as e:
                print(f"  Error reading {f.name}: {e}")
    if known_dir.exists():
        for f in sorted(known_dir.glob("exist_planet_ZS-T-*.json")):
            try:
                with open(f) as fh:
                    data = json.load(fh)
                tic_id = data.get("tic_id", "0")
                known.append({"file": f.name, "tic_id": tic_id})
            except Exception:
                pass
    return discoveries, known

def write_aggregates(sector_dir: Path, sector: int, discoveries, known):
    discoveries_list = [d for d in discoveries if d.get("period_days", 0) > 0]
    now_utc = datetime.now(timezone.utc).isoformat()

    discoveries_json = {
        "sector": sector,
        "total_discoveries": len(discoveries_list),
        "scan_date": now_utc,
        "total_targets_scanned": len(discoveries) + len(known),
        "elapsed_minutes": 0.0,
        "planets": []
    }
    for i, d in enumerate(discoveries_list, 1):
        cvs = d.get("cvs", 0)
        if cvs >= 0.80:
            v = "PLANET CANDIDATE"
        elif cvs >= 0.55:
            v = "LIKELY PLANET CANDIDATE"
        elif cvs >= 0.35:
            v = "AMBIGUOUS"
        else:
            v = "FALSE POSITIVE"
        discoveries_json["planets"].append({
            "#": i,
            "tic_id": d["tic_id"],
            "zspace_id": d["zspace_id"],
            "period_days": round(d["period_days"], 5),
            "cvs": round(cvs, 4),
            "verdict": v,
        })

    with open(sector_dir / "discoveries.json", "w") as f:
        json.dump(discoveries_json, f, indent=2)
    print(f"  Wrote discoveries.json ({len(discoveries_list)} entries)")

    disc_dir = sector_dir / "discoveries"
    summary_json = {
        "sector": sector,
        "timestamp_utc": now_utc,
        "total_targets": len(discoveries) + len(known),
        "processed": len(discoveries) + len(known),
        "new_discoveries": len(discoveries_list),
        "known_planets": len(known),
        "false_positives": 0,
        "failed": 0,
        "discoveries": [{
            "zspace_id": d["zspace_id"],
            "tic_id": d["tic_id"],
            "period_days": d["period_days"],
            "cvs_score": d["cvs"],
            "output_file": str(disc_dir / d["file"]),
        } for d in discoveries_list],
        "errors": [],
        "elapsed_minutes": 0.0,
        "rate_per_minute": 0,
    }
    with open(sector_dir / "summary.json", "w") as f:
        json.dump(summary_json, f, indent=2)
    print(f"  Wrote summary.json ({len(discoveries_list)} discoveries)")

def main():
    for sector in SECTORS:
        sector_dir = BASE_DIR / f"sector_{sector}"
        if not sector_dir.exists():
            print(f"Sector {sector}: directory not found, skipping")
            continue
        print(f"\nSector {sector}:")
        discoveries, known = collect_discoveries(sector_dir)
        print(f"  Discovery files: {len(discoveries)}, Known files: {len(known)}")
        write_aggregates(sector_dir, sector, discoveries, known)

if __name__ == "__main__":
    main()
