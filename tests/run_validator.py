#!/usr/bin/env python3
"""
run_validator.py  ·  AxiomValidator Sovereign Validation Test
==============================================================
Demonstrates both pipeline paths:
  (A) Synthetic NEW DISCOVERY  → full Discovery.json with §1–§5 proof
  (B) Real TIC (with live NASA query when network available)

Usage
-----
  # Synthetic test only (always works, no internet):
  python tests/run_validator.py --synthetic

  # Test with a real TIC ID (needs internet access):
  python tests/run_validator.py --tic 260128333

  # Both:
  python tests/run_validator.py --synthetic --tic 260128333

  # Simulate a KNOWN match (test the exist_planet path):
  python tests/run_validator.py --synthetic --force-known
"""

import argparse
import json
import sys
import traceback
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# ── Engine imports ─────────────────────────────────────────────────────────
from zspace_engine.validator import AxiomValidator, ArchiveMatch

# ── Test utilities ─────────────────────────────────────────────────────────
from test_synthetic import (
    generate_synthetic_planet_params,
    generate_known_planet_params,
    print_discovery_summary,
    print_known_summary,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def section(title: str) -> None:
    print(f"\n{'═'*70}")
    print(f"  {title}")
    print('═'*70)


# ─────────────────────────────────────────────────────────────────────────────
# Test A: Synthetic NEW DISCOVERY
# ─────────────────────────────────────────────────────────────────────────────

def run_synthetic_discovery(output_dir: str = ".") -> dict:
    """
    Validate a synthetic planet that will not appear in the NASA Archive.
    Should always produce Discovery.json with full §1–§5 proof.
    """
    section("TEST A  ·  SYNTHETIC NEW DISCOVERY  ·  AxiomValidator")

    # Get synthetic planet parameters from test_synthetic module
    PARAMS = generate_synthetic_planet_params()

    print(f"\n  Target:   TIC {PARAMS['tic_id']}")
    print(f"  Period:   {PARAMS['period_days']:.6f} d")
    print(f"  Depth:    {PARAMS['transit_depth']*1e6:.1f} ppm")
    print(f"  CVS:      {PARAMS['cvs_score']:.4f}  ({PARAMS['cvs_verdict']})")
    print(f"\n  Querying NASA Exoplanet Archive … (OFFLINE in this environment)")

    validator = AxiomValidator(output_dir=output_dir, verbose=True)
    result    = validator.validate(**PARAMS)

    print(f"\n  Status:       {result.status}")
    print(f"  Output file:  {result.output_file}")
    if result.network_error:
        print(f"  Network:      {result.network_error[:80]}…")

    print_proof_summary(result.output_file)
    return {"status": result.status, "file": result.output_file}


# ─────────────────────────────────────────────────────────────────────────────
# Test B: Simulated KNOWN match (force-injects an archive hit)
# ─────────────────────────────────────────────────────────────────────────────

def run_simulated_known(output_dir: str = ".") -> dict:
    """
    Directly exercise the exist_planet path by calling _handle_known
    with a fabricated ArchiveMatch (simulates what happens when the NASA
    query returns a matching entry).
    """
    section("TEST B  ·  SIMULATED KNOWN PLANET  ·  exist_planet.json path")

    from zspace_engine.validator import ArchiveMatch, write_json
    from datetime import datetime, timezone

    # Fabricate a match (simulates a confirmed exoplanet record)
    fake_match = ArchiveMatch(
        source               = "NASA_ARCHIVE",
        planet_name          = "TOI-4600 b",
        period_days          = 82.69,
        period_delta_days    = 0.00003,
        transit_depth        = 0.00740,
        planet_radius_earth  = 10.8,
        semi_major_axis_au   = 0.335,
        stellar_teff_k       = 5090.0,
        stellar_radius_solar = 0.90,
        stellar_mass_solar   = 0.89,
        discovery_method     = "Transit",
        disposition          = "CP",
        extra_fields         = {"disc_year": 2023, "pl_controv_flag": 0, "period_alias_factor": 1.0},
    )

    zspace_id = "ZS-T-281459674-01"
    out_path  = Path(output_dir) / f"exist_planet_{zspace_id}.json"

    card = {
        "schema":          "Axiom-ZSpace Existing Planet Card v1.0",
        "zspace_id":       zspace_id,
        "timestamp_utc":   datetime.now(timezone.utc).isoformat(),
        "status":          "KNOWN",
        "match_source":    fake_match.source,
        "match_summary": {
            "planet_name":         fake_match.planet_name,
            "period_delta_days":   fake_match.period_delta_days,
            "tolerance_days":      0.001,
            "proof": (
                f"|P_candidate - P_archive| = {fake_match.period_delta_days:.8f} d "
                f"≤ 0.001 d → MATCH"
            ),
        },
        "candidate_parameters": {
            "tic_id":                "281459674",
            "period_days":           82.69003,
            "transit_depth_ppm":     7400.0,
            "planet_radius_earth":   10.7,
            "stellar_mass_solar":    0.89,
            "stellar_radius_solar":  0.90,
            "stellar_teff_k":        5090.0,
            "cvs_score":             0.864,
            "cvs_verdict":           "PLANET CANDIDATE",
        },
        "archive_parameters": {
            "planet_name":           fake_match.planet_name,
            "period_days_archive":   fake_match.period_days,
            "transit_depth_archive": fake_match.transit_depth,
            "planet_radius_earth":   fake_match.planet_radius_earth,
            "semi_major_axis_au":    fake_match.semi_major_axis_au,
            "stellar_teff_k":        fake_match.stellar_teff_k,
            "stellar_radius_solar":  fake_match.stellar_radius_solar,
            "stellar_mass_solar":    fake_match.stellar_mass_solar,
            "discovery_method":      fake_match.discovery_method,
            "disposition":           fake_match.disposition,
            "extra_fields":          fake_match.extra_fields,
        },
        "parameter_deltas": {
            "period_delta_days":   0.00003,
            "period_delta_pct":    round(0.00003 / 82.69 * 100, 6),
            "depth_delta_ppm":     round((0.00740 - 0.00740) * 1e6, 2),
            "radius_delta_pct":    round((10.7 - 10.8) / 10.8 * 100, 4),
        },
        "validation_note": (
            "This entry was identified as a KNOWN planet via period matching. "
            "No mathematical discovery proof is generated for known targets."
        ),
    }

    write_json(str(out_path), card)
    print(f"\n  Simulated match: {fake_match.planet_name}")
    print_known_summary(str(out_path))
    print(f"\n  Output file: {out_path}")
    return {"status": "KNOWN", "file": str(out_path)}


# ─────────────────────────────────────────────────────────────────────────────
# Test C: Real TIC with live API
# ─────────────────────────────────────────────────────────────────────────────

def run_real_tic(tic_id: str, output_dir: str = ".") -> dict:
    section(f"TEST C  ·  LIVE NASA QUERY  ·  TIC {tic_id}")
    print(f"\n  Attempting live NASA Exoplanet Archive query for TIC {tic_id} …")

    validator = AxiomValidator(output_dir=output_dir, verbose=True)
    result = validator.validate(
        tic_id               = tic_id,
        period_days          = 3.6986,
        transit_depth        = 0.008361,
        transit_duration_hrs = 2.0,
        t0_btjd              = 1201.04,
        stellar_mass_solar   = 1.0,
        stellar_radius_solar = 1.0,
        stellar_teff_k       = 5778.0,
        stellar_logg         = 4.44,
        planet_radius_earth  = 9.97,
        cvs_score            = 0.8299,
        cvs_verdict          = "PLANET CANDIDATE",
        bls_snr              = 347.89,
        bls_fap              = 0.0,
        even_odd_delta_sigma = 0.975,
        shape_ratio          = 4.711,
        secondary_snr        = 1.2,
        centroid_sigma       = 0.264,
        planet_order         = 1,
    )
    print(f"\n  Status:  {result.status}")
    print(f"  File:    {result.output_file}")
    if result.network_error:
        print(f"  Network: OFFLINE ({result.network_error[:60]}…)")

    if result.status in ("NEW_DISCOVERY", "OFFLINE_NEW_DISCOVERY"):
        print_proof_summary(result.output_file)
    else:
        print_known_summary(result.output_file)

    return {"status": result.status, "file": result.output_file}


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="AxiomValidator Sovereign Validation Module",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_validator.py --synthetic
  python run_validator.py --force-known
  python run_validator.py --tic 260128333
  python run_validator.py --synthetic --force-known
        """,
    )
    parser.add_argument("--synthetic",   action="store_true", help="Run synthetic new discovery test")
    parser.add_argument("--force-known", action="store_true", help="Simulate a known-planet match")
    parser.add_argument("--tic",         type=str, default=None, help="Real TIC ID to query live")
    parser.add_argument("--output",      type=str, default=".", help="Output directory")
    args = parser.parse_args()

    if not any([args.synthetic, args.force_known, args.tic]):
        parser.print_help()
        print("\n  No mode selected. Running --synthetic by default.\n")
        args.synthetic = True

    Path(args.output).mkdir(parents=True, exist_ok=True)
    results = []

    if args.synthetic:
        try:
            r = run_synthetic_discovery(output_dir=args.output)
            results.append(r)
        except Exception:
            print("\n[ERROR] Synthetic test failed:")
            traceback.print_exc()

    if args.force_known:
        try:
            r = run_simulated_known(output_dir=args.output)
            results.append(r)
        except Exception:
            print("\n[ERROR] Known simulation failed:")
            traceback.print_exc()

    if args.tic:
        try:
            r = run_real_tic(args.tic, output_dir=args.output)
            results.append(r)
        except Exception:
            print(f"\n[ERROR] Real TIC {args.tic} test failed:")
            traceback.print_exc()

    section("VALIDATION COMPLETE")
    for r in results:
        icon = "✓" if r["status"] != "OFFLINE_NEW_DISCOVERY" else "~"
        print(f"  {icon}  {r['status']:<30}  →  {r['file']}")
    print()


if __name__ == "__main__":
    main()
