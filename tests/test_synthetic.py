#!/usr/bin/env python3
"""
test_synthetic.py  ·  Synthetic Planet Test Generation
=======================================================
Provides synthetic test data generation for validating the Axiom-ZSpace
discovery pipeline without requiring network access or real TESS data.

This module contains reusable synthetic test generation functions that can
be used by various test scripts.
"""

import json
from typing import Dict, Any


def generate_synthetic_planet_params() -> Dict[str, Any]:
    """
    Generate synthetic planet parameters for testing.
    
    Returns a dictionary of parameters that simulate a realistic exoplanet
    candidate that will not appear in the NASA Archive (for testing NEW_DISCOVERY path).
    
    Returns:
        dict: Synthetic planet parameters including orbital, physical, and detection metrics
    """
    return {
        "tic_id":                "SYNTHETIC-001",
        "period_days":           3.69864,
        "transit_depth":         0.008361,
        "transit_duration_hrs":  2.0,
        "t0_btjd":               1201.04,
        "stellar_mass_solar":    1.0,
        "stellar_radius_solar":  1.0,
        "stellar_teff_k":        5778.0,
        "stellar_logg":          4.44,
        "planet_radius_earth":   9.97,
        "cvs_score":             0.8299,
        "cvs_verdict":           "PLANET CANDIDATE",
        "cvs_proof_chain": [
            "  [periodicity] w=0.97 · S=1.0000 = 0.9700 | SNR=347.89 > 5.5 → PASS",
            "  [depth]       w=0.83 · S=0.7955 = 0.6603 | CV=0.0205 → S_δ=0.7955",
            "  [limb]        w=0.61 · S=1.0000 = 0.6100 | ratio=4.711 → U-shape",
            "  [stellar]     w=0.31 · S=0.0554 = 0.0172 | secondary flagged → penalty",
            "  CVS = 2.2575 / 2.7200 = 0.8299",
        ],
        "bls_snr":               347.89,
        "bls_fap":               0.0,
        "even_odd_delta_sigma":  0.975,
        "shape_ratio":           4.711,
        "secondary_snr":         1.2,
        "centroid_sigma":        0.264,
        "limb_dark_u1":          0.30,
        "limb_dark_u2":          0.10,
        "zspace_id":             "ZS-T-SYNTHETIC-001-01",
        "planet_order":          1,
    }


def generate_known_planet_params() -> Dict[str, Any]:
    """
    Generate parameters for a known planet (simulates TOI-4600 b).
    
    Returns:
        dict: Parameters matching a known exoplanet for testing KNOWN path
    """
    return {
        "tic_id":                "281459674",
        "period_days":           82.69003,
        "transit_depth":         0.00740,
        "transit_duration_hrs":  3.5,
        "t0_btjd":               1500.0,
        "stellar_mass_solar":    0.89,
        "stellar_radius_solar":  0.90,
        "stellar_teff_k":        5090.0,
        "stellar_logg":          4.50,
        "planet_radius_earth":   10.7,
        "cvs_score":             0.864,
        "cvs_verdict":           "PLANET CANDIDATE",
        "bls_snr":               250.0,
        "bls_fap":               0.0,
        "even_odd_delta_sigma":  0.95,
        "shape_ratio":           4.5,
        "secondary_snr":         0.8,
        "centroid_sigma":        0.2,
        "planet_order":          1,
    }


def print_discovery_summary(discovery_path: str) -> None:
    """
    Print a formatted summary of a Discovery Card JSON file.
    
    Args:
        discovery_path: Path to the Discovery JSON file
    """
    with open(discovery_path) as f:
        card = json.load(f)

    sections = card.get("proof_sections", {})
    print(f"\n  ┌{'─'*66}┐")
    print(f"  │  SOVEREIGN LOGIC CARD  ·  {card.get('zspace_id',''):<37}│")
    print(f"  │  Status: {card.get('status',''):<56}│")
    print(f"  │  Verdict: {card.get('sovereign_verdict',''):<55}│")
    print(f"  └{'─'*66}┘")

    for sec_key, sec in sections.items():
        label = sec.get("section", sec_key)
        verdict = sec.get("verdict", sec.get("overall_verdict", "?"))
        result = sec.get("result", "")
        unit = sec.get("unit", "")
        proof = sec.get("proof", "")
        # Truncate proof for display
        proof_short = proof[:120] + "…" if len(proof) > 120 else proof
        icon = "✓" if "PASS" in verdict or "SOVEREIGN" in verdict else ("✗" if "FAIL" in verdict else "⚠")
        print(f"\n  {icon} {label}  [{verdict}]")
        if result:
            print(f"      Result: {result} {unit}")
        print(f"      Proof:  {proof_short}")

    phys = card.get("physical_summary", {})
    print(f"\n  ─── Physical Summary ───────────────────────────────────────────")
    for k, v in phys.items():
        print(f"      {k:<40} {v}")

    fp = sections.get("section_5_false_positive_ruling", {})
    if fp:
        print(f"\n  ─── False-Positive Test Results ────────────────────────────────")
        for t in fp.get("tests", []):
            icon2 = "✓" if t["verdict"] == "PASS" else "✗"
            print(f"      {icon2}  {t['test']:<32} {t['comparison']}")
        print()
        for line in fp.get("logical_closure", []):
            print(f"      └▶ {line}")


def print_known_summary(known_path: str) -> None:
    """
    Print a formatted summary of an exist_planet JSON file.
    
    Args:
        known_path: Path to the exist_planet JSON file
    """
    with open(known_path) as f:
        card = json.load(f)
    print(f"\n  ┌{'─'*66}┐")
    print(f"  │  EXISTING PLANET CARD  ·  {card.get('zspace_id',''):<37}│")
    print(f"  │  Match: {card['match_summary']['planet_name']:<57}│")
    print(f"  └{'─'*66}┘")
    ms = card["match_summary"]
    print(f"\n      Archive match: {ms['planet_name']}")
    print(f"      Period delta:  {ms['period_delta_days']:.8f} d  (tolerance ±{ms['tolerance_days']} d)")
    print(f"      Proof:         {ms['proof']}")
    ap = card.get("archive_parameters", {})
    print(f"\n      Archive parameters:")
    for k, v in ap.items():
        if v is not None and k not in ("extra_fields",):
            print(f"        {k:<35} {v}")
    deltas = card.get("parameter_deltas", {})
    print(f"\n      Candidate vs Archive deltas:")
    for k, v in deltas.items():
        if v is not None:
            print(f"        {k:<35} {v}")
