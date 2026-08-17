#!/usr/bin/env python3
"""
fetch_nea_snapshot.py — regenerate the offline NEA snapshots for the real
benchmark sample.

The real-data benchmark (``benchmarks_real/run_real.py``) must be reproducible
without network access. These files pin the NASA Exoplanet Archive state the
deterministic sample is drawn from:

  real_ps_star.json       Kepler "planet search" host stars with confirmed
                          planets in 1.0 <= P <= 13.5 d + stellar params
                          (the true-target pool).
  real_quiet.json         Kepler DR25 stars with ``nkoi = 0 AND nconfp = 0``
                          and full stellar params (the quiet-star proxy pool).
  real_known_signals.json Cumulative KOI rows per selected kepid, used by the
                          runner to classify found periods against known
                          planets / cataloged FPs.

Run BEFORE re-measuring the real benchmark if the archive may have changed:

    python scripts/fetch_nea_snapshot.py --out benchmarks_real/data
"""

import argparse
import csv
import io
import json
from datetime import datetime, timezone
from pathlib import Path

import requests

TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
NOW = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def tap(sql: str) -> list:
    r = requests.get(TAP, params={"query": sql, "format": "csv"}, timeout=600)
    r.raise_for_status()
    return list(csv.DictReader(io.StringIO(r.text)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "benchmarks_real" / "data"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print("fetching true-target pool (Kepler PS hosts, 1.0<=P<=13.5 d) ...")
    true_rows = tap(
        "SELECT k.kepid, k.kepmag, k.teff, k.logg, k.radius, k.mass, "
        "p.pl_name, p.hostname, p.pl_orbper, p.pl_rade "
        "FROM ps k JOIN pscomppars p ON p.hostname = k.hostname "
        "WHERE p.default_flag = 1 "
        "AND p.pl_orbper >= 1.0 AND p.pl_orbper <= 13.5 "
        "AND k.teff IS NOT NULL AND k.logg IS NOT NULL "
        "AND k.radius IS NOT NULL AND k.mass IS NOT NULL "
    )
    # one row per KIC (dedupe multi-planet hosts the same way run_real does)
    seen, deduped = set(), []
    for r in true_rows:
        if r["kepid"] not in seen:
            seen.add(r["kepid"])
            deduped.append(r)
    print(f"  {len(deduped)} host stars")

    print("fetching quiet-star pool (DR25 nkoi=0 AND nconfp=0) ...")
    quiet_rows = tap(
        "SELECT kepid, kepmag, teff, logg, radius, mass, nkoi "
        "FROM keplerstellar17 "
        "WHERE nkoi = 0 AND nconfp = 0 "
        "AND teff IS NOT NULL AND logg IS NOT NULL "
        "AND radius IS NOT NULL AND mass IS NOT NULL "
    )
    print(f"  {len(quiet_rows)} stars")

    # Plain-row lists: the runner expects json.loads() → iterable of rows, no
    # wrapper dict. Fetch provenance goes in the sidecar below.
    (out / "real_ps_star.json").write_text(json.dumps(deduped, indent=1), encoding="utf-8")
    (out / "real_quiet.json").write_text(json.dumps(quiet_rows, indent=1), encoding="utf-8")

    sidecar = {
        "fetched_utc": NOW,
        "source": "NASA Exoplanet Archive TAP (ps, pscomppars, keplerstellar17)",
        "query": "see scripts/fetch_nea_snapshot.py",
        "note": ("real_known_signals.json is the runner's own per-kepid KOI "
                 "cache and is NOT regenerated here — keep the committed "
                 "snapshot and let run_real.py fill it on first run."),
    }
    (out / "_snapshot_meta.json").write_text(json.dumps(sidecar, indent=1), encoding="utf-8")

    print("done. re-run the real benchmark to refresh runs against this snapshot.")


if __name__ == "__main__":
    main()