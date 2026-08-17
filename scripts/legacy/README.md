# scripts/legacy — pre-v1 tools (not part of the v1.0 pipeline)

This directory holds tools from the internal v2.x lineage that are **not**
part of the calibrated v1.0 Python pipeline and are **not** referenced by any
v1.0 test, benchmark, or documentation. They are shipped for history and for
internal side-tasks only.

Contents
--------

| File | What it was for |
|---|---|
| `axiom_zspace.c` | A standalone C port of the detection engine (v2-era). Never compiled in this repo; its CVS classification threshold is hardcoded (0.80/0.55/0.35). |
| `scan_sector.sh` | Driver that exported sector light curves to CSV and ran the C binary. Expects to be run from this directory (`./axiom_zspace`). |
| `export_sector_csv.py` | Phase-1 CSV exporter used by `scan_sector.sh`. |
| `benchmark_known_planets.py` | Pre-v1 benchmark harness against known planets (superseded by `benchmarks_*`). |
| `benchmark_query_probe.py` | Scratch probe for querying the NASA Exoplanet Archive. |

Notes
-----

- The v1.0 measured evidence lives exclusively in `benchmarks_controlled/`
  and `benchmarks_real/` — nothing in this directory contributes claims.
- v1.0 did not preserve the batch/sector outputs these tools produced;
  see CHANGELOG.txt "LEGACY ENTRIES".
- If you build the C engine, do it here:
  `gcc -O3 -fopenmp -o axiom_zspace axiom_zspace.c -lm`
  (the compiled binary `axiom_zspace` is git-ignored).