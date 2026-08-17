# Contributing

Axiom-ZSpace v1.0 is a *measured* pipeline: its credibility is the evidence
in `benchmarks_*/evidence/`, the honesty of its docs, and the 101-test
regression suite. Every contribution is judged against those three.

## The one rule

> **Nothing in the threshold catalog changes without re-measuring on BOTH
> benchmarks.**

`config/production.yaml` (section `thresholds`) and
`zspace_engine/thresholds.py` are data, not code. A threshold edit implies:

1. edit the catalog (one file, one value, one comment saying *why*);
2. re-run the controlled benchmark (`docs/BENCHMARKS.md` §5);
3. re-run the real benchmark;
4. regenerate `THRESHOLDS_REPORT.md`;
5. update README/docs numbers to the **new measured values** and cite the
   run name; commit evidence JSONs together with the config change.

A PR that ships a threshold change without the two measurements will be
rejected. That is not bureaucracy — it is what makes the repo's numbers
meaningful.

## Claim discipline (documentation)

- Every number in docs/README/CHANGELOG must trace to an **artifact in this
  repo**: `benchmarks_*/evidence/*/results_*.json`, a test assertion, or the
  generated `THRESHOLDS_REPORT.md`.
- No "historical" numbers from `archive/` or from memory. The CHANGELOG v2.x
  "LEGACY ENTRIES" section exists precisely because unverifiable numbers were
  once claimed — do not reintroduce that class of claim.
- If you cannot attach an artifact, you may not ship the number.

## Hygiene rules (hard)

- `archive/` is **never** committed, pushed, or referenced as evidence.
- Secrets (API keys, tokens, credential-bearing scripts such as
  `get_models*.ps1`) never enter the repo; `.gitignore` enforces the known
  patterns — extend it when you add another.
- Run outputs (`runs/`, caches, `axiom_output/`, `Discovery_*.json`) stay
  git-ignored. Only evidence runs get versioned (under `benchmarks_*/evidence/`).
- No scratch files at the repo root. Utilities live in `scripts/`; pre-v1
  tools live in `scripts/legacy/` (with its README updated).

## Workflow

1. Open an issue first for behavioral changes (gates, thresholds, verdict
   semantics) — the discussion must state the expected effect on FPR/recall.
2. Branch from `main`; commits one-logical-change each, message style:
   `fix:`, `feat:`, `docs:`, `test:`, `bench:` prefixes
   (e.g. `fix: unify cvs_planet_threshold at 0.80`).
3. Run the gate before push:

```bash
python -m pytest tests/ -q                    # 101 passed
python -m zspace_engine.thresholds_report     # if thresholds touched
```

4. Update the docs that the change touches:
   - threshold catalog / semantics → `THRESHOLDS_REPORT.md` (generated) +
     `docs/ARCHITECTURE.md` (§6–§8) if emission rules changed;
   - benchmark numbers → `docs/BENCHMARKS.md` + root README table + CHANGELOG;
   - new test file → `docs/TESTING.md` (suite map);
   - new CLI surface → `docs/QUICKSTART.md` + `run_pipeline.py --help`;
   - new module → `docs/ARCHITECTURE.md` module register.

## Code style

- Python ≥ 3.10, type hints on public signatures, dataclasses for result
  objects (mirror the existing auditors/context result types).
- No comments that restate code; comments explain *physics* and *evidence*.
- Results flow through dataclasses with an evidence string — never return a
  bare float where a reviewer needs to know *why*.
- Reuse: detectors/auditors take thresholds from the catalog
  (`thresholds.threshold(key, profile=…)`), never hardcode a tunable.

## Review checklist (maintainers)

```
□ 101 tests pass (offline)
□ benchmark artifacts committed with any measurement
□ numbers in the diff trace to artifacts
□ no archive/ secrets run outputs in the diff
□ config/docs/THRESHOLDS_REPORT regenerated & consistent
□ CHANGELOG updated for user-visible change
```