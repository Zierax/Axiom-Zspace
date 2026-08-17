#!/usr/bin/env python3
"""Aggregate chunked benchmark runs into one BIG suite evidence dir.

Large controlled suites are run as independent chunk runs (deterministic,
seed-per-chunk) so they can execute in parallel; this script merges the chunk
results and recomputes every metric with the SAME definitions as the canonical
runner (benchmarks_controlled.run_controlled.compute_metrics), then writes a
combined EVALUATION_REPORT.md plus merged results JSONs.

Usage:
    python scripts/aggregate_benchmark_runs.py \
        --chunks-dir benchmarks_controlled/runs/BIG400_chunks \
        --out benchmarks_controlled/evidence/BIG400 \
        --suite-name BIG400 \
        --suite-seed 20260816 \
        --per-chunk 50 50
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmarks_controlled.run_controlled import compute_metrics

CERTIFIED = ("NEW_DISCOVERY", "OFFLINE_NEW_DISCOVERY")


def _chunk_paths(chunks_dir: Path) -> list[Path]:
    return sorted(p for p in chunks_dir.iterdir() if p.is_dir())


def _load_results(chunk: Path, name: str) -> list[dict]:
    with open(chunk / f"results_{name}.json", "r", encoding="utf-8") as fh:
        return json.load(fh)


def _snr_tables(true_res, false_res) -> tuple[list[str], list[str]]:
    buckets = defaultdict(lambda: [0, 0])
    for r in true_res:
        ts = r.get("target_snr")
        perr = r.get("period_error_pct")
        if ts is None:
            continue
        b = "5.5-8" if ts < 8 else ("8-14" if ts < 14 else "14-30")
        buckets[b][1] += 1
        if r["validation_status"] in CERTIFIED and perr is not None and perr <= 5.0:
            buckets[b][0] += 1
    snr_lines = ["| SNR tier | recovered/n | recall |", "|---|---|---|"]
    for b in ("5.5-8", "8-14", "14-30"):
        rec, tot = buckets.get(b, (0, 0))
        pct = rec / tot * 100 if tot else 0.0
        snr_lines.append(f"| {b} | {rec}/{tot} | {pct:.1f}% |")

    fb = defaultdict(list)
    for r in false_res:
        fb[r["subkind"]].append(r)
    kind_lines = ["| kind | count | certified-FP |", "|---|---|---|"]
    for kind, arr in sorted(fb.items()):
        cert = sum(1 for r in arr if r["validation_status"] in CERTIFIED)
        kind_lines.append(f"| {kind} | {len(arr)} | {cert} |")
    return snr_lines, kind_lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chunks-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--suite-name", default="BIG400")
    ap.add_argument("--suite-seed", type=int, default=20260816)
    ap.add_argument("--per-chunk", nargs=2, type=int, default=[50, 50])
    args = ap.parse_args()

    chunks = _chunk_paths(args.chunks_dir)
    if not chunks:
        sys.exit(f"no chunk dirs found under {args.chunks_dir}")

    true_all, false_all, provenance = [], [], []
    for c in chunks:
        tru, fal = _load_results(c, "true"), _load_results(c, "false")
        for row in tru:
            row.setdefault("chunk", c.name)
        for row in fal:
            row.setdefault("chunk", c.name)
        true_all += tru
        false_all += fal
        provenance.append({"chunk": c.name, "true_targets": len(tru),
                           "false_targets": len(fal)})
        if len(tru) != args.per_chunk[0] or len(fal) != args.per_chunk[1]:
            sys.exit(f"{c.name}: expected {args.per_chunk} results, got {len(tru)}/{len(fal)}")

    m = compute_metrics(true_all, false_all)

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "results_true.json").write_text(json.dumps(true_all, indent=2), encoding="utf-8")
    (args.out / "results_false.json").write_text(json.dumps(false_all, indent=2), encoding="utf-8")
    (args.out / "chunks.json").write_text(json.dumps(
        {"suite": args.suite_name, "suite_seed": args.suite_seed,
         "per_chunk": args.per_chunk, "chunks": provenance}, indent=2), encoding="utf-8")

    snr_lines, kind_lines = _snr_tables(true_all, false_all)
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    report = [
        "# Axiom-ZSpace Aggregated Controlled Evaluation",
        "",
        f"**Generated:** {now}",
        f"**Suite:** {args.suite_name} (aggregated from {len(chunks)} chunk runs)",
        f"**Composition:** {len(true_all)} true + {len(false_all)} false targets; "
        f"chunks of {args.per_chunk[0]}+{args.per_chunk[1]} with seeds "
        f"{args.suite_seed}..{args.suite_seed + len(chunks) - 1} (seed = suite_seed + chunk index).",
        "**Method:** full local pipeline, BLIND search (no period prior, no ground-truth",
        "expected period), validator archive query stubbed OFFLINE so the sovereign",
        "verdict reflects ONLY intrinsic physics/QA gating. Metrics computed with the",
        "canonical compute_metrics() of benchmarks_controlled/run_controlled.py.",
        "",
        "## Key Metrics",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| True (injected) planets | {m['n_true']} |",
        f"| Contamination targets | {m['n_false']} |",
        f"| **Recall@correct period (TP)** | {m['recall_period']*100:.1f}%  ({m['tp_period']}/{m['n_true']}) |",
        f"| **Detection recall (any period)** | {m['recall_any']*100:.1f}%  ({m['tp_any']}/{m['n_true']}) |",
        f"| Detected-at-any-period (raw BLS) | {m['detected_any_true']}/{m['n_true']} |",
        f"| **FPR (contamination certified as planet)** | {m['fpr']*100:.2f}%  ({m['fp_contamination']} FP / {m['n_false']} ) |",
        f"| Wrong-ephemeris certs on TRUE set | {m['fp_wrong_ephem']} (certified at a period outside the 5% window) |",
        f"| **Precision (certified & correct-period)** | {m['precision']*100:.1f}% |",
        f"| **F1 (period-level)** | {m['f1']:.3f} |",
        f"| Confusion | TP={m['tp_period']} FP={m['fp']} FN={m['fn']} TN={m['tn']} |",
        "",
        "## Recall by injected SNR",
        "",
        *snr_lines,
        "",
        "## Per-subkind contamination results",
        "",
        *kind_lines,
        "",
        "## Provenance (chunks)",
        "",
        "| chunk | true | false | seed |",
        "|---|---|---|---|",
    ]
    for i, p in enumerate(provenance):
        report.append(f"| {p['chunk']} | {p['true_targets']} | {p['false_targets']} | {args.suite_seed + i} |")
    report += [
        "",
        "Each chunk is independently reproducible from its seed with run_controlled.py.",
        "",
        "## Reproduce this suite",
        "",
        "```bash",
        f"python benchmarks_controlled/run_controlled.py --true {args.per_chunk[0]} --false {args.per_chunk[1]} \\",
        f"    --out benchmarks_controlled/runs/<cN> --seed <suite_seed + N>    # every chunk N in 0..{len(chunks)-1}",
        f"python scripts/aggregate_benchmark_runs.py --chunks-dir benchmarks_controlled/runs/<dir> \\",
        f"    --out {args.out} --suite-name {args.suite_name} --suite-seed {args.suite_seed}",
        "```",
        "",
        "## Files",
        "",
        "- `results_true.json` — merged per-target results for injected planets",
        "- `results_false.json` — merged per-target results for contamination",
        "- `chunks.json` — chunk/seed provenance",
        "- `EVALUATION_REPORT.md` — this report",
        "",
        "**Note on honesty:** wrong-ephemeris certifications on the true set are counted",
        "separately from contamination false-certifications (the FPR row is contamination-only).",
        "",
    ]
    (args.out / "EVALUATION_REPORT.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Wrote {args.out / 'EVALUATION_REPORT.md'}")
    print(f"TP={m['tp_period']} FP={m['fp_contamination']} WRONG_EPHEM_TRUE={m['fp_wrong_ephem']} "
          f"recall={m['recall_period']*100:.1f}% fpr={m['fpr']*100:.2f}%")


if __name__ == "__main__":
    main()