"""
thresholds_report.py — Generate THRESHOLDS_REPORT.md from the threshold
catalog (config/production.yaml + zspace_engine/thresholds.py).

Usage:
    python -m zspace_engine.thresholds_report [--out THRESHOLDS_REPORT.md]
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

from zspace_engine import thresholds as T

PROFILE_ORDER = ("conservative", "balanced", "sensitive")

BENCH_CMD = (
    "python benchmarks_controlled/run_controlled.py "
    "--true 100 --false 80 --out benchmarks_controlled/runs/<NAME> --seed 20260814"
)
METRICS_KNOWN = {
    "conservative": "Contamination FPR 0/80 (0.0%) on the fixed-seed suite · 4.25% (17/400) across the 800-target BIG400 suite · Recall 52/100 fixed-seed / 41.2% BIG400 — measured",
    "balanced":     "Contamination FPR 0/80 (0.0%) on the fixed-seed suite · 4.25% (17/400) across the 800-target BIG400 suite · Recall 52/100 fixed-seed / 41.2% BIG400 — measured",
    "sensitive":    "UNMEASURED — every deviation requires a full benchmark re-run",
}


def _table_rows() -> str:
    lines = ["| key | name | unit | weight | conservative | balanced | sensitive |", "|---|---|---|---|---|---|---|"]
    for k, meta in T.catalog().items():
        vals = [str(T.profile_values(p).get(k)) for p in PROFILE_ORDER]
        lines.append(
            f"| `{k}` | {meta['name']} | {meta.get('unit','')} | {meta.get('weight','')} "
            f"| {vals[0]} | {vals[1]} | {vals[2]} |"
        )
    return "\n".join(lines)


def _profile_sections() -> str:
    out = []
    for p in PROFILE_ORDER:
        out.append(
            f"### `{p}`\n\n"
            f"**Mode:** {('FPR=0 priority, measured' if p == 'conservative' else 'default, measured' if p == 'balanced' else 'recall priority, EXPERIMENTAL')}  \n"
            f"**Measured metrics:** {METRICS_KNOWN[p]}  \n\n"
            "| key | value |\n|---|---|\n"
            + "\n".join(f"| `{k}` | {T.profile_values(p).get(k)} |" for k in T.catalog())
        )
    return "\n\n".join(out)


def _per_threshold_sections() -> str:
    out = []
    for k, m in T.catalog().items():
        vals = {p: T.profile_values(p)[k] for p in PROFILE_ORDER}
        out.append(
            f"### {m['name']}  (`{k}`)\n\n"
            f"- **Unit:** {m.get('unit','—')} · **Direction:** {m.get('direction','—')} · "
            f"**Weight:** {m.get('weight','—')}\n"
            f"- **Purpose:** {m['purpose']}\n"
            f"- **Measured evidence:** {m['evidence']}\n"
            f"- **Values:** conservative `{vals['conservative']}` · balanced `{vals['balanced']}` · "
            f"sensitive `{vals['sensitive']}`\n"
            f"- **Pros of tightening:** {m.get('pros_tight','—')}\n"
            f"- **Cons of tightening:** {m.get('cons_tight','—')}\n"
            f"- **Pros of loosening:** {m.get('pros_loose','—')}\n"
            f"- **Cons of loosening:** {m.get('cons_loose','—')}\n"
            f"- **FPR risk when loosened:** {m.get('fpr_risk_loose','—')}\n"
        )
    return "\n".join(out)


def generate() -> str:
    title = "# Axiom-ZSpace Threshold Catalog — Reference Report\n"
    meta = (
        f"\n_Generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} from "
        "`config/production.yaml` + `zspace_engine/thresholds.py` — single source of truth._\n"
    )
    overview = (
        "\n## 1. What this is\n\n"
        "Every tunable number of the detection/validation pipeline lives in one place. "
        "No code edits are needed to change a threshold: edit "
        "`config/production.yaml` → section `thresholds`, then re-run the benchmark.\n\n"
        "## 2. Three pre-defined profiles\n\n"
        "Profiles trade **false-positive rate (FPR)** against **recall**. "
        "Changing the active profile: `python -m zspace_engine.thresholds --set balanced` "
        "(or edit `thresholds.profile` in the YAML).\n\n"
        "**Always re-measure after any change** — the values below marked [M] were "
        f"measured on the controlled benchmark (command below):\n\n    {BENCH_CMD}\n"
    )
    howto = (
        "\n## 3. Decision rules used by the validator\n\n"
        "- All **critical** gates must pass (`critical_passed`), otherwise the candidate is "
        "FALSE_POSITIVE (circuit breaker).\n"
        "- With all critical gates passing: ≤`verdict_max_fail_pass` total fails → "
        "SOVEREIGN_PASS; ≤`verdict_max_fail_conditional` → CONDITIONAL_PASS; more → FALSE_POSITIVE.\n"
        "- `coherent_override_enabled` lets repeated-observation evidence override FP-2 "
        "(the power-FAP firewall). It is OFF in measured profiles because probe runs "
        "measured contamination FPR up to 62.5% (5/8) when ON (PROBE_FPR68 series).\n"
    )
    table = "\n## 4. All thresholds at a glance\n\n" + _table_rows() + "\n"
    profile_detail = "\n## 5. Profile values in full\n\n" + _profile_sections() + "\n"
    detail = "\n## 6. Per-threshold rationale (measured evidence, pros/cons)\n\n" + _per_threshold_sections() + "\n"
    workflow = (
        "\n## 7. Change workflow\n\n"
        "1. Edit `config/production.yaml` (one value, or build a new profile block).\n"
        "2. Verify thresholds load: `python -m zspace_engine.thresholds --show`.\n"
        "3. Re-measure: " + BENCH_CMD + "\n"
        "4. Compare FPR/recall against the measured baselines in section 5; "
        "update this report with the new numbers.\n"
    )
    return title + meta + overview + howto + table + profile_detail + detail + workflow


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="THRESHOLDS_REPORT.md")
    args = parser.parse_args()
    path = Path(args.out)
    path.write_text(generate(), encoding="utf-8")
    print(f"Wrote {path.resolve()}")


if __name__ == "__main__":
    main()