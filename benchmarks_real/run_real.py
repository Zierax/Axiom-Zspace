#!/usr/bin/env python3
"""
run_real.py  ·  Real-Data Honesty Benchmark (Kepler, NEA-confirmed truth)
==========================================================================
Runs the SAME blind pipeline as the controlled benchmark (ingestion → BLS →
auditors → sovereign validator, archive stub OFFLINE) on REAL Kepler long-
cadence light curves:

  * true  subset: hosts of NEA-CONFIRMED transit planets (period 1.0-13.5 d,
                  Kepler mission), truth = NEA period
  * false subset: Kepler stars with ZERO KOIs / zero confirmed planets in the
                  DR25 stellar delivery (proxy for quiet stars; the real FPR
                  is "rate of certifying stars with no KNOWN planet")

Scoring is identical to the controlled run (Recall@period ≤5%, FPR, precision).

Honest caveats (written into the report):
  * quiet-star FPR is a PROXY: unknown planets may exist below detection limits
  * 5% period window is an acceptance rule, not a truth claim at that precision
  * long cadence (29.4 min) + quarter gaps + real systematics, unlike synthetic

Usage:
  python benchmarks_real/run_real.py [--n-true 12] [--n-false 12] [--out DIR]
"""
import argparse
import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from zspace_engine.logging_config import setup_logging, get_logger
import benchmarks_controlled.run_controlled as RC

setup_logging()
logger = get_logger(__name__)

CACHE_DIR = Path(__file__).resolve().parent / "cache"
DATA_DIR = Path(__file__).resolve().parent / "data"
QUARTERS = list(range(4, 10))  # Q4..Q9 ≈ 540 d baseline, plenty of transits


class RealTarget:
    """Adapter making a real Kepler light curve look like a SyntheticTarget:
    evaluate_target() is reused verbatim (blind search)."""

    def __init__(self, tic_id: str, kind: str, subkind: str, label_period: Optional[float],
                 time: np.ndarray, flux: np.ndarray, quality: np.ndarray,
                 stellar: dict, white_ppm: float, meta_extra: dict):
        self.tic_id = tic_id
        self.kind = kind
        self.subkind = subkind
        self.label_period = label_period
        self.time = np.asarray(time)
        self.flux = np.asarray(flux)
        self.quality = np.asarray(quality)
        self.stellar = stellar
        self.injected_depth = None
        self.meta = {"white_ppm": white_ppm, "target_snr": None, **meta_extra}


def load_or_fetch(kepid: int) -> Optional[Dict]:
    """Download Kepler Q4-Q9 long-cadence, return {time, flux, quality}."""
    import lightkurve as lk

    cache_file = CACHE_DIR / f"KIC{kepid}.npz"
    if cache_file.exists():
        d = np.load(cache_file)
        return {"time": d["time"], "flux": d["flux"], "quality": d["quality"]}

    try:
        s = lk.search_lightcurve(target=int(kepid), mission="Kepler", cadence="long", quarter=QUARTERS)
        if len(s) == 0:
            logger.warning(f"KIC{kepid}: no Q4-Q9 long-cadence data, skipping")
            return None
        lcs = s.download_all(download_dir=str(CACHE_DIR / "dl"), quality_bitmask="default")
        if lcs is None or len(lcs) == 0:
            return None

        times, fluxes, quals = [], [], []
        for lc in lcs:
            f = np.asarray(lc.flux.value)
            t = np.asarray(lc.time.value)
            q = np.asarray(lc.quality.value, dtype=int)
            good = np.isfinite(f) & np.isfinite(t) & (q == 0) & (f > 0)
            t, f, q = t[good], f[good], q[good]
            if len(t) < 100:
                continue
            m = np.median(f)
            times.append(t)
            fluxes.append(f / m)
            quals.append(q)

        if not times:
            return None
        t = np.concatenate(times)
        f = np.concatenate(fluxes)
        q = np.concatenate(quals)
        t = np.asarray(t)
        f = np.asarray(f)
        q = np.asarray(q)
        order = np.argsort(t)
        t, f, q = t[order], f[order], q[order]

        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache_file, time=t, flux=f, quality=q)
        logger.info(f"KIC{kepid}: cached {len(t)} pts, baseline {t[-1]-t[0]:.1f} d")
        return {"time": t, "flux": f, "quality": q}
    except Exception as e:
        logger.error(f"KIC{kepid}: fetch failed: {e}")
        return None


def to_stellar(row: dict) -> dict:
    """KEPLERSTELLAR/ps row -> pipeline stellar dict."""
    return {
        "st_mass": float(row["mass"] or 0.6),
        "st_rad": float(row["radius"] or 0.6),
        "st_teff": float(row["teff"] or 3900.0),
        "st_logg": float(row["logg"] or 4.66),
    }


def white_ppm(flux: np.ndarray) -> float:
    return float(np.nanstd(flux) * 1e6) if np.isfinite(flux).any() else 0.0


def known_signals(kepids) -> dict:
    """NEA cumulative rows for the sample: known planets (non-FP) and cataloged
    false positives (EB etc.) per kepid. Cached snapshot, offline afterwards."""
    cache = DATA_DIR / "real_known_signals.json"
    if cache.exists():
        d = json.loads(cache.read_text())
    else:
        d = {}
    missing = [k for k in kepids if str(k) not in d]
    if missing:
        import requests, csv, io
        TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
        for i in range(0, len(missing), 12):
            chunk = missing[i:i + 12]
            sql = ("SELECT kepid, koi_period, koi_pdisposition, kepler_name, koi_kepmag "
                   "FROM cumulative WHERE kepid IN (%s)" % ",".join(map(str, chunk)))
            r = requests.get(TAP, params={"query": sql, "format": "csv"}, timeout=240)
            if r.status_code == 200:
                for row in csv.DictReader(io.StringIO(r.text)):
                    d.setdefault(str(row["kepid"]), []).append(row)
            else:
                logger.warning(f"known_signals chunk failed: {r.status_code}")
        cache.write_text(json.dumps(d), encoding="utf-8")
    return {int(k): v for k, v in d.items()}


def match_known(period_found: Optional[float], rows) -> str:
    """Classify a certified period against the NEA cumulative catalog of the
    SAME host: 'target' handled by caller; here returns a tag for any match:
      'planet' (any non-FP KOI period within 5%), 'fp' (cataloged FP within 5%),
      'none' (no known signal at this period)."""
    if period_found is None:
        return "none"
    for row in rows:
        try:
            p = float(row["koi_period"])
        except (TypeError, ValueError):
            continue
        if abs(period_found - p) / p <= 0.05:
            return "fp" if row.get("koi_pdisposition") == "FALSE POSITIVE" else "planet"
    return "none"


def real_metrics(true_results, false_results) -> dict:
    """Honest real-data metrics. For true targets a certification is:
      tp_target    — matches the NEA period of the target planet (<=5%)
      tp_other     — matches ANOTHER known (non-FP) planet of the same host
      fp_unknown   — certified but no known signal at the found period
    For quiet stars any certification is a fp (no known signals by sample
    construction)."""
    CERT = RC.CERTIFIED
    tp_target = tp_other = fp_unknown = fp_quiet = tn_quiet = 0
    for r in true_results:
        cert = r["validation_status"] in CERT
        if not cert:
            continue
        tag = match_known(r.get("period_found"), r.get("known_rows", []))
        if tag == "planet":
            if r.get("matched_target"):
                tp_target += 1
            else:
                tp_other += 1
        else:
            fp_unknown += 1
    for r in false_results:
        if r["validation_status"] in CERT:
            fp_quiet += 1
        else:
            tn_quiet += 1
    n_true = len(true_results)
    n_false = len(false_results)
    recall_target = tp_target / n_true if n_true else 0.0
    recall_known = (tp_target + tp_other) / n_true if n_true else 0.0
    fpr_quiet = fp_quiet / n_false if n_false else 0.0
    fp = fp_unknown + fp_quiet
    tn = tn_quiet
    fpr_total = fp / (fp + tn) if (fp + tn) else 0.0
    precision = tp_target / (tp_target + fp) if (tp_target + fp) else 0.0
    return {
        "n_true": n_true, "n_false": n_false,
        "tp_target": tp_target, "tp_other": tp_other,
        "fp_quiet": fp_quiet, "fp_unknown": fp_unknown,
        "fp": fp, "tn": tn,
        "recall_target": round(recall_target, 4),
        "recall_known": round(recall_known, 4),
        "fpr_quiet": round(fpr_quiet, 4),
        "fpr_total": round(fpr_total, 4),
        "precision": round(precision, 4),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-true", type=int, default=12)
    ap.add_argument("--n-false", type=int, default=12)
    ap.add_argument("--out", default=None)
    ap.add_argument("--seed", type=int, default=20260816)
    args = ap.parse_args()

    run_dir = Path(args.out) if args.out else (
        Path(__file__).resolve().parent / "runs" /
        datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S"))
    (run_dir / "validation").mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    # ── sample selection (deterministic from the NEA snapshot) ──────────────
    confirmed = json.loads((DATA_DIR / "real_ps_star.json").read_text())
    quiet = json.loads((DATA_DIR / "real_quiet.json").read_text())

    # drop entries without all stellar params
    confirmed = [c for c in confirmed if c.get("radius") and c.get("mass")
                 and c.get("teff") and c.get("logg")]
    quiet = [q for q in quiet if q.get("radius") and q.get("mass")
             and q.get("teff") and q.get("logg")]

    # spread uniformly across planet radius (depth proxy) — mimics a balanced
    # recall test rather than an easy-depth-dominated one
    confirmed.sort(key=lambda c: float(c["pl_rade"] or 0))
    step = max(1, len(confirmed) // max(1, args.n_true)) if args.n_true else len(confirmed) + 1
    true_picks = [confirmed[i] for i in range(0, len(confirmed), step)][:args.n_true]

    quiet.sort(key=lambda q: int(q["kepid"]))
    step2 = max(1, len(quiet) // max(1, args.n_false)) if args.n_false else len(quiet) + 1
    false_picks = [quiet[i] for i in range(0, len(quiet), step2)][:args.n_false]

    logger.info(f"true picks: {len(true_picks)}, false picks: {len(false_picks)}")
    for p in true_picks:
        logger.info("  T  KIC%s P=%.4f d rade=%s" % (p["kepid"], float(p["pl_orbper"]), p["pl_rade"]))
    for p in false_picks:
        logger.info("  F  KIC%s" % p["kepid"])

    all_kepids = [int(p["kepid"]) for p in true_picks + false_picks]
    known = known_signals(all_kepids)
    for p in true_picks + false_picks:
        krows = known.get(int(p["kepid"]), [])
        names = {r.get("kepler_name") for r in krows if r.get("kepler_name")}
        if names:
            p["kepler_name"] = " / ".join(sorted(names))

    results_true, results_false = [], []
    for i, p in enumerate(true_picks, 1):
        kepid = int(p["kepid"])
        data = load_or_fetch(kepid)
        if data is None:
            logger.warning(f"KIC{kepid}: SKIPPED (no data)")
            results_true.append({"tic_id": str(kepid), "error": "no data",
                                 "kind": "true_planet", "subkind": "real_kepler",
                                 "label_period": float(p["pl_orbper"]),
                                 "detected": False, "validation_status": "NO_DATA"})
            continue
        target = RealTarget(
            tic_id=f"KIC{kepid}", kind="true_planet", subkind="real_kepler_confirmed",
            label_period=float(p["pl_orbper"]),
            time=data["time"], flux=data["flux"], quality=data["quality"],
            stellar=to_stellar(p), white_ppm=white_ppm(data["flux"]),
            meta_extra={"grp": "confirmed"},
        )
        r = RC.evaluate_target(target, run_dir)
        if r.get("error"):
            logger.warning(f"KIC{kepid}: first eval error '{r['error']}', retrying once")
            r = RC.evaluate_target(target, run_dir)
        r["kepid"] = kepid
        r["nea_period"] = float(p["pl_orbper"])
        r["pl_rade"] = p.get("pl_rade")
        rows = known.get(kepid, [])
        r["known_rows"] = [{"koi_period": x["koi_period"],
                            "koi_pdisposition": x.get("koi_pdisposition"),
                            "kepler_name": x.get("kepler_name")} for x in rows]
        pf = r.get("period_found")
        r["matched_target"] = bool(pf and abs(pf - r["nea_period"]) / r["nea_period"] <= 0.05)
        r["match_known"] = match_known(pf, rows)
        results_true.append(r)
        json.dump(results_true, open(run_dir / "results_true.json", "w"), indent=2)
        logger.info(f"[true {i}/{len(true_picks)}] KIC{kepid} -> {r['validation_status']} "
                    f"snr={r.get('snr')} perr={r.get('period_error_pct')}")

    for i, p in enumerate(false_picks, 1):
        kepid = int(p["kepid"])
        data = load_or_fetch(kepid)
        if data is None:
            logger.warning(f"KIC{kepid}: SKIPPED (no data)")
            results_false.append({"tic_id": str(kepid), "error": "no data",
                                  "kind": "false_quiet", "subkind": "real_quiet",
                                  "label_period": None, "detected": False,
                                  "validation_status": "NO_DATA"})
            continue
        target = RealTarget(
            tic_id=f"KIC{kepid}", kind="false_quiet", subkind="real_quiet_star",
            label_period=None,
            time=data["time"], flux=data["flux"], quality=data["quality"],
            stellar=to_stellar(p), white_ppm=white_ppm(data["flux"]),
            meta_extra={"grp": "quiet"},
        )
        r = RC.evaluate_target(target, run_dir)
        r["kepid"] = kepid
        results_false.append(r)
        json.dump(results_false, open(run_dir / "results_false.json", "w"), indent=2)
        logger.info(f"[false {i}/{len(false_picks)}] KIC{kepid} -> {r['validation_status']} "
                    f"snr={r.get('snr')}")

    for r in results_true + results_false:
        r.setdefault("period_error_pct", None)
        r.setdefault("detected", False)
        r.setdefault("validation_status", "SKIPPED")
    rmetrics = real_metrics(results_true, results_false)
    write_report(run_dir, results_true, results_false, None, rmetrics, started,
                 true_picks, false_picks)

    print("\n" + "=" * 60)
    print("REAL-DATA HONESTY EVALUATION (Kepler / NEA truth)")
    print("=" * 60)
    print(f"Recall@target period   : {rmetrics['recall_target']*100:.1f}%  ({rmetrics['tp_target']}/{rmetrics['n_true']})")
    print(f"Recall@any known planet: {rmetrics['recall_known']*100:.1f}%  ({(rmetrics['tp_target']+rmetrics['tp_other'])}/{rmetrics['n_true']})")
    print(f"Quiet-star FPR (proxy) : {rmetrics['fpr_quiet']*100:.1f}%  ({rmetrics['fp_quiet']} FP / {rmetrics['n_false']})")
    print(f"Total FP rate          : {rmetrics['fpr_total']*100:.1f}%  ({rmetrics['fp']} FP / {rmetrics['fp']+rmetrics['tn']})")
    print(f"Precision (target)     : {rmetrics['precision']*100:.1f}%")
    print(f"Report: {run_dir / 'EVALUATION_REPORT.md'}")


def write_report(run_dir: Path, true_res, false_res, metrics, rmetrics, started_at: str,
                 true_picks, false_picks) -> None:
    CERT = RC.CERTIFIED
    lines = ["# Axiom-ZSpace Real-Data Honesty Evaluation (Kepler / NEA)", ""]
    lines.append(f"**Generated:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}")
    lines.append(f"**Run:** {run_dir.name}")
    lines.append(f"**Started:** {started_at}")
    lines.append("**Method:** identical blind pipeline to the controlled bench (BLS +")
    lines.append("ladder + auditors + sovereign validator, archive stub OFFLINE); input =")
    lines.append("real Kepler long-cadence light curves (Q4-Q9 ~540 d baseline).")
    lines += ["", "## Key Metrics", "",
              "| Metric | Value |", "|---|---|",
              f"| Confirmed planets (NEA truth, 1.0-13.5 d) | {rmetrics['n_true']} |",
              f"| Quiet stars (zero KOI proxy) | {rmetrics['n_false']} |",
              f"| **Recall@target period** (certified & P within 5% of NEA period) | {rmetrics['recall_target']*100:.1f}% ({rmetrics['tp_target']}/{rmetrics['n_true']}) |",
              f"| **Recall@any known planet of host** (certified & P matches another known planet) | {rmetrics['recall_known']*100:.1f}% ({(rmetrics['tp_target']+rmetrics['tp_other'])}/{rmetrics['n_true']}) |",
              f"| **Quiet-star certification rate** (proxy FPR) | {rmetrics['fpr_quiet']*100:.1f}% ({rmetrics['fp_quiet']}/{rmetrics['n_false']}) |",
              f"| **Total FP rate** (incl. wrong-ephemeris certs on true hosts) | {rmetrics['fpr_total']*100:.1f}% ({rmetrics['fp']} FP / {rmetrics['fp']+rmetrics['tn']}) |",
              f"| **Precision (target-level)** | {rmetrics['precision']*100:.1f}% |",
              f"| Confusion | TP(target)={rmetrics['tp_target']} TP(other known)={rmetrics['tp_other']} FP(quiet)={rmetrics['fp_quiet']} FP(unknown)={rmetrics['fp_unknown']} TN={rmetrics['tn']} |",
              "",
              "*Controlled-bench comparison (for context): recall 52%, FPR 0% on",
              "  synthetic noise/contamination targets. Real data is harder:",
              "  long cadence, real systematics, and multiple real signals per",
              "  system (e.g. Kepler-37's outer planets alias into the band).*",
              ""]
    lines += ["## Target list (true)", "",
              "| KIC | star | NEA P (d) | R_p (Re) | found P | perr% | match | status |",
              "|---|---|---|---|---|---|---|---|"]
    for r, p in zip(true_res, true_picks):
        perr = r.get("period_error_pct")
        perr_s = f"{perr:.1f}" if perr is not None else "-"
        mt = "TARGET" if r.get("matched_target") else (r.get("match_known") or "-")
        lines.append(f"| {r.get('kepid')} | {p.get('kepler_name') or ''} | {p['pl_orbper']} | {p.get('pl_rade')} | "
                     f"{r.get('period_found')} | {perr_s} | {mt} | {r.get('validation_status')} |")
    lines += ["", "## Target list (false / quiet)", "", "| KIC | found P | snr | status |", "|---|---|---|---|"]
    for r, p in zip(false_res, false_picks):
        lines.append(f"| {r.get('kepid')} | {r.get('period_found')} | {r.get('snr')} | {r.get('validation_status')} |")
    lines += ["", "## Known signals per true host (NEA cumulative)", "", "| KIC | P (d) | disp | name |", "|---|---|---|---|"]
    for r in true_res:
        for k in r.get("known_rows", []):
            lines.append(f"| {r.get('kepid')} | {k['koi_period']} | {k['koi_pdisposition']} | {k.get('kepler_name') or ''} |")
    lines += ["", "## Honesty caveats (REQUIRED READING)", "",
              "1. **The quiet-star FPR is a proxy.** 'Quiet' = zero KOIs / zero",
              "   confirmed planets in DR25. Undiscovered planets or undetected",
              "   eclipsing binaries below the Kepler detection limit may still",
              "   exist; a certification on such a star is flagged FP but is not",
              "   proven to be noise.",
              "2. **Truth = NEA period (5% window).** A certification at a period",
              "   matching ANOTHER known planet of the same host is scored as",
              "   'other known' (real signal, different planet), not as a false",
              "   positive. A certification matching a cataloged FALSE POSITIVE",
              "   (e.g. an EB) is scored as FP — correctly.",
              "3. **Band aliases:** signals with true period >13.5 d can appear as",
              "   in-band sub-harmonics (e.g. Kepler-37 d: P=39.8 d, 5th harmonic",
              "   at 7.96 d certified). The alias resolver tests 2x/3x only.",
              "4. Long cadence (29.4 min) + quarter gaps + real systematics make",
              "   real-data recall lower than the controlled synthetic recall.",
              "5. Periods were restricted to the pipeline search band (0.5-13.5 d);",
              "   the 5% window is an acceptance rule, not a precision claim.",
              ""]
    (run_dir / "EVALUATION_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()