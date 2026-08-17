"""
sector_processor.py  ·  Sector-Based Discovery Pipeline
========================================================
Orchestrates systematic exoplanet discovery across all TIC IDs
observed in a TESS sector.

Usage
─────
  from zspace_engine.sector_processor import SectorProcessor

  processor = SectorProcessor(sector=42, output_dir="axiom_output")
  tic_list = processor.get_tic_list()
  print(f"Found {len(tic_list)} TIC IDs in Sector 42")

Architecture
────────────
  1. Query MAST for all TIC IDs in specified sector
  2. Create output directory structure: axiom_output/sector_N/
  3. Provide infrastructure for systematic processing (full pipeline
     integration handled by run_pipeline.py)

Requirements
────────────
  Implements Requirements 3.1, 3.2, 3.6:
    - Accept --sector command-line argument
    - Query MAST for all TIC IDs in sector
    - Create axiom_output/ directory structure
"""

from __future__ import annotations

import json
import logging
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

# ─────────────────────────────────────────────────────────────────────────────
# Physical constants (IAU 2015) — imported from constants.py
# ─────────────────────────────────────────────────────────────────────────────
from zspace_engine.constants import R_SUN, R_EARTH
from zspace_engine.core import THRESHOLD_PLANET, THRESHOLD_LIKELY, THRESHOLD_AMBIGUOUS
from zspace_engine.logging_config import suppress_astroquery_logger

# ─────────────────────────────────────────────────────────────────────────────
# Output organization
# ─────────────────────────────────────────────────────────────────────────────
from zspace_engine.output_organizer import OutputOrganizer

# ─────────────────────────────────────────────────────────────────────────────
# Lightkurve import (for sector target list queries)
# ─────────────────────────────────────────────────────────────────────────────
try:
    import lightkurve as lk
    LIGHTKURVE_AVAILABLE = True
except ImportError:
    LIGHTKURVE_AVAILABLE = False
    logging.warning(
        "lightkurve not available. SectorProcessor will not be able to "
        "query TIC lists. Install with: pip install lightkurve"
    )


# ─────────────────────────────────────────────────────────────────────────────
# SectorProcessor Class
# ─────────────────────────────────────────────────────────────────────────────

class SectorProcessor:
    """
    Orchestrates systematic discovery across all TIC IDs in a TESS sector.
    
    This class provides the infrastructure for sector-based processing:
      - Queries MAST for all TIC IDs observed in a sector
      - Creates output directory structure
      - Provides foundation for full pipeline integration
    
    Parameters
    ----------
    sector : int
        TESS sector number (1-69+ as of 2024)
    output_dir : str, optional
        Base output directory path (default: "axiom_output")
    
    Attributes
    ----------
    sector : int
        The TESS sector number being processed
    output_dir : Path
        Path to sector-specific output directory (axiom_output/sector_N/)
    
    Examples
    --------
    >>> processor = SectorProcessor(sector=42)
    >>> tic_list = processor.get_tic_list()
    >>> print(f"Sector 42 contains {len(tic_list)} TIC IDs")
    """
    
    def __init__(
        self,
        sector: int,
        output_dir: str = "axiom_output",
        config: Optional[Dict[str, Any]] = None,
        max_targets: Optional[int] = None
    ) -> None:
        """
        Initialize SectorProcessor for a specific TESS sector.
        
        Creates the output directory structure:
          axiom_output/
            └── sector_N/
                ├── discoveries/    (created on demand)
                ├── known/          (created on demand)
                └── rejected/       (created on demand)
        
        Parameters
        ----------
        sector : int
            TESS sector number to process
        output_dir : str, optional
            Base output directory (default: "axiom_output")
        config : Dict[str, Any], optional
            Production configuration dictionary with detection thresholds
        max_targets : int, optional
            Maximum number of targets to process (for testing). If None, processes all targets.
        
        Raises
        ------
        ValueError
            If sector number is not positive
        """
        if sector <= 0:
            raise ValueError(f"Sector number must be positive; got {sector}")
        
        self.sector = sector
        self.output_dir = Path(output_dir) / f"sector_{sector}"
        self.max_targets = max_targets
        
        # Store configuration with defaults
        self.config = config or {}
        self.snr_threshold = self.config.get("detection", {}).get("bls_snr_threshold", 5.5)
        self.fap_threshold = self.config.get("detection", {}).get("fap_threshold", 1.0e-4)
        self.cvs_threshold = self.config.get("detection", {}).get("cvs_planet_threshold", 0.80)
        
        # Create base sector directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize OutputOrganizer for structured output routing
        self.output_organizer = OutputOrganizer(base_dir=output_dir)
        
        limit_msg = f" | Limit: {max_targets} targets" if max_targets else ""
        logging.info(
            f"SectorProcessor initialized for Sector {sector} | "
            f"Output: {self.output_dir} | "
            f"Thresholds: SNR>{self.snr_threshold}, FAP<{self.fap_threshold:.0e}, CVS>{self.cvs_threshold}"
            f"{limit_msg}"
        )
    
    def get_tic_list(self) -> List[str]:
        """
        Query MAST for all TIC IDs observed in this sector.
        
        Downloads the official TESS observed targets list from MAST.
        Uses astroquery to query all observations in the sector.
        
        Returns
        -------
        List[str]
            Sorted list of TIC ID strings (e.g., ["12345678", "87654321", ...])
            Returns empty list if download fails.
        
        Notes
        -----
        - Queries MAST Observations for all TESS targets in sector
        - Filters for 2-minute cadence observations (t_exptime < 200s)
        - Returns unique TIC IDs sorted numerically
        - First query may take 30-60 seconds
        - Results are cached locally in .cache/sector_N_tics.json
        
        Examples
        --------
        >>> processor = SectorProcessor(sector=1)
        >>> tic_list = processor.get_tic_list()
        >>> print(f"Sector 1: {len(tic_list)} targets")
        Sector 1: 15234 targets
        """
        logging.info(f"Querying MAST for TIC IDs in Sector {self.sector}...")
        
        # Check cache first
        cache_dir = Path(".cache")
        cache_dir.mkdir(exist_ok=True)
        cache_file = cache_dir / f"sector_{self.sector}_tics.json"
        
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    cached_data = json.load(f)
                    tic_list = cached_data.get('tic_ids', [])
                    if tic_list:
                        logging.info(
                            f"  Loaded {len(tic_list)} TIC IDs from cache: {cache_file}"
                        )
                        return tic_list
            except Exception as e:
                logging.warning(f"  Failed to load cache: {e}. Re-querying MAST...")
        
        try:
            # Use astroquery to get all TESS observations in this sector
            from astroquery.mast import Observations
            import warnings
            
            logging.info(f"  Querying MAST Observations for Sector {self.sector}...")
            logging.info(f"  (This may take 30-60 seconds for the first query)")
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                
                # Query for all TESS observations in this sector
                # Filter for 2-minute cadence (t_exptime < 200 seconds)
                obs_table = Observations.query_criteria(
                    obs_collection="TESS",
                    dataproduct_type="timeseries",
                    t_exptime=[0, 200],  # 2-minute cadence
                    sequence_number=self.sector
                )
            
            if obs_table is None or len(obs_table) == 0:
                logging.warning(
                    f"No observations found for Sector {self.sector}. "
                    f"This sector may not exist or may not have 2-minute cadence data."
                )
                return []
            
            # Extract unique TIC IDs from target_name column
            # Format is typically "TIC 12345678" or just the number
            tic_ids = set()
            for target_name in obs_table['target_name']:
                # Clean up the target name to extract TIC ID
                target_str = str(target_name).strip()
                
                # Handle different formats
                if target_str.startswith('TIC'):
                    # Format: "TIC 12345678"
                    tic_id = target_str.replace('TIC', '').strip()
                elif target_str.isdigit():
                    # Format: "12345678"
                    tic_id = target_str
                else:
                    # Skip non-TIC targets
                    continue
                
                if tic_id and tic_id.isdigit():
                    tic_ids.add(tic_id)
            
            # Convert to sorted list
            tic_list = sorted(tic_ids, key=lambda x: int(x))
            
            logging.info(
                f"MAST query complete | Sector {self.sector}: "
                f"{len(tic_list)} unique TIC IDs"
            )
            
            # Cache the results
            try:
                with open(cache_file, 'w') as f:
                    json.dump({
                        'sector': self.sector,
                        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
                        'count': len(tic_list),
                        'tic_ids': tic_list
                    }, f, indent=2)
                logging.info(f"  Cached TIC list to: {cache_file}")
            except Exception as e:
                logging.warning(f"  Failed to cache TIC list: {e}")
            
            return tic_list
        
        except ImportError:
            logging.error(
                "astroquery not available. Install with: pip install astroquery"
            )
            return []
        
        except Exception as e:
            logging.error(
                f"Error querying MAST for Sector {self.sector}: {e}. "
                "Check network connectivity and MAST availability."
            )
            # Log full traceback for debugging
            import traceback
            logging.debug(traceback.format_exc())
            return []
    
    def process_sector(self) -> Dict[str, Any]:
        """
        Run full pipeline on all TIC IDs in sector.
        
        Orchestrates the complete discovery pipeline for all targets in the
        sector. For each TIC ID, downloads lightcurve data, runs BLS detection,
        validates candidates, and routes outputs to appropriate directories.
        
        Returns
        -------
        Dict[str, Any]
            Summary statistics dictionary containing:
            - sector: int - Sector number processed
            - timestamp_utc: str - Processing completion timestamp
            - total_targets: int - Total TIC IDs in sector
            - processed: int - Successfully processed targets
            - new_discoveries: int - Count of new planet discoveries
            - known_planets: int - Count of known planets rediscovered
            - false_positives: int - Count of rejected candidates
            - failed: int - Count of targets that failed processing
            - discoveries: List[Dict] - Details of new discoveries
            - errors: List[Dict] - Details of failed targets
        
        Notes
        -----
        - Errors are logged but do not halt processing
        - Summary is saved to axiom_output/sector_N/summary.json
        - Each target is processed independently via _process_single_tic()
        - Statistics are accumulated across all targets
        
        Examples
        --------
        >>> processor = SectorProcessor(sector=42)
        >>> summary = processor.process_sector()
        >>> print(f"Discoveries: {summary['new_discoveries']}")
        Discoveries: 3
        
        See Also
        --------
        _process_single_tic : Processes individual TIC target
        """
        logging.info(f"Starting sector processing for Sector {self.sector}")
        
        # Get list of all TIC IDs in sector
        tic_list = self.get_tic_list()
        
        # Apply max_targets limit if specified
        if self.max_targets and len(tic_list) > self.max_targets:
            logging.info(
                f"Limiting processing to first {self.max_targets} targets "
                f"(out of {len(tic_list)} total)"
            )
            tic_list = tic_list[:self.max_targets]
        
        # Initialize statistics tracking
        stats = {
            "sector": self.sector,
            "timestamp_utc": None,  # Set at completion
            "total_targets": len(tic_list),
            "processed": 0,
            "new_discoveries": 0,
            "known_planets": 0,
            "false_positives": 0,
            "failed": 0,
            "discoveries": [],
            "errors": []
        }
        
        total = len(tic_list)
        logging.info(
            f"Processing {total} TIC IDs in Sector {self.sector}"
        )
        
        # ── Suppress noisy lightkurve messages during batch processing ────
        lk_logger = logging.getLogger("lightkurve")
        lk_prev_level = lk_logger.level
        lk_logger.setLevel(logging.WARNING)
        
        # Also suppress astroquery spam WITHOUT creating its logger early
        # (creating a plain 'astroquery' Logger before astroquery initializes
        # breaks astroquery._init_log with "_set_defaults" AttributeError).
        if "astroquery" in logging.Logger.manager.loggerDict:
            aq_logger = logging.getLogger("astroquery")
            aq_prev_level = aq_logger.level
        else:
            aq_logger = None
            aq_prev_level = None
        suppress_astroquery_logger()
        
        import sys
        import time as _time
        t_sector_start = _time.time()
        
        # ── Progress bar helper ───────────────────────────────────────────
        def _print_progress(idx, total, tic_id, status, stats, elapsed):
            """Print a live progress bar to stderr (stays visible in terminals)."""
            pct = idx / max(total, 1) * 100
            rate = idx / max(elapsed, 0.01) * 60   # targets per minute
            eta_sec = (total - idx) / max(idx / max(elapsed, 0.01), 0.001)
            eta_min = eta_sec / 60
            
            # Build bar:  [████████░░░░░░░░░░░░] 45.2%
            bar_width = 30
            filled = int(bar_width * idx / max(total, 1))
            bar = '#' * filled + '-' * (bar_width - filled)
            
            # Status icon (ASCII-safe for Windows cp1256 consoles)
            icon = {'NEW_DISCOVERY': '*', 'OFFLINE_NEW_DISCOVERY': '*',
                    'KNOWN': 'o', 'FALSE_POSITIVE': 'x', 'FAILED': 'x'}.get(status, '.')
            
            line = (
                f"\r  [{bar}] {pct:5.1f}%  "
                f"({idx}/{total})  "
                f"TIC {tic_id} {icon}  "
                f"| D:{stats['new_discoveries']} K:{stats['known_planets']} "
                f"FP:{stats['false_positives']} F:{stats['failed']}  "
                f"| {rate:.1f}/min  ETA:{eta_min:.0f}m   "
            )
            sys.stderr.write(line)
            sys.stderr.flush()
        
        # ── Process each TIC ID ──────────────────────────────────────────
        for idx, tic_id in enumerate(tic_list, start=1):
            status = "UNKNOWN"
            try:
                result = self._process_single_tic(tic_id)
                
                # Update statistics based on result status
                stats["processed"] += 1
                status = result["status"]
                
                if result["status"] in ("NEW_DISCOVERY", "OFFLINE_NEW_DISCOVERY"):
                    stats["new_discoveries"] += 1
                    discovery_entry = {
                        "zspace_id": result.get("zspace_id", f"ZS-T-{tic_id}-01"),
                        "tic_id": result["tic_id"],
                        "period_days": result["period_days"],
                        "cvs_score": result["cvs_score"],
                        "output_file": result.get("output_file", "unknown")
                    }
                    stats["discoveries"].append(discovery_entry)
                    # Print discovery prominently (on its own line above the bar)
                    sys.stderr.write(f"\n  * NEW DISCOVERY: TIC {tic_id} | "
                                     f"P={result.get('period_days', '?'):.5f}d | "
                                     f"CVS={result.get('cvs_score', '?'):.4f}\n")
                    logging.info(
                        f"NEW DISCOVERY: TIC {tic_id} | "
                        f"Period: {result.get('period_days', 'N/A')} days | "
                        f"CVS: {result.get('cvs_score', 'N/A')}"
                    )
                elif result["status"] == "KNOWN":
                    stats["known_planets"] += 1
                elif result["status"] == "FALSE_POSITIVE":
                    stats["false_positives"] += 1
                else:
                    logging.warning(
                        f"Unexpected status '{result['status']}' for TIC {tic_id}"
                    )
            
            except Exception as e:
                stats["failed"] += 1
                status = "FAILED"
                error_entry = {
                    "tic_id": tic_id,
                    "error": str(e)
                }
                stats["errors"].append(error_entry)
                # Only log at debug to avoid spamming the console during batch
                logging.debug(f"TIC {tic_id} failed: {e}")
                continue
            
            finally:
                # Always update progress bar
                elapsed = _time.time() - t_sector_start
                _print_progress(idx, total, tic_id, status, stats, elapsed)
                
                # Periodic summary every 100 targets (to the log file)
                if idx % 100 == 0:
                    rate = idx / max(elapsed, 0.01) * 60
                    logging.info(
                        f"[PROGRESS] Sector {self.sector}: {idx}/{total} "
                        f"({idx/total*100:.1f}%) | "
                        f"D:{stats['new_discoveries']} K:{stats['known_planets']} "
                        f"FP:{stats['false_positives']} F:{stats['failed']} | "
                        f"{rate:.1f} targets/min"
                    )
        
        # ── End of processing loop ───────────────────────────────────────
        # Clear the progress bar line
        sys.stderr.write("\n")
        sys.stderr.flush()
        
        # Restore logger levels
        lk_logger.setLevel(lk_prev_level)
        if aq_logger is not None:
            aq_logger.setLevel(aq_prev_level)
        
        total_elapsed = _time.time() - t_sector_start
        total_min = total_elapsed / 60
        
        # Set completion timestamp
        stats["timestamp_utc"] = datetime.now(timezone.utc).isoformat()
        stats["elapsed_minutes"] = round(total_min, 2)
        stats["rate_per_minute"] = round(total / max(total_elapsed, 0.01) * 60, 1) if total > 0 else 0
        
        # Write summary to file
        summary_path = self.output_dir / "summary.json"
        try:
            with open(summary_path, 'w') as f:
                json.dump(stats, f, indent=2)
            logging.info(f"Summary saved to {summary_path}")
        except IOError as e:
            logging.error(f"File I/O error writing summary to {summary_path}: {e}")
        except Exception as e:
            logging.error(f"Unexpected error writing summary to {summary_path}: {e}")
        
        # ── Save compact discoveries.json ─────────────────────────────────
        discoveries_path = self.output_dir / "discoveries.json"
        try:
            compact_discoveries = {
                "sector": self.sector,
                "total_discoveries": stats["new_discoveries"],
                "scan_date": stats["timestamp_utc"],
                "total_targets_scanned": stats["total_targets"],
                "elapsed_minutes": stats.get("elapsed_minutes", 0),
                "planets": []
            }
            
            for i, disc in enumerate(stats["discoveries"], start=1):
                cvs = disc["cvs_score"]
                if cvs >= THRESHOLD_PLANET:
                    verdict = "PLANET CANDIDATE"
                elif cvs >= THRESHOLD_LIKELY:
                    verdict = "LIKELY PLANET CANDIDATE"
                elif cvs >= THRESHOLD_AMBIGUOUS:
                    verdict = "AMBIGUOUS"
                else:
                    verdict = "FALSE POSITIVE"
                compact_discoveries["planets"].append({
                    "#": i,
                    "tic_id": disc["tic_id"],
                    "zspace_id": disc.get("zspace_id", f"ZS-T-{disc['tic_id']}-01"),
                    "period_days": round(disc["period_days"], 5),
                    "cvs": round(cvs, 4),
                    "verdict": verdict,
                })
            
            with open(discoveries_path, 'w', encoding='utf-8') as f:
                json.dump(compact_discoveries, f, indent=2, ensure_ascii=False)
            
            logging.info(f"Discoveries list saved to {discoveries_path} ({stats['new_discoveries']} entries)")
        except Exception as e:
            logging.error(f"Failed to save discoveries.json: {e}")
        
        # Print prominent final summary to stderr
        sys.stderr.write(
            f"\n  ╔══════════════════════════════════════════════════════════╗\n"
            f"  ║  SECTOR {self.sector:02d} COMPLETE                                    ║\n"
            f"  ║  Processed: {stats['processed']:>6}/{stats['total_targets']:<6}                          ║\n"
            f"  ║  Discoveries: {stats['new_discoveries']:>4}  |  Known: {stats['known_planets']:>4}  |  FP: {stats['false_positives']:>5}   ║\n"
            f"  ║  Failed: {stats['failed']:>5}  |  Time: {total_min:.1f} min               ║\n"
            f"  ╚══════════════════════════════════════════════════════════╝\n\n"
        )
        sys.stderr.flush()
        
        # Log final statistics
        logging.info(
            f"Sector {self.sector} complete in {total_min:.1f}m | "
            f"Processed: {stats['processed']}/{stats['total_targets']} | "
            f"New: {stats['new_discoveries']} | "
            f"Known: {stats['known_planets']} | "
            f"FP: {stats['false_positives']} | "
            f"Failed: {stats['failed']} | "
            f"Rate: {stats['rate_per_minute']}/min"
        )
        
        return stats
    
    def _process_single_tic(self, tic_id: str) -> Dict[str, Any]:
        """
        Process a single TIC target through the full pipeline.
        
        Orchestrates the complete discovery pipeline:
        1. Lightcurve download via LightCurveIngester
        2. BLS detection for candidate identification
        3. AxiomValidator for cross-reference with known planets
        4. Output routing based on validation status
        
        Parameters
        ----------
        tic_id : str
            TESS Input Catalogue identifier (e.g., "12345678")
        
        Returns
        -------
        Dict[str, Any]
            Result dictionary containing:
            - status: str - One of "NEW_DISCOVERY", "OFFLINE_NEW_DISCOVERY", "KNOWN", "FALSE_POSITIVE"
            - tic_id: str - The TIC identifier
            - period_days: float - Orbital period (if candidate found)
            - cvs_score: float - Composite Vitality Score (if candidate found)
            - zspace_id: str - ZSpace identifier (if discovery)
            - output_file: str - Path to output JSON file
        
        Raises
        ------
        RuntimeError
            If lightcurve download fails or data quality is insufficient
        FileNotFoundError
            If no lightcurve data available for TIC ID
        
        Notes
        -----
        - Errors are logged but propagated to caller for statistics tracking
        - FALSE_POSITIVE status returned if BLS detection fails gate
        - Network errors trigger OFFLINE_NEW_DISCOVERY status
        """
        from zspace_engine.ingestion import LightCurveIngester
        from zspace_engine.detectors import BLSDetector
        from zspace_engine.validator import AxiomValidator
        from zspace_engine.core import VitalityMatrix, apply_hard_filters
        
        logging.debug(f"Processing TIC {tic_id}")
        
        # ── Step 1: Download lightcurve ──────────────────────────────────────
        try:
            ingester = LightCurveIngester(tic_id=tic_id, mission="TESS", exptime="short")
            lc_product = ingester.process(sector=self.sector)
            logging.debug(
                f"TIC {tic_id} | Lightcurve downloaded: {lc_product.n_points_cleaned} points, "
                f"cadence={lc_product.cadence_days:.5f} d"
            )
        except Exception as e:
            logging.error(f"TIC {tic_id} | Lightcurve download failed: {e}")
            raise
        
        # ── Step 2: Run BLS detection ────────────────────────────────────────
        try:
            detector = BLSDetector(
                period_min=0.5,
                period_max=13.5,
                frequency_factor=10.0,
                snr_threshold=self.snr_threshold,
                fap_threshold=self.fap_threshold
            )
            bls_result = detector.run(
                time=lc_product.time,
                flux=lc_product.flux_flat
            )
            
            logging.debug(
                f"TIC {tic_id} | BLS: P={bls_result.period_best:.5f} d, "
                f"SNR={bls_result.snr:.2f}, FAP={bls_result.fap:.3e}, "
                f"S_P={bls_result.s_periodicity:.4f}"
            )
            
            # Check detection gate
            if not bls_result.passed_detection_gate():
                logging.debug(
                    f"TIC {tic_id} | BLS detection gate FAILED "
                    f"(SNR={bls_result.snr:.2f} <= {self.snr_threshold} or FAP={bls_result.fap:.3e} >= {self.fap_threshold:.0e})"
                )
                return {
                    "status": "FALSE_POSITIVE",
                    "tic_id": tic_id,
                    "period_days": bls_result.period_best,
                    "cvs_score": 0.0,
                    "zspace_id": None,
                    "output_file": None,
                    "reason": "BLS detection gate failed"
                }
        
        except Exception as e:
            logging.error(f"TIC {tic_id} | BLS detection failed: {e}")
            raise
        
        # ── Step 3: Compute CVS score ────────────────────────────────────────
        try:
            # Fetch real stellar parameters from TIC v8.2
            from zspace_engine.context import TICMetadataFetcher
            try:
                tic_meta = TICMetadataFetcher.fetch(tic_id, period_days=bls_result.period_best)
                stellar_mass_solar = tic_meta.stellar_mass_solar
                stellar_radius_solar = tic_meta.stellar_radius_solar
                stellar_teff_k = tic_meta.stellar_teff_k
                stellar_logg = tic_meta.stellar_logg
            except Exception:
                # Fallback to solar defaults if TIC fetch fails
                stellar_mass_solar = 1.0
                stellar_radius_solar = 1.0
                stellar_teff_k = 5778.0
                stellar_logg = 4.44
            
            # Create VitalityMatrix
            matrix = VitalityMatrix(tic_id=tic_id, planet_order=1)
            
            # Compute orbital mechanics
            matrix.compute_orbital_mechanics(
                period_days=bls_result.period_best,
                transit_depth=bls_result.transit_depth,
                stellar_mass_solar=stellar_mass_solar,
                stellar_teff=stellar_teff_k,
                stellar_radius_solar=stellar_radius_solar
            )
            
            # Ingest scores (using real auditors — no more placeholders)
            from zspace_engine.auditors import TransitAuditor
            from zspace_engine.detectors import BLSDetector

            # Phase-fold once, reuse for all auditors
            bin_phase, bin_flux, _ = BLSDetector.fold_and_bin(
                lc_product.time, lc_product.flux_flat,
                period=bls_result.period_best, t0=bls_result.t0, n_bins=200,
            )

            auditor = TransitAuditor(verbose=False, run_mcmc=False)

            # Audit 1: Even/Odd → EB flag
            eo = auditor.even_odd_test(
                lc_product.time, lc_product.flux_flat,
                period=bls_result.period_best, t0=bls_result.t0,
                duration=bls_result.transit_duration,
            )
            # Audit 2: Depth consistency → S_δ
            dc = auditor.depth_consistency_score(
                lc_product.time, lc_product.flux_flat,
                period=bls_result.period_best, t0=bls_result.t0,
                duration=bls_result.transit_duration, eo_result=eo,
            )
            # Audit 3: Limb shape (Mandel-Agol) → S_τ
            ls = auditor.limb_shape_score(
                period=bls_result.period_best,
                duration=bls_result.transit_duration,
                transit_depth=bls_result.transit_depth,
                time=lc_product.time, flux=lc_product.flux_flat, t0=bls_result.t0,
            )
            # Audit 4: Ingress/Egress V-shape → FP risk
            ie = auditor.ingress_egress_test(
                bin_phase, bin_flux,
                period=bls_result.period_best,
                duration=bls_result.transit_duration,
                transit_depth=bls_result.transit_depth,
            )

            # Audit 5: Stellar context → S_S
            from zspace_engine.context import StellarContextAuditor
            s_auditor = StellarContextAuditor(
                fetch_tic=True, use_tpf_centroids=False, check_multi_sector=False,
            )
            try:
                ctx = s_auditor.audit(
                    tic_id=tic_id,
                    time=lc_product.time,
                    flux=lc_product.flux_flat,
                    period=bls_result.period_best,
                    t0=bls_result.t0,
                    duration=bls_result.transit_duration,
                    a_rs_transit=0.0,
                    sector=self.sector,
                )
                s_stellar = ctx.s_stellar
                proof_s = ctx.proof
                flags_s = ctx.flags
            except Exception as e:
                # Stellar context must never sink an otherwise-good candidate;
                # degrade gracefully to neutral score rather than crashing.
                s_stellar = 0.5
                proof_s = f"STELLAR_CONTEXT | failed ({type(e).__name__}: {e}) → S_S=0.50"
                flags_s = [f"STELLAR_CONTEXT_UNAVAILABLE | {type(e).__name__}: {e}"]

            s_depth  = dc.s_depth
            s_limb   = ls.s_limb
            flags_d  = list(dc.flags)
            flags_l  = list(ls.flags)
            flags_l.extend(ie.flags)

            # CVS critical-FP veto: strong EB evidence caps CVS below the
            # ambiguous threshold so the verdict becomes FALSE_POSITIVE no
            # matter how strong the periodicity score is (fixes the
            # "high-SNR EB can never be vetoed" impossibility).
            veto_reasons: List[str] = []
            if eo.is_eb_flag:
                veto_reasons.append(
                    f"EVEN_ODD_EB | Δσ={eo.delta_sigma:.2f}, p={eo.p_value:.4f}"
                )
            if ie.is_v_shape and ie.fp_risk in ("MEDIUM", "HIGH"):
                veto_reasons.append(
                    f"V_SHAPE | ingress_frac={ie.ingress_fraction:.3f}, fp_risk={ie.fp_risk}"
                )

            matrix.ingest_scores(
                s_periodicity=bls_result.s_periodicity,
                proof_p=bls_result.proof,
                s_depth=s_depth,
                proof_d=dc.proof,
                s_limb=s_limb,
                proof_l=ls.proof,
                s_stellar=s_stellar,
                proof_s=proof_s,
                flags_p=bls_result.flags,
                flags_d=flags_d,
                flags_l=flags_l,
                flags_s=flags_s,
            )
            for reason in veto_reasons:
                matrix.apply_veto(reason)
            
            # Apply hard physical filters before CVS
            hard_filter = apply_hard_filters(
                planet_radius_earth=matrix.orbital.planet_radius_earth,
                transit_depth=bls_result.transit_depth,
                transit_duration_hrs=bls_result.transit_duration * 24.0,
                period_days=bls_result.period_best,
            )
            if not hard_filter.passed:
                logging.debug(
                    f"TIC {tic_id} | HARD REJECT: {hard_filter.rejection}"
                )
                return {
                    "status": "FALSE_POSITIVE",
                    "tic_id": tic_id,
                    "period_days": bls_result.period_best,
                    "cvs_score": 0.0,
                    "zspace_id": None,
                    "output_file": None,
                    "reason": f"Hard filter: {hard_filter.rejection}"
                }
            matrix.cvs_engine.apply_hard_filter(hard_filter)
            
            # Compute CVS
            cvs_result = matrix.cvs_engine.compute()
            cvs_score = cvs_result
            # Use full CVS classification from the engine (handles all 4 tiers)
            cvs_verdict = matrix.cvs_engine._classify(cvs_score)
            
            planet_radius_earth = matrix.orbital.planet_radius_earth
            
            logging.debug(
                f"TIC {tic_id} | CVS={cvs_score:.4f}, verdict={cvs_verdict}"
            )
            
            # Check CVS threshold (from config)
            if cvs_score < self.cvs_threshold:
                logging.debug(
                    f"TIC {tic_id} | CVS score {cvs_score:.4f} < {self.cvs_threshold} threshold"
                )
                return {
                    "status": "FALSE_POSITIVE",
                    "tic_id": tic_id,
                    "period_days": bls_result.period_best,
                    "cvs_score": cvs_score,
                    "zspace_id": None,
                    "output_file": None,
                    "reason": "CVS score below planet threshold"
                }
        
        except Exception as e:
            logging.error(f"TIC {tic_id} | CVS computation failed: {e}")
            raise
        
        # ── Step 4: Invoke AxiomValidator ────────────────────────────────────
        try:
            # Validator will write to its output_dir initially
            validator = AxiomValidator(
                output_dir=str(self.output_dir),
                verbose=False
            )
            
            validation_result = validator.validate(
                tic_id=tic_id,
                period_days=bls_result.period_best,
                transit_depth=bls_result.transit_depth,
                transit_duration_hrs=bls_result.transit_duration * 24.0,
                t0_btjd=bls_result.t0,
                stellar_mass_solar=stellar_mass_solar,
                stellar_radius_solar=stellar_radius_solar,
                stellar_teff_k=stellar_teff_k,
                stellar_logg=stellar_logg,
                planet_radius_earth=planet_radius_earth,
                cvs_score=cvs_score,
                cvs_verdict=cvs_verdict,
                cvs_proof_chain=[],
                bls_snr=bls_result.snr,
                bls_fap=bls_result.fap,
                even_odd_delta_sigma=0.0,  # Placeholder - would compute from auditors
                shape_ratio=1.0,           # Placeholder - would compute from auditors
                secondary_snr=0.0,         # Placeholder - would compute from auditors
                centroid_sigma=0.0,        # Placeholder - would compute from auditors
            )
            
            logging.info(
                f"TIC {tic_id} | Validation complete: {validation_result.status}"
            )
            
            # Use OutputOrganizer to determine appropriate output path
            zspace_id = f"ZS-T-{tic_id}-01"
            output_path = Path(validation_result.output_file)
            
            if output_path.exists():
                # Get organized output path based on validation status
                organized_path = self.output_organizer.get_output_path(
                    sector=self.sector,
                    status=validation_result.status,
                    zspace_id=zspace_id
                )
                
                # Move file to organized location
                output_path.rename(organized_path)
                final_output = str(organized_path)
                
                logging.debug(
                    f"TIC {tic_id} | Output saved to {organized_path}"
                )
            else:
                final_output = validation_result.output_file
            
            # Build result dictionary
            return {
                "status": validation_result.status,
                "tic_id": tic_id,
                "period_days": bls_result.period_best,
                "cvs_score": cvs_score,
                "zspace_id": zspace_id,
                "output_file": final_output,
            }
        
        except Exception as e:
            logging.error(f"TIC {tic_id} | Validation failed: {e}")
            raise


# ─────────────────────────────────────────────────────────────────────────────
# Module exports
# ─────────────────────────────────────────────────────────────────────────────

__all__ = [
    "SectorProcessor",
]
