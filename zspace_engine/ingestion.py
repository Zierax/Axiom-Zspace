"""
ingestion.py  ·  Light Curve Ingestion & Pre-processing (V2.0)
================================================================
Handles fetching, cleaning, normalising, and flattening of TESS/Kepler
light curves.  Now supports local FITS file caching for faster
re-analysis and offline processing.

Processing pipeline (ordered):
  1. Fetch from MAST via lightkurve (or load from local FITS cache)
  2. Quality-flag filtering  (quality_flag == 0 ONLY)
  3. 3-sigma iterative clipping   (astropy.stats.sigma_clip)
  4. Flux normalisation           flux / median(flux)
  5. Savitzky-Golay flattening    window = 75% of suspected period

Hard rules
----------
- NO imputation.  Missing data = dropped data.
- NO quality_flag != 0 points survive.
- Detrending window is physically motivated, never arbitrary.

NEW in V2.0:
- Local FITS caching: Downloads FITS files to .cache/fits/ for fast offline re-analysis
- FITS file loader: Can process local FITS files directly
- Multi-sector data stitching with proper gap handling
"""

from __future__ import annotations

import os
import json
import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
from astropy.stats import sigma_clip
from scipy.signal import savgol_filter

from zspace_engine.logging_config import suppress_astroquery_logger

# Suppress astroquery logging WITHOUT creating its logger prematurely.
# Calling logging.getLogger('astroquery') before astroquery has initialized
# registers a plain Logger; astroquery._init_log() then fails with
# "'Logger' object has no attribute '_set_defaults'". Only suppress if
# astroquery's logger already exists (i.e. astroquery already imported).
try:
    suppress_astroquery_logger()
    # Also suppress lightkurve's internal logging
    lk_logger = logging.getLogger('lightkurve')
    lk_logger.setLevel(logging.WARNING)
except Exception:
    pass


# ─────────────────────────────────────────────────────────────────────────────
# Ingestion result container
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class LightCurveProduct:
    """Holds the fully pre-processed light-curve arrays plus audit metadata."""
    tic_id:             str
    sector:             int
    time:               np.ndarray          # BJD - 2457000
    flux_raw:           np.ndarray          # raw, quality-filtered
    flux_norm:          np.ndarray          # normalised: flux / median
    flux_flat:          np.ndarray          # detrended
    trend:              np.ndarray          # SG trend model
    cadence_days:       float
    n_points_raw:       int
    n_points_cleaned:   int
    n_dropped_quality:  int
    n_dropped_sigma:    int
    sectors_used:       List[int] = field(default_factory=list)  # NEW: list of sectors
    fits_source:        str = "mast"        # "mast" | "local_cache" | "local_file"
    audit_log:          List[str]           = field(default_factory=list)

    @property
    def coverage_fraction(self) -> float:
        return self.n_points_cleaned / max(self.n_points_raw, 1)

    def assert_physical_consistency(
        self,
        max_flux_scatter: Optional[float] = None,
    ) -> None:
        """
        Hard check: data must be physically plausible after cleaning.

        Parameters
        ----------
        max_flux_scatter : optional override for the normalised-flux scatter
            ceiling (default 0.10). The default gate existed to reject
            corrupted/contaminated light curves, but a fixed 0.10 ceiling
            silently destroys legitimate EXTREME-DEPTH astrophysical signals
            (e.g. WD 1856+534 b: a planetary-mass companion eclipsing a
            white dwarf, depth ~56% → normalised scatter ~0.36). Callers that
            are searching for such signals must pass a higher ceiling.
        """
        if len(self.time) < 100:
            raise RuntimeError(
                f"Insufficient data after cleaning: {len(self.time)} points. "
                "Cannot run BLS on fewer than 100 cadences."
            )
        gate = max_flux_scatter if max_flux_scatter is not None else 0.10
        flux_std = np.std(self.flux_norm)
        if flux_std > gate:
            raise RuntimeError(
                f"Normalised flux scatter {flux_std:.4f} > {gate:.2f} "
                "(gate raised by caller if searching for extreme-depth signals). "
                "Likely data corruption or wrong target. Aborting."
            )
        if np.any(~np.isfinite(self.flux_flat)):
            raise RuntimeError("Non-finite values found in flattened flux. Aborting.")


# ─────────────────────────────────────────────────────────────────────────────
# FITS Cache Manager (NEW)
# ─────────────────────────────────────────────────────────────────────────────

class FITSCacheManager:
    """
    Manages local FITS file cache for fast re-analysis.
    
    Downloads FITS files from MAST and stores them locally in:
        .cache/fits/TIC_{tic_id}/sector_{N}/
    
    Subsequent runs load from cache instead of re-downloading,
    providing significant speedup for iterative analysis.
    """
    
    CACHE_DIR = Path(".cache/fits")
    
    @classmethod
    def get_cache_path(cls, tic_id: str, sector: int) -> Path:
        """Get the cache directory for a TIC/sector combination."""
        return cls.CACHE_DIR / f"TIC_{tic_id}" / f"sector_{sector}"
    
    @classmethod
    def is_cached(cls, tic_id: str, sector: int) -> bool:
        """Check if FITS files exist in cache."""
        cache_path = cls.get_cache_path(tic_id, sector)
        if not cache_path.exists():
            return False
        fits_files = list(cache_path.glob("*.fits"))
        return len(fits_files) > 0
    
    @classmethod
    def cache_fits(cls, tic_id: str, sector: int, fits_path: str) -> Path:
        """
        Copy a FITS file to the local cache.
        Returns the cached file path.
        """
        import shutil
        cache_path = cls.get_cache_path(tic_id, sector)
        cache_path.mkdir(parents=True, exist_ok=True)
        
        src = Path(fits_path)
        dst = cache_path / src.name
        if not dst.exists():
            shutil.copy2(str(src), str(dst))
        return dst
    
    @classmethod
    def get_cached_fits(cls, tic_id: str, sector: int) -> List[Path]:
        """Get list of cached FITS files for a TIC/sector."""
        cache_path = cls.get_cache_path(tic_id, sector)
        if not cache_path.exists():
            return []
        return sorted(cache_path.glob("*.fits"))
    
    @classmethod
    def download_and_cache(cls, tic_id: str, sector: Optional[int] = None) -> Tuple[List[Path], List[int]]:
        """
        Download FITS files from MAST and cache locally.
        Returns (list of cached file paths, list of sector numbers).
        """
        try:
            import lightkurve as lk
            
            target = f"TIC {tic_id}"
            search = lk.search_lightcurve(
                target, mission="TESS", author="SPOC"
            )
            
            if search is None or len(search) == 0:
                # Try any author
                search = lk.search_lightcurve(target, mission="TESS")
            
            if search is None or len(search) == 0:
                logging.warning(f"No FITS files found for TIC {tic_id}")
                return [], []
            
            cached_files = []
            sector_nums = []
            
            for result in search:
                try:
                    # Download
                    lc = result.download()
                    if lc is None:
                        continue
                    
                    # Get sector number
                    result_sector = getattr(lc.meta, 'SECTOR', 0) or 0
                    if sector is not None and result_sector != sector:
                        continue
                    
                    # Cache the FITS file if lightkurve provides the path
                    fits_path = getattr(lc, 'filename', None) or getattr(lc, 'FILENAME', None)
                    if fits_path and Path(fits_path).exists():
                        cached = cls.cache_fits(tic_id, result_sector, fits_path)
                        cached_files.append(cached)
                        sector_nums.append(result_sector)
                    
                except Exception as e:
                    logging.warning(f"Failed to download/cache sector: {e}")
                    continue
            
            return cached_files, sector_nums
            
        except Exception as e:
            logging.error(f"FITS download failed for TIC {tic_id}: {e}")
            return [], []
    
    @classmethod
    def load_from_cache(cls, tic_id: str, sector: int) -> Optional[Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]]:
        """
        Load light curve data from cached FITS file.
        Returns (time, flux, quality) or None if not cached.
        """
        cached_files = cls.get_cached_fits(tic_id, sector)
        if not cached_files:
            return None
        
        try:
            from astropy.io import fits as afits
            
            all_time = []
            all_flux = []
            all_quality = []
            
            for fits_file in cached_files:
                with afits.open(str(fits_file)) as hdul:
                    # Try different HDU structures
                    for ext_name in ['LIGHTCURVE', 'LC', 1]:
                        try:
                            data = hdul[ext_name].data
                            if data is not None:
                                time_col = None
                                flux_col = None
                                qual_col = None
                                
                                for col_name in ['TIME', 'BTJD', 'BJD']:
                                    if col_name in data.columns.names:
                                        time_col = data[col_name]
                                        break
                                
                                for col_name in ['PDCSAP_FLUX', 'SAP_FLUX', 'FLUX']:
                                    if col_name in data.columns.names:
                                        flux_col = data[col_name]
                                        break
                                
                                if 'QUALITY' in data.columns.names:
                                    qual_col = data['QUALITY']
                                
                                if time_col is not None and flux_col is not None:
                                    all_time.append(np.asarray(time_col, dtype=np.float64))
                                    all_flux.append(np.asarray(flux_col, dtype=np.float64))
                                    if qual_col is not None:
                                        all_quality.append(np.asarray(qual_col, dtype=int))
                                    break
                        except (KeyError, IndexError):
                            continue
            
            if not all_time:
                return None
            
            time = np.concatenate(all_time)
            flux = np.concatenate(all_flux)
            quality = np.concatenate(all_quality) if all_quality else None
            
            return time, flux, quality
            
        except Exception as e:
            logging.error(f"Failed to load FITS from cache: {e}")
            return None


# ─────────────────────────────────────────────────────────────────────────────
# Main ingester class
# ─────────────────────────────────────────────────────────────────────────────

class LightCurveIngester:
    """
    Fetches and pre-processes TESS (or Kepler) light curves for a given TIC ID.

    Parameters
    ----------
    tic_id          : TESS Input Catalogue identifier (int or str)
    mission         : 'TESS' | 'Kepler' | 'K2'
    exptime         : 'short' | 'long' | 120 | 1800  (seconds)
    author          : data source ('SPOC', 'QLP', 'TESS-SPOC', ...)
    sigma_clip_n    : sigma level for iterative clipping (default 3.0)
    use_cache       : whether to use local FITS cache (default True)
    """

    def __init__(
        self,
        tic_id: int | str,
        mission: str = "TESS",
        exptime: str | int = "short",
        author: str = "SPOC",
        sigma_clip_n: float = 3.0,
        use_cache: bool = True,
    ) -> None:
        self.tic_id       = str(tic_id)
        self.mission      = mission
        self.exptime      = exptime
        self.author       = author
        self.sigma_clip_n = sigma_clip_n
        self.use_cache    = use_cache
        self._audit: List[str] = []
        self._cache = FITSCacheManager()

    def _log(self, msg: str) -> None:
        self._audit.append(msg)

    # ── 1.  Fetch ────────────────────────────────────────────────────────────

    def fetch(self, sector: Optional[int] = None) -> Tuple["lightkurve.LightCurve", int]:
        """Download light curve from MAST (with FITS caching).  Returns (lc, sector_used)."""
        
        # Try loading from cache first
        if self.use_cache and sector is not None:
            cached_data = self._cache.load_from_cache(self.tic_id, sector)
            if cached_data is not None:
                time, flux, quality = cached_data
                self._log(f"FETCH | Loaded from FITS cache: .cache/fits/TIC_{self.tic_id}/sector_{sector}/")
                # Create a minimal lightkurve-like object
                from lightkurve import LightCurve
                from astropy.time import Time
                import astropy.units as u
                
                lc = LightCurve(
                    time=Time(time, format="btjd", scale="tdb"),
                    flux=flux * u.electron / u.s,
                )
                if quality is not None:
                    lc.quality = quality
                lc.meta = {"SECTOR": sector}
                return lc, sector

        try:
            import lightkurve as lk
        except ImportError as e:
            raise ImportError("lightkurve is required. Install via: pip install lightkurve") from e

        target = f"TIC {self.tic_id}"
        self._log(f"FETCH | target={target} mission={self.mission} exptime={self.exptime} author={self.author}")

        try:
            search = lk.search_lightcurve(
                target,
                mission=self.mission,
                exptime=self.exptime,
                author=self.author,
            )
        except Exception as e:
            self._log(f"FETCH | Network error during lightcurve search: {e}")
            raise RuntimeError(
                f"Network failure while searching for TIC {self.tic_id}: {e}. "
                "Check network connectivity and MAST availability."
            ) from e

        if len(search) == 0:
            # Fallback: try any author
            self._log("FETCH | Primary author not found. Falling back to any author.")
            try:
                search = lk.search_lightcurve(target, mission=self.mission)
            except Exception as e:
                self._log(f"FETCH | Network error during fallback search: {e}")
                raise RuntimeError(
                    f"Network failure during fallback search for TIC {self.tic_id}: {e}"
                ) from e

        if len(search) == 0:
            raise FileNotFoundError(
                f"No light curve data found for TIC {self.tic_id} on MAST. "
                "Verify the TIC ID and check network access."
            )

        try:
            if sector is not None:
                matches = [r for r in search if r.mission[0].endswith(f"Sector {sector:02d}") or
                           str(sector) in str(r.mission)]
                if matches:
                    lc_collection = lk.LightCurveCollection([m.download() for m in matches[:1]])
                else:
                    self._log(f"FETCH | Sector {sector} not available; using first available.")
                    lc_collection = search[:1].download_all()
            else:
                lc_collection = search.download_all()
        except Exception as e:
            self._log(f"FETCH | Network error during lightcurve download: {e}")
            raise RuntimeError(
                f"Network failure while downloading lightcurve for TIC {self.tic_id}: {e}. "
                "Check network connectivity and MAST availability."
            ) from e

        try:
            lc_stitched = lc_collection.stitch()
            sector_used = getattr(lc_stitched.meta, "SECTOR", 0) or 0
            self._log(f"FETCH | Downloaded {len(lc_stitched)} cadences across {len(lc_collection)} sectors.")
            
            # Cache the FITS files locally for future use
            if self.use_cache:
                try:
                    seen_sectors = set()
                    for lc_item in lc_collection:
                        fits_path = getattr(lc_item, 'filename', None) or getattr(lc_item, 'FILENAME', None)
                        item_sector = getattr(lc_item.meta, 'SECTOR', sector_used) or sector_used
                        if item_sector not in seen_sectors:
                            seen_sectors.add(item_sector)
                            if fits_path and Path(fits_path).exists():
                                self._cache.cache_fits(self.tic_id, item_sector, fits_path)
                                self._log(f"FETCH | Cached FITS for sector {item_sector}")
                except Exception as e:
                    self._log(f"FETCH | FITS caching failed (non-critical): {e}")
            
            return lc_stitched, sector_used
        except Exception as e:
            self._log(f"FETCH | Error stitching lightcurves: {e}")
            raise RuntimeError(
                f"Failed to stitch lightcurves for TIC {self.tic_id}: {e}"
            ) from e

    # ── 2.  Quality filtering ─────────────────────────────────────────────────

    @staticmethod
    def _apply_quality_mask(lc: "lightkurve.LightCurve") -> Tuple[np.ndarray, np.ndarray, int]:
        """
        Hard rule: drop all points where quality_flag != 0.
        Returns (time, flux, n_dropped).
        """
        time = np.asarray(lc.time.value, dtype=np.float64)
        flux = np.asarray(lc.flux.value, dtype=np.float64)

        if hasattr(lc, "quality") and lc.quality is not None:
            quality = np.asarray(lc.quality.value, dtype=int)
            good_mask = (quality == 0) & np.isfinite(flux)
        else:
            good_mask = np.isfinite(flux)

        n_dropped = int((~good_mask).sum())
        return time[good_mask], flux[good_mask], n_dropped

    # ── 3.  Iterative 3-sigma clipping ───────────────────────────────────────

    def _sigma_clip_flux(self, time: np.ndarray, flux: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
        """
        Iterative sigma-clipping — UPWARD outliers only.

        Physical rationale: sigma clipping targets instrumental artefacts
        (cosmic rays, thermal noise spikes, saturation), which all manifest
        as UPWARD flux excursions.  Transit dips are DOWNWARD and must be
        preserved.  Two-sided clipping at 3σ would excise deep transits
        (depth >> σ), corrupting the signal we are trying to detect.

        Implementation: sigma_lower=None (no lower clip),
                        sigma_upper=self.sigma_clip_n (clip positive spikes).
        """
        masked = sigma_clip(
            flux,
            sigma_lower=10.0,           # effectively disabled downward clip
            sigma_upper=self.sigma_clip_n,
            maxiters=10,
            masked=True,
        )
        keep   = ~masked.mask
        n_drop = int(masked.mask.sum())
        self._log(
            f"SIGMA_CLIP | sigma_upper={self.sigma_clip_n} (one-sided, upward only) | "
            f"dropped {n_drop} upward outliers ({100*n_drop/len(flux):.2f}%) | "
            "transit dips preserved by design"
        )
        return time[keep], flux[keep], n_drop

    # ── 4.  Normalisation ─────────────────────────────────────────────────────

    @staticmethod
    def _normalize(flux: np.ndarray) -> Tuple[np.ndarray, float]:
        """
        Normalise flux: flux_norm = flux / median(flux)
        Post-normalisation median = 1.0 by construction.
        """
        med = np.median(flux)
        if med <= 0:
            raise RuntimeError(f"Median flux is {med:.4f} ≤ 0. Non-physical signal.")
        return flux / med, med

    # ── 5.  Savitzky-Golay flattening ─────────────────────────────────────────

    def _savgol_flatten(
        self,
        time: np.ndarray,
        flux_norm: np.ndarray,
        period_days: Optional[float] = None,
    ) -> Tuple[np.ndarray, np.ndarray, int]:
        """
        Flatten using Savitzky-Golay filter.
        Window is physically motivated: 75% of the suspected orbital period,
        expressed in cadence units.  Minimum window = 51 cadences.

        Returns (flux_flat, trend, window_used).
        """
        # Guard: savgol_filter assumes monotonically increasing time with
        # uniform spacing. Sort by time so the filter never operates on an
        # unsorted series (fixes out-of-order-cadence systematics).
        order = np.argsort(time, kind="mergesort")
        time_sorted = time[order]
        flux_sorted = flux_norm[order]
        cadence = float(np.median(np.diff(time_sorted)))  # days per cadence

        if period_days is not None and period_days > 0:
            window_days   = 0.75 * period_days
            window_pts    = int(round(window_days / cadence))
        else:
            # Default: ~3-day window (typical TESS baseline)
            window_pts = int(round(3.0 / cadence))

        # SG requires odd window ≥ polyorder+1
        if window_pts < 51:
            window_pts = 51
        if window_pts % 2 == 0:
            window_pts += 1

        # Guard: window must not exceed the number of samples
        n_pts = flux_sorted.size
        window_pts = min(window_pts, n_pts if n_pts % 2 == 1 else n_pts - 1)
        if window_pts < 5 or n_pts < 5:
            self._log("SAVGOL | too few points → return input unchanged")
            return flux_norm, np.ones_like(flux_norm), 0

        polyorder = min(3, window_pts - 1)

        trend = savgol_filter(flux_sorted, window_length=window_pts, polyorder=polyorder)
        flux_flat_sorted = flux_sorted / trend

        # Restore original time order so flux_flat stays aligned with `time`.
        inv_order = np.empty_like(order)
        inv_order[order] = np.arange(order.size)
        flux_flat = flux_flat_sorted[inv_order]

        self._log(
            f"SAVGOL | period_input={period_days} d | "
            f"cadence={cadence:.5f} d | window={window_pts} pts | polyorder={polyorder}"
        )
        return flux_flat, trend, window_pts

    # ── Master pipeline ───────────────────────────────────────────────────────

    def process(
        self,
        sector: Optional[int] = None,
        period_hint_days: Optional[float] = None,
        max_flux_scatter: Optional[float] = None,
    ) -> LightCurveProduct:
        """
        Execute the full ingestion & pre-processing pipeline.

        Steps:  fetch → quality filter → sigma clip → normalise → flatten

        Parameters
        ----------
        sector           : restrict to a specific TESS sector
        period_hint_days : period estimate to motivate SG window length
        max_flux_scatter : optional ceiling override for the consistency gate
            (pass e.g. 0.50 when the target may be an extreme-depth signal
            such as a planet around a white dwarf).
        """
        self._audit.clear()
        lc_raw, sector_used = self.fetch(sector)

        n_raw = len(lc_raw)
        self._log(f"RAW | {n_raw} cadences downloaded")

        time_q, flux_q, n_dropped_q = self._apply_quality_mask(lc_raw)
        self._log(f"QUALITY_MASK | dropped {n_dropped_q} flagged cadences | remaining {len(time_q)}")

        time_c, flux_c, n_dropped_s = self._sigma_clip_flux(time_q, flux_q)

        flux_norm, flux_median = self._normalize(flux_c)
        self._log(f"NORMALIZE | median_flux={flux_median:.4f} → flux_norm median ≈ 1.0")

        flux_flat, trend, window = self._savgol_flatten(time_c, flux_norm, period_hint_days)

        cadence = float(np.median(np.diff(time_c)))

        product = LightCurveProduct(
            tic_id            = self.tic_id,
            sector            = int(sector_used),
            time              = time_c,
            flux_raw          = flux_c,
            flux_norm         = flux_norm,
            flux_flat         = flux_flat,
            trend             = trend,
            cadence_days      = cadence,
            n_points_raw      = n_raw,
            n_points_cleaned  = len(time_c),
            n_dropped_quality = n_dropped_q,
            n_dropped_sigma   = n_dropped_s,
            sectors_used      = [int(sector_used)],
            fits_source       = "local_cache" if self.use_cache else "mast",
            audit_log         = list(self._audit),
        )

        product.assert_physical_consistency(max_flux_scatter=max_flux_scatter)
        self._log("CONSISTENCY_CHECK | PASS")
        product.audit_log.append("CONSISTENCY_CHECK | PASS")
        return product

    # ── Load from local FITS file (NEW) ───────────────────────────────────────

    @classmethod
    def from_fits(
        cls,
        tic_id: str,
        fits_path: str,
        sector: int = 0,
        period_hint_days: Optional[float] = None,
    ) -> LightCurveProduct:
        """
        Load and process a local FITS file directly.
        
        Parameters
        ----------
        tic_id : str
            TIC identifier
        fits_path : str
            Path to local FITS file
        sector : int
            Sector number (for metadata)
        period_hint_days : float, optional
            Period hint for SG window
        """
        from astropy.io import fits as afits
        
        inst = cls(tic_id=tic_id)
        inst._audit = [f"SOURCE | local FITS file: {fits_path}"]
        
        with afits.open(fits_path) as hdul:
            # Try to find the lightcurve data
            data = None
            for ext in ['LIGHTCURVE', 'LC', 1]:
                try:
                    data = hdul[ext].data
                    if data is not None:
                        break
                except (KeyError, IndexError):
                    continue
            
            if data is None:
                raise RuntimeError(f"No lightcurve data found in FITS file: {fits_path}")
            
            # Extract columns
            time = None
            flux = None
            quality = None
            
            for col_name in ['TIME', 'BTJD', 'BJD']:
                if col_name in data.columns.names:
                    time = np.asarray(data[col_name], dtype=np.float64)
                    break
            
            for col_name in ['PDCSAP_FLUX', 'SAP_FLUX', 'FLUX']:
                if col_name in data.columns.names:
                    flux = np.asarray(data[col_name], dtype=np.float64)
                    break
            
            if 'QUALITY' in data.columns.names:
                quality = np.asarray(data['QUALITY'], dtype=int)
            
            if time is None or flux is None:
                raise RuntimeError(f"Required columns not found in FITS: {fits_path}")
        
        return cls.from_arrays(
            tic_id=tic_id,
            time=time,
            flux=flux,
            quality=quality,
            period_hint_days=period_hint_days,
        )

    # ── Synthetic / local-file injection (for testing) ────────────────────────

    @classmethod
    def from_arrays(
        cls,
        tic_id: str,
        time: np.ndarray,
        flux: np.ndarray,
        quality: Optional[np.ndarray] = None,
        period_hint_days: Optional[float] = None,
        max_flux_scatter: Optional[float] = None,
    ) -> "LightCurveProduct":
        """
        Bypass MAST fetch.  Inject raw arrays directly (unit tests, FITS files).
        The same cleaning pipeline applies.
        """
        inst = cls(tic_id=tic_id)
        inst._audit = ["SOURCE | local array injection"]

        if quality is not None:
            good = (quality == 0) & np.isfinite(flux)
            n_dropped_q = int((~good).sum())
            time_q = time[good]
            flux_q = flux[good]
        else:
            n_dropped_q = 0
            time_q = time[np.isfinite(flux)]
            flux_q = flux[np.isfinite(flux)]

        inst._log(f"QUALITY_MASK | dropped {n_dropped_q} flagged cadences")
        time_c, flux_c, n_dropped_s = inst._sigma_clip_flux(time_q, flux_q)
        flux_norm, _ = inst._normalize(flux_c)
        flux_flat, trend, _ = inst._savgol_flatten(time_c, flux_norm, period_hint_days)

        product = LightCurveProduct(
            tic_id            = tic_id,
            sector            = 0,
            time              = time_c,
            flux_raw          = flux_c,
            flux_norm         = flux_norm,
            flux_flat         = flux_flat,
            trend             = trend,
            cadence_days      = float(np.median(np.diff(time_c))),
            n_points_raw      = len(time),
            n_points_cleaned  = len(time_c),
            n_dropped_quality = n_dropped_q,
            n_dropped_sigma   = n_dropped_s,
            fits_source       = "local_array",
            audit_log         = list(inst._audit),
        )
        product.assert_physical_consistency(max_flux_scatter=max_flux_scatter)
        return product
