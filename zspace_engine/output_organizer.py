"""
Output organization system for Axiom-ZSpace discovery pipeline.

This module manages the structured directory hierarchy for discovery outputs,
organizing results by sector and discovery status (new discoveries, known planets,
false positives).

Structure:
    axiom_output/
    ├── sector_1/
    │   ├── discoveries/
    │   │   └── Discovery_ZS-T-12345678-01_20240115_184500.json
    │   ├── known/
    │   │   └── exist_planet_ZS-T-87654321-01_20240115_184500.json
    │   ├── rejected/
    │   │   └── false_positive_ZS-T-11111111-01_20240115_184500.json
    │   └── summary.json
    └── sector_2/
        └── ...

Requirements: 5.1, 5.2
"""

from pathlib import Path
from datetime import datetime
from typing import Optional


class OutputOrganizer:
    """
    Manages structured output directory hierarchy for discovery pipeline.
    
    Organizes outputs by sector and discovery status, ensuring proper directory
    structure and preventing file overwrites through timestamping.
    """
    
    def __init__(self, base_dir: str = "axiom_output"):
        """
        Initialize OutputOrganizer with base directory.
        
        Args:
            base_dir: Root directory for all outputs (default: "axiom_output")
        
        Creates the base directory if it doesn't exist.
        
        Requirements: 5.1, 5.2
        """
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
    
    def get_output_path(
        self,
        sector: int,
        status: str,
        zspace_id: str
    ) -> Path:
        """
        Returns appropriate output path based on discovery status.
        
        Routes Discovery Cards to appropriate subdirectories:
        - NEW_DISCOVERY → axiom_output/sector_N/discoveries/
        - KNOWN → axiom_output/sector_N/known/
        - FALSE_POSITIVE → axiom_output/sector_N/rejected/
        
        Adds timestamps to prevent overwrites and creates subdirectories as needed.
        
        Args:
            sector: TESS sector number
            status: Discovery status ("NEW_DISCOVERY", "KNOWN", "FALSE_POSITIVE")
            zspace_id: ZSpace identifier (e.g., "ZS-T-12345678-01")
        
        Returns:
            Path object for the output file with timestamp
        
        Requirements: 5.3, 5.4, 5.5, 5.8
        """
        sector_dir = self.base_dir / f"sector_{sector}"
        
        # Route to appropriate subdirectory based on status
        if status in ("NEW_DISCOVERY", "OFFLINE_NEW_DISCOVERY"):
            subdir = sector_dir / "discoveries"
            filename = f"Discovery_{zspace_id}.json"
        elif status == "KNOWN":
            subdir = sector_dir / "known"
            filename = f"exist_planet_{zspace_id}.json"
        else:  # FALSE_POSITIVE
            subdir = sector_dir / "rejected"
            filename = f"false_positive_{zspace_id}.json"
        
        # Create subdirectories if they don't exist
        subdir.mkdir(parents=True, exist_ok=True)
        
        # Add timestamp to prevent overwrites
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        stem = filename.rsplit('.', 1)[0]
        return subdir / f"{stem}_{timestamp}.json"
