# Tests Directory

This directory contains test and validation scripts for the Axiom-ZSpace exoplanet discovery engine.

## Contents

### run_validator.py
Validation test script that demonstrates the AxiomValidator module with both synthetic and real TIC targets.

**Usage:**
```bash
# Run synthetic test (no internet required)
python tests/run_validator.py --synthetic

# Test with real TIC ID (requires internet)
python tests/run_validator.py --tic 260128333

# Simulate known planet match
python tests/run_validator.py --force-known
```

### test_synthetic.py
Reusable synthetic test data generation module. Provides functions for:
- Generating synthetic planet parameters
- Generating known planet parameters
- Printing formatted discovery summaries

**Usage:**
```python
from test_synthetic import generate_synthetic_planet_params

params = generate_synthetic_planet_params()
# Use params for testing...
```

## Production Entry Point

The production entry point for the pipeline is `run_pipeline.py` in the root directory, not in this tests/ directory.

```bash
# Production usage
python run_pipeline.py --tic 260128333
python run_pipeline.py --sector 42
```
