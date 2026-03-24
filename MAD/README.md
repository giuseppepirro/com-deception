# MAD: Multilevel Community Deception (Release)

This directory contains a clean release snapshot of the codebase (without internal logs, configs, or private modules). It supports:
- CSD (Community Structure Deception)
- SCD (Single Community Deception)
- IND (Individual Node Deception)

## Structure
- `MAD/`: core algorithms and utilities
- `experiment_configs/`: JSON configs for batch runs
- `scripts/`: helper scripts (including IND report)
- `datasets/`: placeholder directory 

## Setup
Python 3.8+ recommended.

Install core dependencies:
```bash
pip install -r requirements.txt
```

## Datasets
Datasets are not bundled. Download directed graphs from SNAP and place them under `datasets/`.
See `datasets/README.md` for expected file naming and layout. A tiny `sample.edgelist` is included for smoke tests.

SNAP: https://snap.stanford.edu/data/

## Run Experiments
Batch run (CSD/SCD/IND) using a JSON config:
```bash
python run_experiments.py experiment_configs/experiments.json
```
The default config targets the bundled `sample` dataset.

Single MAD run (example):
```bash
PYTHONPATH="$PWD/MAD" python MAD/experimentsMAD.py --dataset anybeat --detector leiden --budget 5
```

## Results
Results are written under `results/<method>/` (created at runtime). Each run stores its parameters and metrics.

## Notes
- GEMSEC source is included under `third_party/GEMSEC`, but it requires a legacy TensorFlow setup to run.
