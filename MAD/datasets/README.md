Datasets are not bundled in this release, except for a tiny `sample.edgelist` for smoke tests.

Expected layout in this directory:
- *.edgelist or *.obj graph files
- *.comms community files will be produced after the first usage

Recommended sources (SNAP datasets used in the paper):
- https://snap.stanford.edu/data/

Place each dataset using the same naming convention referenced by configs in `experiment_configs/`.

Sample usage:
```bash
PYTHONPATH="$PWD/MAD" python MAD/experimentsMAD.py --dataset sample --detector leiden --budget 5
```
