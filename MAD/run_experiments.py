import json
import sys
from datetime import datetime
from pathlib import Path
import importlib.util


ROOT = Path(__file__).resolve().parent
RESULTS_LOG = ROOT / "results" / "driver_runs.log"


def load_config(path: Path):
    with open(path, "r") as f:
        return json.load(f)


def log_run(entry: str):
    RESULTS_LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_LOG, "a") as f:
        f.write(entry + "\n")
    print(entry)


def save_params(path: Path, params: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(params, f, indent=2)


def run_mad(run_cfg, dataset_root):
    print(f"Running MAD with config: {run_cfg}")
    mad_root = ROOT / "MAD"
    sys.path.insert(0, str(mad_root))
    exp_path = mad_root / "experimentsMAD.py"
    if not exp_path.exists():
        raise FileNotFoundError(f"Missing {exp_path}")
    spec = importlib.util.spec_from_file_location("experimentsMAD", exp_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader, "Failed to load experimentsMAD module spec"
    spec.loader.exec_module(module)
    expMADDirected = module.expMADDirected

    dataset = run_cfg.get("dataset")
    detector = run_cfg.get("detector", "leiden")
    budget = run_cfg.get("budget", 100)
    budget_percentage = bool(run_cfg.get("budget_percentage", False))
    mode = run_cfg.get("mode", "structure")
    target_criterion = run_cfg.get("target_criterion", "degree")
    rounds = int(run_cfg.get("rounds", 1))

    experiments = expMADDirected(
        dataset_root,
        [dataset],
        [detector],
        [budget],
        budget_percentage,
        rounds,
        mode,
        target_criterion,
        allowed_del_ratio=1.0,
    )
    experiments.runAllExperiments()
    params_path = ROOT / "results" / "mad" / f"params_mad_{dataset}_{detector}_{datetime.now().strftime('%Y%m%d-%H%M%S')}.json"
    save_params(params_path, run_cfg)
    log_run(
        f"{datetime.now().isoformat()} method=mad mode={mode} dataset={dataset} detector={detector} "
        f"budget={budget} budget_pct={budget_percentage} target_criterion={target_criterion}"
    )


def main():
    cfg_path = Path(sys.argv[1]) if len(sys.argv) > 1 else (ROOT / "experiment_configs" / "experiments.json")
    cfg = load_config(cfg_path)
    dataset_root = str((ROOT / cfg.get("dataset_root", "datasets")).resolve())

    for run_cfg in cfg.get("runs", []):
        method = run_cfg.get("method", "").lower()
        if method == "mad":
            run_mad(run_cfg, dataset_root + "/")
        else:
            log_run(f"{datetime.now().isoformat()} method={method} status=skipped reason=unsupported")


if __name__ == "__main__":
    main()
