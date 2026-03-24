import os
import timeit
from math import log
from datetime import datetime
from pathlib import Path

import igraph as ig
import numpy as np
from cdlib import evaluation, viz  # viz unused but kept to match original imports

import MAD
from Utilities_DIRECTED import *
from ndre_counter import count_pre_security_index


def _configure_cache_dirs():
    """Ensure matplotlib/fontconfig caches point to writable local paths to avoid warnings."""
    base = Path(__file__).resolve().parent / ".cache"
    mpl_dir = base / "matplotlib"
    xdg_cache = base
    for d in (mpl_dir, xdg_cache):
        d.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))


_configure_cache_dirs()


def compute_security_index(graph: ig.Graph, communities):
    """Compute Hp/H security index given a graph and community partition."""
    try:
        name_to_idx = {n: i for i, n in enumerate(graph.vs["name"])} if "name" in graph.vs.attributes() else None
        # normalize communities to vertex indices
        norm_coms = []
        for com in communities:
            norm = []
            for node in com:
                if isinstance(node, int):
                    norm.append(node)
                elif name_to_idx is not None and node in name_to_idx:
                    norm.append(name_to_idx[node])
            norm_coms.append(norm)

        m = graph.ecount()
        if m <= 0:
            return None
        out_degs = graph.strength(graph.vs, mode='out', weights=None)
        # H(G): structural entropy of out-degree distribution
        H = 0.0
        for d in out_degs:
            if d > 0:
                p = d / m
                H -= p * log(p, 2)
        if H <= 0:
            return None
        # H(G) - H^P(G): sum over communities using intra out-edges
        gap = 0.0
        for com in norm_coms:
            if not com:
                continue
            vol = sum(out_degs[i] for i in com)
            if vol <= 0:
                continue
            inner = graph.induced_subgraph(com, implementation="auto").ecount()
            if inner <= 0:
                continue
            gap += (inner / m) * log(m / vol, 2)
        return gap / H
    except Exception:
        return None


class expMADDirected:
    def __init__(self, dataset_path, datasets, detection_algos, budgets, is_budget_percentage, rounds, level, node_criterion, allowed_del_ratio=0.5):
        self.dataset_path = dataset_path
        self.datasets = datasets
        self.detection_algos = detection_algos
        self.budget_updates = budgets
        self.del_ratio = allowed_del_ratio
        self.is_budget_percentage = is_budget_percentage
        self.rounds = rounds
        self.runner = Utils_DIRECTED(dataset_path=self.dataset_path)
        self.deception_level = level
        self.target_criterion = node_criterion
        self.run_type = {"structure": "csd", "community": "scd", "individual": "ind"}.get(self.deception_level, "csd")
        suffix = f"{self.run_type}_directed_"
        if self.run_type == "ind":
            safe_criterion = str(self.target_criterion).replace(" ", "_").replace("+", "plus")
            suffix = f"{self.run_type}_{safe_criterion}_directed_"
        results_dir = Path("results/mad")
        results_dir.mkdir(parents=True, exist_ok=True)
        exp_file_name = str(results_dir / f"{suffix}{datetime.today().strftime('%Y-%m-%d-%H:%M:%S')}.experiments")
        self.file_experiments = open(exp_file_name, "w+")
        self.current_algo = None
        self.experiment_configuration = False

    def runExperiment(self, network, detection_algo, budget, is_budget_percentage):
        current_experiment_configuration = network + "_" + detection_algo
        self.runner.read_directed_graph(network)
        if current_experiment_configuration != self.experiment_configuration:
            recompute_communities = True
            self.runner.initializeCommunityData(detection_algo=detection_algo, dataset=network, recompute_communities=recompute_communities)
            self.experiment_configuration = current_experiment_configuration

        coms_before = self.runner.communities_object
        num_coms_before = len(coms_before.communities)
        deception_before, member_for_community = self.runner.getDeceptionScore(coms_before.communities)

        ebd2n = MAD.MAD(self.runner)
        graph_before = self.runner.graph.copy()
        step_metrics = []

        def capture_metrics(g, coms):
            if len(coms) == 0:
                return
            try:
                ref = self.runner.communities_object
                nmi = evaluation.normalized_mutual_information(ref, self.runner.communities_after_object)[0]
                jacc, rec = self.runner.getJaccardAndRecallScores(ref.communities, self.runner.communities_after_object.communities)
                sec = compute_security_index(g, coms)
                step_metrics.append((nmi, jacc, rec, sec))
            except Exception:
                pass

        ebd2n.metrics_callback = capture_metrics

        start = timeit.default_timer()
        ebd2n.initializeDataStructuresDeceptionALLCommunities()
        if self.deception_level == "structure":
            updated_graph = ebd2n.runEBD2NALL(budget=budget, budget_percentage=is_budget_percentage, allowed_del=self.del_ratio, lev="structure")
        elif self.deception_level == "community":
            ebd2n.initializeDataStructures()
            updated_graph = ebd2n.runEBD2NALL(budget=budget, budget_percentage=is_budget_percentage, allowed_del=self.del_ratio, lev="community")
        elif self.deception_level == "individual":
            ebd2n.initializeDataStructures()
            ebd2n.initializeIndividualDeceptionDataStructures(criterion=self.target_criterion)
            self.runner.getPreDeceptionTargetIndividualMainAttributes(ebd2n.target_node_id)
            updated_graph = ebd2n.runEBD2NALL(budget=budget, budget_percentage=is_budget_percentage, allowed_del=self.del_ratio, lev="individual", criterion=self.target_criterion)
        stop = timeit.default_timer()
        time = stop - start

        coms_after = self.runner.computeCommunitiesAfterUpdateCDLIB(detection_algo, updated_graph)
        coms_after = self.runner.communities_after_object
        if self.deception_level == "individual":
            self.runner.getPostDeceptionTargetIndividualMainAttributes(updated_graph, coms_after.communities, ebd2n.target_node_id)

        num_coms_after = len(coms_after.communities)
        deception_after, member_for_community = self.runner.getDeceptionScore(coms_after.communities)
        member_for_community = np.count_nonzero(np.array(member_for_community) > 0)

        nmi_before, jaccard_before, recall_before = "n/a", "n/a", "n/a"
        nmi_after, jaccard_after, recall_after = "n/a", "n/a", "n/a"
        if len(coms_before.communities) > 0:
            nmi_before, jaccard_before, recall_before = 1.0, 1.0, 1.0  # self-comparison
        if len(coms_before.communities) > 0 and len(coms_after.communities) > 0:
            nmi_after = evaluation.normalized_mutual_information(coms_before, coms_after)[0]
            jaccard_after, recall_after = self.runner.getJaccardAndRecallScores(coms_before.communities, coms_after.communities)

        sec_before_val = compute_security_index(graph_before, coms_before.communities)
        sec_after_val = compute_security_index(updated_graph, coms_after.communities)
        sec_before = sec_before_val if sec_before_val is not None else "n/a"
        sec_after = sec_after_val if sec_after_val is not None else "n/a"

        if self.deception_level == "structure":
            line = [
                self.run_type,
                network,
                detection_algo,
                str(budget),
                f"{time}",
                str(num_coms_before),
                str(num_coms_after),
                f"{nmi_before}",
                f"{jaccard_before}",
                f"{recall_before}",
                f"{nmi_after}",
                f"{jaccard_after}",
                f"{recall_after}",
                f"{sec_before}",
                f"{sec_after}",
            ]
        elif self.deception_level == "community":
            line = [
                self.run_type,
                network,
                detection_algo,
                str(budget),
                f"{time}",
                str(deception_before),
                str(deception_after),
                f"{nmi_before}",
                f"{jaccard_before}",
                f"{recall_before}",
                f"{nmi_after}",
                f"{jaccard_after}",
                f"{recall_after}",
                f"{sec_before}",
                f"{sec_after}",
                str(self.runner.target_community_id),
                str(self.runner.community_size),
                str(member_for_community),
            ]
        else:
            label_free_score = getattr(ebd2n, "label_free_score", 0.0)
            max_deg = getattr(ebd2n, "original_max_outdegree", 1.0)
            # Max ego-score per edit occurs at sim=1 and max out-degree
            per_edit_cap = log(1.0 + max_deg, 2)
            max_possible = max(per_edit_cap * max(1, getattr(ebd2n, "beta", budget)), 1e-9)
            label_free_norm = max(0.0, min(1.0, label_free_score / max_possible))
            line = [
                self.run_type,
                network,
                detection_algo,
                str(budget),
                f"{time}",
                str(ebd2n.target_node_id),
                self.target_criterion,
                str(label_free_score),
                f"{label_free_norm}",
            ]

        out_str = "\t".join(line)
        self.file_experiments.write("\n" + out_str)
        print(out_str)

        if self.deception_level in ("structure", "community") and step_metrics:
            steps_path = f"results/mad/{self.run_type}_{network}_{detection_algo}_steps.txt"
            with open(steps_path, "w") as sf:
                sf.write("Step\tNMI\tJaccard\tRecall\tSecurity\n")
                for idx, (snmi, sjacc, srec, ssec) in enumerate(step_metrics, 1):
                    sf.write(f"{idx}\t{snmi}\t{sjacc}\t{srec}\t{ssec}\n")

    def runAllExperiments(self):
        if self.deception_level == "structure":
            header = [
                "RunType",
                "Network",
                "DetectionAlgorithm",
                "BudgetUpdates",
                "Time(s)",
                "#CommunitiesBefore",
                "#CommunitiesAfter",
                "NMI_Before",
                "Jaccard_Before",
                "Recall_Before",
                "NMI_After",
                "Jaccard_After",
                "Recall_After",
                "SecurityBefore",
                "SecurityAfter",
            ]
        elif self.deception_level == "community":
            header = [
                "RunType",
                "Network",
                "DetectionAlgorithm",
                "BudgetUpdates",
                "Time(s)",
                "InitialHidingScore",
                "FinalHidingScore",
                "NMI_Before",
                "Jaccard_Before",
                "Recall_Before",
                "NMI_After",
                "Jaccard_After",
                "Recall_After",
                "SecurityBefore",
                "SecurityAfter",
                "TargetComID",
                "SizeTargetCom",
                "CommunitySplit",
            ]
        else:
            header = [
                "RunType",
                "Network",
                "DetectionAlgorithm",
                "BudgetUpdates",
                "Time(s)",
                "TargetNode",
                "TargetCriterion",
                "LabelFreeScore",
                "LabelFreeScoreNorm01",
            ]

        header_line = "\t".join(header)
        self.file_experiments.write(header_line)
        print(header_line)

        for network in self.datasets:
            for com_algo in self.detection_algos:
                for i in range(0, self.rounds):
                    for budget in self.budget_updates:
                        self.runExperiment(network, com_algo, budget, self.is_budget_percentage)
        self.file_experiments.close()


def main():
    dataset_path = "datasets/"
    datasets = ["aca"]
    com_algos = ["gem"]
    level = "community"
    target_node_criterion = "degree"
    budget_updates = [100]
    is_budget_percentage = False
    allowed_del_ratio = 1.0
    num_rounds = 1

    experiments = expMADDirected(
        dataset_path,
        datasets,
        com_algos,
        budget_updates,
        is_budget_percentage,
        num_rounds,
        level,
        target_node_criterion,
        allowed_del_ratio,
    )
    experiments.runAllExperiments()


if __name__ == "__main__":
    main()
