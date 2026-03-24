import argparse
import math
import sys
from pathlib import Path
import signal
from MAD import Utilities_DIRECTED

import igraph as ig
from cdlib import algorithms

from MAD.Utilities_DIRECTED import Utils_DIRECTED

ROOT = Path(__file__).resolve().parents[1]
MAD_PATH = ROOT / "MAD"
if str(MAD_PATH) not in sys.path:
    sys.path.append(str(MAD_PATH))

import MAD as MAD_mod  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Run IND deception and report DR/LE across detectors")
    p.add_argument("--dataset", default="anybeat", help="Dataset name (without extension)")
    p.add_argument("--base-detector", default="leiden", help="Detector used to initialize MAD")
    p.add_argument("--budget", type=int, default=30, help="Edit budget")
    p.add_argument("--budget-percentage", action="store_true", help="Treat budget as percentage")
    p.add_argument("--criterion", default="degree", help="Target node criterion")
    p.add_argument("--target-node", default=None, help="Override target node id or name")
    p.add_argument(
        "--targets",
        default="T1,T2,T3",
        help="Comma-separated target types (T1,T2,T3) for full IND report",
    )
    p.add_argument(
        "--detectors",
        default="leiden,infomap,dm,surprise,ds,gem,walk,tc",
        help="Comma-separated detector list for DR/LE",
    )
    p.add_argument("--timeout", type=int, default=30, help="Seconds per detector computation")
    return p.parse_args()


def compute_partition(g: ig.Graph, name: str):
    name = name.lower()
    if name == "leiden":
        return algorithms.leiden(g)
    if name == "infomap":
        return algorithms.infomap(g)
    if name == "dm":
        return algorithms.rb_pots(g)
    if name == "surprise":
        return algorithms.surprise_communities(g)
    if name == "ds":
        return algorithms.gdmp2(g)
    if name == "gem" or name == "gemsec":
        if hasattr(algorithms, "gemsec"):
            return algorithms.gemsec(g)
        try:
            from karateclub import GEMSEC
            from cdlib.classes import NodeClustering
        except Exception as exc:
            raise RuntimeError(
                "GEMSEC not available. Install karateclub or use a different detector."
            ) from exc
        import networkx as nx
        H = nx.Graph()
        H.add_nodes_from(range(g.vcount()))
        H.add_edges_from(g.get_edgelist())
        n = max(H.number_of_nodes(), 1)
        clusters = max(2, int(n ** 0.5))
        model = GEMSEC(dimensions=64, clusters=clusters)
        model.fit(H)
        memberships = model.get_memberships()
        comm_map = {}
        for node, cid in memberships.items():
            comm_map.setdefault(cid, []).append(node)
        communities = list(comm_map.values())
        return NodeClustering(communities, graph=g, method_name="gemsec-karateclub")
    if name == "walk":
        return algorithms.walktrap(g)
    if name == "tc":
        return algorithms.threshold_clustering(g)
    raise ValueError(f"Unknown detector: {name}")


def compute_partition_timeout(g: ig.Graph, name: str, timeout: int):
    if timeout <= 0:
        return compute_partition(g, name)

    def _handler(signum, frame):
        raise TimeoutError(f"Timeout after {timeout}s")

    old_handler = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(timeout)
    try:
        return compute_partition(g, name)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)


def find_community_id(communities, node_name):
    for idx, com in enumerate(communities):
        if node_name in com:
            return idx
    return None


def ensure_weights(g: ig.Graph):
    if "weight" not in g.es.attributes():
        g.es["weight"] = [1.0] * g.ecount()


def resolve_target_override(runner, mad, target_node):
    # Resolve override target (int id or vertex name)
    target_id = None
    try:
        target_id = int(target_node)
    except Exception:
        pass
    if target_id is None:
        try:
            target_id = runner.graph.vs.find(name=target_node).index
        except Exception as exc:
            raise SystemExit(f"Unknown target node: {target_node}") from exc
    mad.target_node_id = target_id
    # rebuild ego caches for the new target
    ego_nodes = {target_id}
    first_hop = set(runner.graph.neighbors(target_id, mode="out")) | set(
        runner.graph.neighbors(target_id, mode="in")
    )
    ego_nodes |= first_hop
    second_hop = set()
    for n in first_hop:
        second_hop |= set(runner.graph.neighbors(n, mode="out")) | set(
            runner.graph.neighbors(n, mode="in")
        )
    ego_nodes |= second_hop
    mad.ego_nodes = ego_nodes
    mad.ego_neighbor_cache = {
        n: set(runner.graph.neighbors(n, mode="out")).intersection(ego_nodes) for n in ego_nodes
    }
    mad.ego_outdegree = {n: len(neigh) for n, neigh in mad.ego_neighbor_cache.items()}
    return target_id


def main():
    args = parse_args()

    detectors = [d.strip() for d in args.detectors.split(",") if d.strip()]
    if args.target_node is not None:
        targets = [args.criterion]
    else:
        targets = [t.strip() for t in args.targets.split(",") if t.strip()]

    per_target_results = []
    global_failures = {}
    dr_agg = {d: [] for d in detectors}
    le_list = []

    for criterion in targets:
        runner = Utils_DIRECTED(dataset_path="datasets/")
        runner.read_directed_graph(args.dataset)
        runner.initializeCommunityData(detection_algo=args.base_detector, dataset=args.dataset, recompute_communities=True)

        mad = MAD_mod.MAD(runner)
        mad.initializeDataStructures()
        mad.initializeIndividualDeceptionDataStructures(criterion=criterion)

        if args.target_node is not None:
            resolve_target_override(runner, mad, args.target_node)

        target_id = mad.target_node_id
        target_name = runner.graph.vs[target_id]["name"] if "name" in runner.graph.vs.attributes() else target_id

        pre_graph = runner.graph.copy()

        # Run IND deception
        updated_graph = mad.runEBD2NALL(
            budget=args.budget,
            budget_percentage=args.budget_percentage,
            allowed_del=mad.allowed_num_deletions / mad.beta if mad.beta else 0.5,
            lev="individual",
            criterion=criterion,
        )

        ensure_weights(pre_graph)
        ensure_weights(updated_graph)

        results = []
        failures = []
        for det in detectors:
            try:
                pre = compute_partition_timeout(pre_graph, det, args.timeout)
                post = compute_partition_timeout(updated_graph, det, args.timeout)
                pre_id = find_community_id(pre.communities, target_name)
                post_id = find_community_id(post.communities, target_name)
                changed = 1 if (pre_id is not None and post_id is not None and pre_id != post_id) else 0
                results.append((det, pre_id, post_id, changed))
                dr_agg[det].append(changed)
            except Exception as exc:
                failures.append((det, str(exc)))
                global_failures.setdefault(det, []).append(str(exc))

        # Label Entropy across detectors (post assignments)
        post_ids = [r[2] for r in results if r[2] is not None]
        le = 0.0
        if post_ids:
            from collections import Counter

            counts = Counter(post_ids)
            total = sum(counts.values())
            for c in counts.values():
                p = c / total
                le -= p * math.log(p, 2)
        le_list.append(le)

        per_target_results.append(
            {
                "criterion": criterion,
                "target_id": target_id,
                "target_name": target_name,
                "results": results,
                "le": le,
                "failures": failures,
            }
        )

    print("IND Report")
    print(f"Dataset: {args.dataset}")
    print(f"Base detector: {args.base_detector}")
    print(f"Budget: {args.budget} (percentage={args.budget_percentage})")
    print(f"Targets: {', '.join(targets)}")

    for block in per_target_results:
        print(f"\nTarget {block['criterion']} -> id={block['target_id']} name={block['target_name']}")
        print("Per-detector DR (single run):")
        for det, pre_id, post_id, changed in block["results"]:
            print(f"{det}\tpre={pre_id}\tpost={post_id}\tDR={changed}")
        print(f"Label Entropy (LE): {block['le']:.4f}")
        if block["failures"]:
            print("Detectors failed for this target:")
            for det, err in block["failures"]:
                print(f"{det}: {err}")

    # Aggregated DR per detector across targets
    print("\nAggregated DR per detector (avg over targets):")
    for det in detectors:
        if not dr_agg[det]:
            print(f"{det}\tDR=NA")
            continue
        dr_avg = sum(dr_agg[det]) / len(dr_agg[det])
        print(f"{det}\tDR={dr_avg:.4f}")

    avg_le = sum(le_list) / len(le_list) if le_list else 0.0
    print(f"\nAvg Label Entropy (LE) over targets: {avg_le:.4f}")

    if global_failures:
        print("\nDetectors failed (any target):")
        for det, errs in global_failures.items():
            last = errs[-1] if errs else ""
            print(f"{det}: {last}")


if __name__ == "__main__":
    main()
