import copy
import random
import heapq
# import timeit
# import time
from math import log
from typing import List

import numpy as np
import sys
from igraph import *

from ndre_counter import count_security_index_by_pre_addition, count_security_index_by_pre_deletion
from ndre_counter import count_pre_security_index

INCOMING = "IN"
OUTGOING = "OUT"


class MAD:
    def __init__(self, runner_utilities):
        self.runner_utilities = runner_utilities
        self.graph = self.runner_utilities.graph
        self.induced_subgraph = None
        self.beta = None
        self.allowed_num_deletions = 0
        self.graph_to_induced_mapping = {}
        self.induced_to_graph_mapping = {}

        self.target_community = self.runner_utilities.target_community
        self.communities = self.runner_utilities.communities
        self.target_community_id = self.runner_utilities.target_community_id

        # related to Individual node deception:
        self.target_node_id = None
        self.dest_community_id = None
        self.dest_community = None

        ## For each i, number of OUTGOING edges toward other nodes in comH
        self.internal_outgoing_edges = {}
        ## For each i, number of INCOMING edges from other nodes in comH
        self.internal_incoming_edges = {}
        ## For each i, number of OUTGOING edges toward nodes outside comH
        self.external_outgoing_edges = {}
        ## For each i, number of INCOMING edges from nodes outside comH
        self.external_incoming_edges = {}
        #
        ###Degree only considering outgoing edges
        self.outgoing_degree = {}
        ###Degree only considering incoming edges
        self.incoming_degree = {}
        ###

        ###Community degree only considering outgoing edges
        self.community_outgoing_degree = {}
        ###Community degree only considering incoming edges
        self.community_incoming_degree = {}
        ###

        ###For each community, the number of edges involving only its members
        self.community_internal_edges = {}

        ###For each community, the number of edges where only one member is in the edge
        self.community_external_edges = {}

        ###Intra-community edges
        self.intra_community_edges = 0

        ###Inter community edges
        self.inter_community_edges = 0

        ##Number of edges
        self.number_of_edges = 0

        self.__degree_distribute = {}
        self.__betweenness_distribute = {}
        self.__outdegree_distribute = {}
        self.__indegree_distribute = {}

        # Weights for label-free/degree-based scoring
        self.omega_d = 1.0
        self.omega_s = 1.0
        self.label_free_score = 0.0
        self.original_max_outdegree = 1.0
        # Local caches for performance
        self._neighbor_cache = {}
        self.ego_nodes = set()
        self.ego_neighbor_cache = {}
        self.ego_outdegree = {}
        self._low_degree_pool = []
        self._neighbor_cache = {}
        # Cache for determining whether communities use vertex names or indices
        self._community_uses_names = None

        # CSD optimizations (top-K candidates + lazy evaluation)
        self._csd_del_heap = []
        self._csd_add_heap = []
        self._csd_pre_count_version_del = 0
        self._csd_pre_count_version_add = 0
        self._csd_candidates_per_source = 3

    def initializeDataStructuresDeceptionALLCommunities(self):
        ###Inter community edges
        # print("initializeDataStructuresDeceptionALLCommunities(self)")
        self.number_of_edges = self.graph.ecount()
        #####The following data structures are needes for EBD2N
        self.__total_edge_set = set(self.graph.get_edgelist())

        self.__partitions_inner_edge_count: List[int] = list()
        self.__partitions_out_volume: List[int] = list()

        outdegree_distribute = self.graph.strength(self.graph.vs, mode='out', weights=None)
        self.__outdegree_distribute = {i: j for i, j in enumerate(outdegree_distribute)}

        indegree_distribute = self.graph.strength(self.graph.vs, mode='in', weights=None)
        self.__indegree_distribute = {i: j for i, j in enumerate(indegree_distribute)}
        # Add name-based keys for robust lookups (cdlib communities often store names)
        if "name" in self.graph.vs.attributes():
            for i, name in enumerate(self.graph.vs["name"]):
                if i in self.__outdegree_distribute:
                    self.__outdegree_distribute[name] = self.__outdegree_distribute[i]
                if i in self.__indegree_distribute:
                    self.__indegree_distribute[name] = self.__indegree_distribute[i]

        # ascending order
        self.__sorted_partitions_by_min_outdegree_node_asc: List[List[int]] = list()
        # self.__sorted_partitions_by_min_indegree_node_asc: List[List[int]] = list()

        self.__sorted_partitions_by_max_outdegree_node_desc: List[List[int]] = list()
        # self.__sorted_partitions_by_max_indegree_node_desc: List[List[int]] = list()
        # self.__sorted_dict = list()

        self.__partitions_num = len(self.communities)
        self.__available_addition_edges = list()
        self.__available_deletion_edges = list()
        # Precompute a low-degree vertex pool for quick sampling (global)
        all_out = self.graph.strength(self.graph.vs, mode='out', weights=None)
        order = sorted(range(len(all_out)), key=lambda i: all_out[i])
        self._low_degree_pool = order[: min(500, len(order))]

        self.__set_necessary_info()

    def initializeDataStructures(self):
        ###### Maps the ids of the nodes in the induced subgraph to the ids of the nodes in the original graph
        self.induced_subgraph = self.runner_utilities.getInducedSubgraph()
        #
        for v in self.induced_subgraph.vs.indices:
            self.graph_to_induced_mapping[self.getInducedGraphNodeId(self.induced_subgraph, v)] = v
            self.induced_to_graph_mapping[v] = self.getInducedGraphNodeId(self.induced_subgraph, v)
        #
        self.internal_outgoing_edges = self.getEdgesInternal(OUTGOING)
        self.internal_incoming_edges = self.getEdgesInternal(INCOMING)
        #
        self.external_outgoing_edges = self.getEdgesExternal(OUTGOING)
        self.external_incoming_edges = self.getEdgesExternal(INCOMING)

        #
        ###Node degree only considering outgoing edges
        self.outgoing_degree = self.getDegree(OUTGOING)
        ###Node degree only considering incoming edges
        self.incoming_degree = self.getDegree(INCOMING)

        ##Community degree
        self.community_outgoing_degree = self.getCommunityDegree(OUTGOING)
        self.community_incoming_degree = self.getCommunityDegree(INCOMING)

        ###Parameter eta (intra-community edges)
        self.intra_community_edges = self.getIntraCommunityEdges()

        self.inter_community_edges = self.number_of_edges - self.intra_community_edges

        ##For each community find the number of edges inside the community
        self.community_internal_edges = self.getCommunityInternalEdges()

        ##For each community find the number of edges that point outside the community
        self.community_external_edges = self.getCommunityExternalEdges()

        ###Inter community edges
        self.number_of_edges = self.graph.ecount()

    def initializeIndividualDeceptionDataStructures(self, criterion="degree"):

        degree_distribute = self.graph.strength(self.graph.vs, mode='all', weights=None)
        self.__degree_distribute = {i: j for i, j in enumerate(degree_distribute)}
        self.original_max_outdegree = max(degree_distribute) if len(degree_distribute) else 1.0

        outdegree_distribute = self.graph.strength(self.graph.vs, mode='out', weights=None)
        self.__outdegree_distribute = {i: j for i, j in enumerate(outdegree_distribute)}
        indegree_distribute = self.graph.strength(self.graph.vs, mode='in', weights=None)
        self.__indegree_distribute = {i: j for i, j in enumerate(indegree_distribute)}
        if "name" in self.graph.vs.attributes():
            for i, name in enumerate(self.graph.vs["name"]):
                self.__outdegree_distribute[name] = self.__outdegree_distribute.get(i, 0)
                self.__indegree_distribute[name] = self.__indegree_distribute.get(i, 0)

        if criterion != "degree":
            betweenness_list = Graph.betweenness(self.graph)
            self.__betweenness_distribute = {i: j for i, j in enumerate(betweenness_list)}

            self.betweenness_degree_distribute = {i: self.__betweenness_distribute[i] + self.__degree_distribute[i] for
                                                  i in self.__betweenness_distribute.keys()}

        self.setPopularTargetNode(basedOn=criterion)

        # Build ego network (2-hop) caches for IND-only operations using target_node_id
        ego_nodes = {self.target_node_id}
        first_hop = set(self.graph.neighbors(self.target_node_id, mode='out')) | set(
            self.graph.neighbors(self.target_node_id, mode='in'))
        ego_nodes |= first_hop
        second_hop = set()
        for n in first_hop:
            second_hop |= set(self.graph.neighbors(n, mode='out')) | set(self.graph.neighbors(n, mode='in'))
        ego_nodes |= second_hop
        self.ego_nodes = ego_nodes
        self.ego_neighbor_cache = {n: set(self.graph.neighbors(n, mode='out')).intersection(self.ego_nodes) for n in ego_nodes}
        self.ego_outdegree = {n: len(neigh) for n, neigh in self.ego_neighbor_cache.items()}

        self.setupDestinationCommunity()

    def __set_necessary_info(self):
        for index, part in enumerate(self.communities):
            part_edge_count = self.graph.induced_subgraph(part, implementation="auto").ecount()

            # out degree equals outVolume - subgraph.edgeCount
            self.__partitions_inner_edge_count.append(part_edge_count)
            self.__partitions_out_volume.append(sum(self.graph.strength(part, mode='out', weights=None)))

        self.sort_partitions()

    def sort_partitions(self):
        """sort self.communities based on min/max nodes using their in/out degree"""

        self.__sorted_partitions_by_min_outdegree_node_asc = sorted(self.communities, key=lambda x: min(
            self.graph.strength(x, mode='out', weights=None)[:len(x)]))
        # self.__sorted_partitions_by_min_indegree_node_asc = sorted(self.communities, key=lambda x: min(self.graph.strength(x, mode='in', weights=None)[:len(x)]))

        self.__sorted_partitions_by_max_outdegree_node_desc = sorted(self.communities, key=lambda x: max(
            self.graph.strength(x, mode='out', weights=None)[:len(x)]), reverse=True)
        # self.__sorted_partitions_by_max_indegree_node_desc = sorted(self.communities, key=lambda x: max(self.graph.strength(x,mode='in', weights=None )[:len(x)]), reverse=True)

    # --------------------------- helper utilities --------------------------- #
    def _uses_name_ids(self) -> bool:
        """Return True if communities store node IDs as names (strings)."""
        if self._community_uses_names is None:
            try:
                first = next(iter(self.target_community))
                self._community_uses_names = isinstance(first, str)
            except Exception:
                self._community_uses_names = False
        return self._community_uses_names

    def _neighbors_in_id_space(self, node, mode="out") -> list:
        """Return neighbors using same ID representation as communities (names or indices)."""
        neigh_idx = self.graph.neighbors(node, mode=mode)
        if self._uses_name_ids() and "name" in self.graph.vs.attributes():
            return [self.graph.vs[i]["name"] for i in neigh_idx]
        return neigh_idx

    def _all_nodes_in_id_space(self) -> list:
        if self._uses_name_ids() and "name" in self.graph.vs.attributes():
            return list(self.graph.vs["name"])
        return list(range(self.graph.vcount()))

    def _update_outdegree(self, node, delta: int) -> None:
        """Update outdegree distribution for both index- and name-based keys."""
        node_id = self._resolve_vertex_id(node)
        new_val = self.__outdegree_distribute.get(node_id, 0) + delta
        self.__outdegree_distribute[node_id] = new_val
        if "name" in self.graph.vs.attributes() and isinstance(node_id, int):
            name = self.graph.vs[node_id]["name"]
            self.__outdegree_distribute[name] = new_val

    def _update_indegree(self, node, delta: int) -> None:
        """Update indegree distribution for both index- and name-based keys."""
        node_id = self._resolve_vertex_id(node)
        new_val = self.__indegree_distribute.get(node_id, 0) + delta
        self.__indegree_distribute[node_id] = new_val
        if "name" in self.graph.vs.attributes() and isinstance(node_id, int):
            name = self.graph.vs[node_id]["name"]
            self.__indegree_distribute[name] = new_val
    def _out_neighbors(self, vid: int) -> set:
        """Return set of out-neighbor vertex ids."""
        # IND mode: restrict to ego network if available
        if self.ego_neighbor_cache and vid in self.ego_neighbor_cache:
            return self.ego_neighbor_cache[vid]
        if vid in self._neighbor_cache:
            return self._neighbor_cache[vid]
        neigh = set(self.graph.neighbors(vid, mode='out'))
        self._neighbor_cache[vid] = neigh
        return neigh

    def _out_degree(self, vid: int) -> float:
        if self.ego_outdegree and vid in self.ego_outdegree:
            return float(self.ego_outdegree[vid])
        if vid in self.__outdegree_distribute:
            return float(self.__outdegree_distribute[vid])
        return float(self.graph.strength(vid, mode='out', weights=None))

    def _jaccard_out(self, a: int, b: int) -> float:
        """Jaccard similarity between out-neighborhoods of a and b."""
        na = self._out_neighbors(a)
        nb = self._out_neighbors(b)
        if not na and not nb:
            return 0.0
        inter = len(na & nb)
        uni = len(na | nb)
        return inter / uni if uni else 0.0

    def _ego_score(self, v: int) -> float:
        """Ego-NDRE score component: sim(u,v) * psi(d^+(v)), psi(d)=log2(1+d)."""
        psi = log(1.0 + self._out_degree(v), 2)
        return self._jaccard_out(self.target_node_id, v) * psi

    def _target_gap_after_edge(self, community_id: int, add: bool, is_intra: bool) -> float:
        """Compute ΔH(C_t) after a single edge edit per paper (Eq. 10)."""
        vol = self.__partitions_out_volume[community_id]
        inner = self.__partitions_inner_edge_count[community_id]
        g = vol - inner
        if add:
            vol_p = vol + 1
            g_p = g if is_intra else g + 1
        else:
            vol_p = vol - 1
            g_p = g if is_intra else g - 1
        if vol_p <= 0:
            return 0.0
        ratio = g_p / vol_p
        if ratio <= 0 or ratio >= 1:
            return 0.0
        return -ratio * log(ratio, 2)

    # --------------------------- CSD optimization helpers --------------------------- #
    def _csd_top_k(self) -> int:
        n = max(self.graph.vcount(), 1)
        return max(1, min(100, max(1, int(0.1 * n))))

    def _csd_candidate_sources(self, mode: str) -> list:
        nodes = list(range(self.graph.vcount()))
        if mode == "del":
            nodes = [n for n in nodes if self._out_degree(n) > 0]
            nodes.sort(key=lambda n: self._out_degree(n), reverse=True)
        else:
            nodes.sort(key=lambda n: self._out_degree(n))
        return nodes[: self._csd_top_k()]

    def _csd_build_candidates(self, mode: str) -> list:
        candidates = set()
        all_nodes = list(range(self.graph.vcount()))
        if not all_nodes:
            return []
        sources = self._csd_candidate_sources(mode)
        max_per_source = max(1, self._csd_candidates_per_source)
        max_tries = min(50, len(all_nodes))
        for u in sources:
            neigh = self._out_neighbors(u)
            if mode == "del":
                if not neigh:
                    continue
                neigh_list = list(neigh)
                random.shuffle(neigh_list)
                for v in neigh_list[:max_per_source]:
                    candidates.add((u, v))
            else:
                picked = 0
                tries = 0
                while picked < max_per_source and tries < max_tries:
                    v = random.choice(all_nodes)
                    tries += 1
                    if v == u or v in neigh:
                        continue
                    candidates.add((u, v))
                    picked += 1
        return list(candidates)

    def _csd_score_edge(self, edge, pre_count, total_degree, addition: bool):
        src = self.runner_utilities.getNodeCommunity(self._resolve_vertex_id(edge[0]))
        des = self.runner_utilities.getNodeCommunity(self._resolve_vertex_id(edge[1]))
        if src == -1 or des == -1:
            return None, None
        src_des = (src, des)
        if addition:
            score = count_security_index_by_pre_addition(
                pre_count,
                edge,
                src_des,
                total_degree + 1,
                self.__partitions_inner_edge_count,
                self.__partitions_out_volume,
                self.__outdegree_distribute,
            )
        else:
            score = count_security_index_by_pre_deletion(
                pre_count,
                edge,
                src_des,
                total_degree - 1,
                self.__partitions_inner_edge_count,
                self.__partitions_out_volume,
                self.__outdegree_distribute,
            )
        return score, src_des

    def _csd_build_heap(self, mode: str):
        addition = mode == "add"
        pre_count = count_pre_security_index(
            self.graph,
            self.communities,
            self.__partitions_inner_edge_count,
            self.__partitions_out_volume,
            self.__outdegree_distribute,
            addition=addition,
        )
        total_degree = self.graph.ecount()
        if addition:
            stamp = self._csd_pre_count_version_add
        else:
            stamp = self._csd_pre_count_version_del

        heap = []
        for edge in self._csd_build_candidates(mode):
            if addition and self.graph.are_connected(*edge):
                continue
            if (not addition) and (not self.graph.are_connected(*edge)):
                continue
            score, src_des = self._csd_score_edge(edge, pre_count, total_degree, addition=addition)
            if score is None:
                continue
            heapq.heappush(heap, (score, edge, src_des, stamp))
        if addition:
            self._csd_add_heap = heap
        else:
            self._csd_del_heap = heap

    def _resolve_vertex_id(self, node):
        """Normalize vertex identifiers across int/name mismatches."""
        if node in self.__outdegree_distribute:
            return node
        try:
            node_int = int(node)
            if node_int in self.__outdegree_distribute:
                return node_int
        except Exception:
            pass
        try:
            return self.graph.vs.find(name=node).index
        except Exception:
            return node

    ##################################################################################################
    ###################################################################################################
    ############################ MAD for hiding a single community #########################

    def __inner_count(self, value):
        if not value:
            return 0
        else:
            return value * log(value, 2)

    def __count_position_entropy_target_community(self, source, community_id, add=True):
        src = source
        src_now = 0
        if add:
            total_degree = self.__partitions_out_volume[community_id] + 1
            src_now = (self.__outdegree_distribute[src] + 1) / total_degree
        else:

            total_degree = self.__partitions_out_volume[community_id] - 1
            src_now = (self.__outdegree_distribute[src] - 1) / total_degree
        src_now = self.__inner_count(src_now)
        pre_position_entropy = 0
        for node in self.target_community:
            if node == src:
                continue

            var = self.__outdegree_distribute[node] / total_degree
            if var > 0:
                pre_position_entropy -= var * log(var, 2)
        return pre_position_entropy - src_now

    def __count_resistance_target_community(self, community_id, add=True):
        # Note that we consider here intra-community edge modification only
        total_degree = self.graph.ecount()
        part_volume = self.__partitions_out_volume[community_id]
        part_inner_edges = self.__partitions_inner_edge_count[community_id]
        inter_outgoing = part_volume - part_inner_edges
        if add:

            resistance = ((inter_outgoing + 1) / (total_degree + 1)) * log((part_volume + 1) / (total_degree + 1), 2)
            return resistance

        else:
            # notice: since we only delete intra, don't modify inter_outgoing
            resistance = (inter_outgoing / (total_degree - 1)) * log((part_volume - 1) / (total_degree - 1), 2)
            return resistance

    ##INTRA-EDGE DELETION
    def getBestIntraEdgeDeletionSingleCommunityNRE(self):
        """SCD rule: only intra-community deletions; choose edge minimizing rho_P(C_t)."""
        min_rho = sys.maxsize
        optimal_edge = None
        edge_partitions = None

        # must have at least one intra edge to delete
        if self.__partitions_inner_edge_count[self.target_community_id] <= 0:
            return min_rho, None, None

        target_set = set(self.target_community)
        # gap is constant across intra deletions; H(C_t) depends on source
        gap = self._target_gap_after_edge(self.target_community_id, add=False, is_intra=True)

        # prioritize high out-degree sources (Theorem 4, but still score numerically)
        candidates = sorted(target_set, key=lambda n: self._out_degree(n), reverse=True)
        max_candidates = min(100, len(candidates))

        for u in candidates[:max_candidates]:
            neighbors = self._neighbors_in_id_space(u, mode="out")
            intra_neighbors = [v for v in neighbors if v in target_set and v != u]
            if not intra_neighbors:
                continue
            Hc = self.__count_position_entropy_target_community(u, community_id=self.target_community_id, add=False)
            if Hc <= 0:
                continue
            rho = gap / Hc
            if rho < min_rho:
                min_rho = rho
                optimal_edge = (u, intra_neighbors[0])
                edge_partitions = (self.runner_utilities.getNodeCommunity(self._resolve_vertex_id(u)),
                                   self.runner_utilities.getNodeCommunity(self._resolve_vertex_id(intra_neighbors[0])))

        if optimal_edge is None:
            return sys.maxsize, None, None
        return min_rho, optimal_edge, edge_partitions

    ###INTER-EDGE addition
    def getBestInterEdgeAdditionSingleCommunityNRE(self):
        """SCD rule: inter-community additions only; choose edge minimizing rho_P(C_t)."""
        min_rho = sys.maxsize
        optimal_edge = None
        edge_partitions = None

        target_set = set(self.target_community)
        all_nodes = self._all_nodes_in_id_space()
        outside_nodes = [n for n in all_nodes if n not in target_set]
        if not outside_nodes:
            return min_rho, None, None

        # gap is constant across inter additions; H(C_t) depends on source
        gap = self._target_gap_after_edge(self.target_community_id, add=True, is_intra=False)

        # prefer low-degree targets outside (weak ties), but score per source
        outside_sorted = sorted(outside_nodes, key=lambda n: self._out_degree(n))

        for u in target_set:
            neighbors = set(self._neighbors_in_id_space(u, mode="out"))
            # find a feasible outside destination
            dest = None
            for v in outside_sorted:
                if v == u or v in neighbors:
                    continue
                dest = v
                break
            if dest is None:
                continue
            Hc = self.__count_position_entropy_target_community(u, community_id=self.target_community_id, add=True)
            if Hc <= 0:
                continue
            rho = gap / Hc
            if rho < min_rho:
                min_rho = rho
                optimal_edge = (u, dest)
                edge_partitions = (self.runner_utilities.getNodeCommunity(self._resolve_vertex_id(u)),
                                   self.runner_utilities.getNodeCommunity(self._resolve_vertex_id(dest)))

        if optimal_edge is None:
            return sys.maxsize, None, None
        return min_rho, optimal_edge, edge_partitions

    def getOrderedCommunity(self, com_id, outDegree=True, reverse=False):

        if outDegree and (not reverse):
            return sorted(self.communities[com_id], key=lambda x: self.graph.strength(x, mode='out', weights=None))
        elif (not outDegree) and (not reverse):
            return sorted(self.communities[com_id], key=lambda x: self.graph.strength(x, mode='in', weights=None))
        elif outDegree and reverse:
            return sorted(self.communities[com_id], key=lambda x: self.graph.strength(x, mode='out', weights=None),
                          reverse=True)
        elif (not outDegree) and reverse:
            return sorted(self.communities[com_id], key=lambda x: self.graph.strength(x, mode='in', weights=None),
                          reverse=True)

    #############################################################################################
    ###################################################################################################
    ############################ MAD for hiding an individual #########################

    def setPopularTargetNode(self, basedOn="degree"):
        # T1: highest degree node in the network
        # T2: node with the highest betweness
        # T3: node having at the same time the highest degree and betweeness

        if basedOn in ("degree", "T3"):
            self.target_node_id = max(self.__degree_distribute, key=self.__degree_distribute.get)
        elif basedOn in ("betweenness", "T2"):
            self.target_node_id = max(self.__betweenness_distribute, key=self.__betweenness_distribute.get)
        elif basedOn in ("between-degree", "combined", "degree+betweenness", "T1"):
            self.target_node_id = max(self.betweenness_degree_distribute, key=self.betweenness_degree_distribute.get)
        else:
            print("Invalid popular node criterion!")
            exit()

    def setupDestinationCommunity(self):
        """Setup destination community, to which the individual node should be moved."""
        target_node_community = self.runner_utilities.getNodeCommunity(self.target_node_id)
        candidate_communities = list(range(len(self.communities)))
        candidate_communities.remove(target_node_community)
        self.dest_community_id = candidate_communities[random.randrange(len(candidate_communities))]
        self.dest_community = self.communities[self.dest_community_id]

    def __count_position_entropy_destination_community(self, source, add=True):
        """NOTE:: This method is used for hiding an individual"""
        src = source
        src_now = 0
        if add:
            total_degree = self.__partitions_out_volume[self.dest_community_id] + self.__outdegree_distribute[
                self.target_node_id] + 1
            src_now = (self.__outdegree_distribute[src] + 1) / total_degree
        else:
            total_degree = self.__partitions_out_volume[self.dest_community_id] + self.__outdegree_distribute[
                self.target_node_id] - 1
            src_now = (self.__outdegree_distribute[src] - 1) / total_degree
        src_now = self.__inner_count(src_now)
        pre_position_entropy = 0
        # H(G) is computed assuming the target node is part of the destination
        self.dest_community.append(self.target_node_id)
        for node in self.dest_community:
            if node == src:
                continue
            var = self.__outdegree_distribute[node] / total_degree
            if var > 0:
                pre_position_entropy -= var * log(var, 2)
        # return the destination to it's original state
        self.dest_community.remove(self.target_node_id)
        return pre_position_entropy - src_now

    def __count_resistance_destination_community(self, add=True):
        """NOTE:: This method is used for hiding an individual"""
        # Note that we consider here intra-community edge modification only
        total_degree = self.graph.ecount()
        part_volume = self.__partitions_out_volume[self.dest_community_id] + self.__outdegree_distribute[
            self.target_node_id]
        self.dest_community.append(self.target_node_id)
        part_inner_edges = self.graph.induced_subgraph(self.dest_community, implementation="auto").ecount()
        self.dest_community.remove(self.target_node_id)
        inter_outgoing = part_volume - part_inner_edges

        if add:
            resistance = (inter_outgoing / (total_degree + 1)) * log((part_volume + 1) / total_degree, 2)
            return resistance
        else:
            resistance = (inter_outgoing / (total_degree - 1)) * log((part_volume - 1) / total_degree, 2)
            return resistance

    def getBestIntraEdgeDeletion_single_node_NRE(self):
        """IND rule: delete edge maximizing sim(u,v)*psi(d^+(v))."""
        neighbors = list(self.graph.neighbors(self.target_node_id, mode='out'))
        if len(neighbors) <= 1:
            return sys.maxsize, None, None

        best_edge = None
        best_score = -1.0
        edge_partitions = None

        for v in neighbors:
            score = self._ego_score(v)
            if score > best_score:
                best_score = score
                best_edge = (self.target_node_id, v)
                # IND assumes no access to global partitions
                edge_partitions = None

        if best_edge is None:
            return sys.maxsize, None, None
        # lower is better for deletion selection logic in callers
        return -best_score, best_edge, edge_partitions

    def getBestInterEdgeAddition_single_node_NRE(self):
        """IND rule: add edge minimizing sim(u,v)*psi(d^+(v)) within the ego graph."""
        neighbors = self._out_neighbors(self.target_node_id)
        # Stay strictly within 2-hop ego view
        ego_pool = self.ego_nodes if self.ego_nodes else set(range(self.graph.vcount()))
        candidate_nodes = (ego_pool - neighbors) - {self.target_node_id}
        if not candidate_nodes:
            return sys.maxsize, None, None

        best_edge = None
        best_score = sys.maxsize
        edge_partitions = None
        for v in candidate_nodes:
            score = self._ego_score(v)
            if score < best_score:
                best_score = score
                best_edge = (self.target_node_id, v)
                # IND assumes no access to global partitions
                edge_partitions = None

        if best_edge is None:
            return sys.maxsize, None, None
        return best_score, best_edge, edge_partitions

    def perform_individual_node_deception(self):
        beta = self.beta
        delBeta = self.allowed_num_deletions
        del_min_security, optimal_del_edge, edge_del_partitions = self.getBestIntraEdgeDeletion_single_node_NRE()
        while (optimal_del_edge != None) and (delBeta > 0):
            # accumulate label-free gain before deletion
            v = optimal_del_edge[1]
            lf_gain = self._ego_score(v)
            self.label_free_score += lf_gain
            self.graph.delete_edges(self.graph.get_eid(*optimal_del_edge))
            # invalidate ego/global neighbor caches for touched nodes
            self._neighbor_cache.pop(self.target_node_id, None)
            self._neighbor_cache.pop(v, None)
            if self.ego_neighbor_cache:
                if self.target_node_id in self.ego_neighbor_cache:
                    self.ego_neighbor_cache[self.target_node_id].discard(v)
                    self.ego_outdegree[self.target_node_id] = len(self.ego_neighbor_cache[self.target_node_id])
            self.__total_edge_set = set(self.graph.get_edgelist())
            if hasattr(self, "_MAD__partitions_out_volume") and edge_del_partitions:
                self.__partitions_out_volume[edge_del_partitions[0]] -= 1
                if edge_del_partitions[0] == edge_del_partitions[1]:
                    if hasattr(self, "_MAD__partitions_inner_edge_count"):
                        self.__partitions_inner_edge_count[edge_del_partitions[0]] -= 1
            self._update_outdegree(optimal_del_edge[0], -1)
            self._update_indegree(optimal_del_edge[1], -1)
            self.sort_partitions()
            delBeta -= 1
            beta -= 1
            ####DELETIONS
            del_min_security, optimal_del_edge, edge_del_partitions = self.getBestIntraEdgeDeletion_single_node_NRE()
        ####ADDITIONS
        add_max_security, optimal_add_edge, edge_add_partitions = self.getBestInterEdgeAddition_single_node_NRE()

        while (optimal_add_edge != None) and (beta > 0):
            v = optimal_add_edge[1]
            lf_penalty = self._ego_score(v)
            self.label_free_score -= lf_penalty
            self.graph.add_edge(*optimal_add_edge)
            self._neighbor_cache.pop(self.target_node_id, None)
            self._neighbor_cache.pop(v, None)
            if self.ego_neighbor_cache:
                if self.target_node_id in self.ego_neighbor_cache:
                    self.ego_neighbor_cache[self.target_node_id].add(v)
                    self.ego_outdegree[self.target_node_id] = len(self.ego_neighbor_cache[self.target_node_id])
            self.__total_edge_set = set(self.graph.get_edgelist())
            self.__total_edge_set.add(optimal_add_edge)
            if hasattr(self, "_MAD__partitions_out_volume") and edge_add_partitions:
                self.__partitions_out_volume[edge_add_partitions[0]] += 1
                if edge_add_partitions[0] == edge_add_partitions[1]:
                    if hasattr(self, "_MAD__partitions_inner_edge_count"):
                        self.__partitions_inner_edge_count[edge_add_partitions[0]] += 1

            self._update_outdegree(optimal_add_edge[0], 1)
            self._update_indegree(optimal_add_edge[1], 1)
            self.sort_partitions()
            beta -= 1
            add_max_security, optimal_add_edge, edge_add_partitions = self.getBestInterEdgeAddition_single_node_NRE()
        return self.graph

    def runEBD2NALL(self, budget, budget_percentage, allowed_del, lev="structure", criterion=None):
        """
        Backwards-compatible alias; experimentsMAD expects runEBD2NALL, which maps to runMAD.
        """
        return self.runMAD(budget, budget_percentage, allowed_del, lev=lev, criterion=criterion)

    def runMAD(self, budget, budget_percentage, allowed_del, lev="structure", criterion=None):

        ######### In case of individual hiding:

        if lev == "individual":
            # We need to return the VARIATION IN the indegree and outdegree and we need to check
            # whether the node changed community
            ## The budget here must be set automatically as a percentage of the DEGREE: DO NOT PASS THE BUDGET
            if budget_percentage:
                self.beta = int(self.__degree_distribute[self.target_node_id] * budget)
            else:
                self.beta = int(budget)
            self.allowed_num_deletions = int(self.beta * allowed_del)
            return self.perform_individual_node_deception()
        ########### In case of structure or community deception:
        beta = int(budget)
        self.beta = beta

        if budget_percentage:  # % of the nubmer of edges
            beta = round((budget * self.graph.ecount()))
            self.beta = beta
        self.allowed_num_deletions = int(self.beta * allowed_del)
        while beta > 0:
            self.perform_best_update(level=lev)
            beta = beta - 1
            # track metrics over time if enabled
            if hasattr(self, "metrics_callback") and callable(self.metrics_callback):
                self.metrics_callback(self.graph, self.communities)
        return self.graph

    #################################################################################################

    def __get_available_edges_for_deletion(self):
        """CONSTRUCT A SET OF CANDIDATE EDGES FOR DELETION"""
        available_edges = list()
        for si in range(self.__partitions_num):
            if self.__partitions_out_volume[si] <= 1:
                continue
            s_order = self.__sorted_partitions_by_max_outdegree_node_desc[si]
            u = max(s_order, key=lambda x: self.__outdegree_distribute[
                x])  # get the highest out-degree node with at least one out-edge
            u_neighbors = self._neighbors_in_id_space(u, mode="out")
            # destination is selected randomly from u neighbors
            v = u_neighbors[random.randrange(len(u_neighbors))]
            available_edges.append((u, v))
        self.__available_deletion_edges = available_edges


    """CONSTRUCT A SET OF CANDIDATE EDGES FOR ADDITION"""
    def __get_available_edges_for_addition(self):
        available_edges = list()
        for si in range(self.__partitions_num):
            s_order = self.__sorted_partitions_by_min_outdegree_node_asc[si]
            u = min(s_order, key=lambda x: self.__outdegree_distribute[x])
            u_neighbors = set(self._neighbors_in_id_space(u, mode="out"))
            for ti in range(self.__partitions_num):
                t_order = random.sample(self.communities[ti], len(self.communities[ti]))
                for v in t_order:
                    if v == u:
                        continue
                    if v not in u_neighbors:
                        available_edges.append((u, v))
                        break
        self.__available_addition_edges = available_edges

    def __choose_edge_deletion(self):
        """CSD rule: choose deletion minimizing global NDRE (rho_P)."""
        if self.graph.ecount() == 0:
            return sys.maxsize, None, None
        self._csd_pre_count_version_del += 1
        stamp = self._csd_pre_count_version_del
        if not self._csd_del_heap:
            self._csd_build_heap("del")
        if not self._csd_del_heap:
            return sys.maxsize, None, None

        pre_count = count_pre_security_index(
            self.graph,
            self.communities,
            self.__partitions_inner_edge_count,
            self.__partitions_out_volume,
            self.__outdegree_distribute,
            addition=False,
        )
        total_degree = self.graph.ecount()

        while self._csd_del_heap:
            score, edge, src_des, edge_stamp = heapq.heappop(self._csd_del_heap)
            if not self.graph.are_connected(*edge):
                continue
            if edge_stamp != stamp:
                score, src_des = self._csd_score_edge(edge, pre_count, total_degree, addition=False)
                if score is None:
                    continue
                heapq.heappush(self._csd_del_heap, (score, edge, src_des, stamp))
                continue
            return score, edge, src_des

        self._csd_build_heap("del")
        if not self._csd_del_heap:
            return sys.maxsize, None, None
        return heapq.heappop(self._csd_del_heap)[:3]

    def __choose_edge_addition(self):
        """CSD rule: choose addition minimizing global NDRE (rho_P)."""
        if self.graph.vcount() < 2:
            return sys.maxsize, None, None
        self._csd_pre_count_version_add += 1
        stamp = self._csd_pre_count_version_add
        if not self._csd_add_heap:
            self._csd_build_heap("add")
        if not self._csd_add_heap:
            return sys.maxsize, None, None

        pre_count = count_pre_security_index(
            self.graph,
            self.communities,
            self.__partitions_inner_edge_count,
            self.__partitions_out_volume,
            self.__outdegree_distribute,
            addition=True,
        )
        total_degree = self.graph.ecount()

        while self._csd_add_heap:
            score, edge, src_des, edge_stamp = heapq.heappop(self._csd_add_heap)
            if self.graph.are_connected(*edge):
                continue
            if edge_stamp != stamp:
                score, src_des = self._csd_score_edge(edge, pre_count, total_degree, addition=True)
                if score is None:
                    continue
                heapq.heappush(self._csd_add_heap, (score, edge, src_des, stamp))
                continue
            return score, edge, src_des

        self._csd_build_heap("add")
        if not self._csd_add_heap:
            return sys.maxsize, None, None
        return heapq.heappop(self._csd_add_heap)[:3]
    ###################################################################################################

    def perform_best_update(self, level="structure"):
        if level == "structure":
            if self.allowed_num_deletions > 0:
                del_min_security, optimal_del_edge, edge_del_partitions = self.__choose_edge_deletion()
            else:
                del_min_security, optimal_del_edge, edge_del_partitions = sys.maxsize, None, None
            add_min_security, optimal_add_edge, edge_add_partitions = self.__choose_edge_addition()
        elif level == "community":
            ####DELETIONS
            if self.allowed_num_deletions > 0:
                del_min_security, optimal_del_edge, edge_del_partitions = self.getBestIntraEdgeDeletionSingleCommunityNRE()
            else:
                del_min_security, optimal_del_edge, edge_del_partitions = sys.maxsize, None, None
            ####ADDITIONS
            add_min_security, optimal_add_edge, edge_add_partitions = self.getBestInterEdgeAdditionSingleCommunityNRE()

        else:
            print("Error: Unrecognized level value!")
            exit()

        if (optimal_del_edge == None) and (optimal_add_edge == None):
            print("Warning: no more edges to delete or add!")
            ######## print the results up to this point
            exit()

        if (add_min_security < del_min_security):
            best_overall_security = add_min_security
            self.graph.add_edge(*optimal_add_edge)
            if level == "structure" and optimal_add_edge in self.__available_addition_edges:
                self.__available_addition_edges.remove(optimal_add_edge)
            self.__total_edge_set = set(self.graph.get_edgelist())
            self.__total_edge_set.add(optimal_add_edge)
            self.__partitions_out_volume[edge_add_partitions[0]] += 1

            if edge_add_partitions[0] == edge_add_partitions[1]:
                self.__partitions_inner_edge_count[edge_add_partitions[0]] += 1
            self._update_outdegree(optimal_add_edge[0], 1)
            self._update_indegree(optimal_add_edge[1], 1)
            self.sort_partitions()
        else:
            best_overall_security = del_min_security
            if self.graph.are_connected(*optimal_del_edge):
                self.graph.delete_edges(self.graph.get_eid(*optimal_del_edge))
            if level == "structure" and optimal_del_edge in self.__available_deletion_edges:
                self.__available_deletion_edges.remove(optimal_del_edge)
            self.__total_edge_set = set(self.graph.get_edgelist())
            self.__partitions_out_volume[edge_del_partitions[0]] -= 1
            if edge_del_partitions[0] == edge_del_partitions[1]:
                self.__partitions_inner_edge_count[edge_del_partitions[0]] -= 1
            self._update_outdegree(optimal_del_edge[0], -1)
            self._update_indegree(optimal_del_edge[1], -1)
            self.sort_partitions()
            self.allowed_num_deletions -= 1
        return best_overall_security

    ###################################################################################################
    ############################    UTILITY METHODS    ##############################

    def getIntraCommunityEdges(self):
        total_intra_edges = 0
        for com in self.communities:
            sg = self.graph.induced_subgraph(com)
            total_intra_edges = total_intra_edges + sg.ecount()
        return total_intra_edges

    def getCommunityDegree(self, edge_direction):
        result_d = {}
        index = 0
        value = 0
        for community in self.communities:
            for node in community:
                if edge_direction == OUTGOING:
                    node_degree = self.graph.degree(node, mode="out")
                    value = value + node_degree
                else:
                    node_degree = self.graph.degree(node, mode='in')
                    value = value + node_degree
            result_d[index] = value
            index = index + 1
            value = 0
        return result_d

    ##not used for modularity
    def getCommunityBridges(self, indG):
        bridge_edges = []
        bridge_edges_original_ids = []
        numCompConn = len(Graph.decompose(indG))
        for e in indG.es:
            copy_induced = copy.deepcopy(indG)
            copy_induced.delete_edges(e)
            new_number = len(copy_induced.decompose())
            if (new_number) > numCompConn:
                bridge_edges.append((indG.vs[e.source].index, indG.vs[e.target].index))
                bridge_edges_original_ids.append((indG.vs[e.source]["name"], indG.vs[e.target]["name"]))
        return bridge_edges_original_ids, bridge_edges

    def getDegree(self, edge_direction):
        result = []
        result_d = {}
        for node_index in range(0, len(self.target_community)):
            if edge_direction == OUTGOING:
                total_degree = self.graph.degree(self.target_community[
                                                     node_index], mode="out")
                result.append(total_degree)
            else:
                total_degree = self.graph.degree(self.target_community[
                                                     node_index], mode='in')
                result.append(total_degree)

        for index, value in enumerate(result):
            result_d[self.target_community[index]] = value
        return result_d

    def getOriginalGraphNodeLabel(self, graph, node):
        return graph.vs[node]["name"]

    def getInducedGraphNodeId(self, graph, nodeLabel):
        # find the id of the node having that specfic nodeLabel
        return graph.vs[nodeLabel]["name"]

    ###For each i in comH the number of nodes in comH that can be reached by only traversing OUTGOING/INCOMING EDGES
    def getInternalReachable(self, edge_direction):
        result = {}
        induced_subgraph = self.induced_subgraph
        ##
        for node in self.target_community:
            if edge_direction == OUTGOING:
                n_reached = induced_subgraph.subcomponent(self.graph_to_induced_mapping[node], mode="out")
            else:
                n_reached = induced_subgraph.subcomponent(self.graph_to_induced_mapping[node], mode="in")
            result[node] = (len(n_reached) - 1)
        return result

    ###For a node i in comH the number of nodes in comH that can be reached by only traversing OUTGOING/INCOMING EDGES
    def getInternalReachableFromNode(self, induced_subgraph, node, edge_direction):
        if edge_direction == OUTGOING:
            n_reached = induced_subgraph.subcomponent(node, mode="out")
        else:
            n_reached = induced_subgraph.subcomponent(node, mode="in")
        return (len(n_reached) - 1)

    ## For each i, number of OUTGOING/INCOMING edges toward/from other nodes in comH
    def getEdgesInternal(self, edge_direction):
        result = np.zeros(len(self.target_community))
        result_d = {}
        for source_node_index in range(0, len(self.target_community)):
            for target_node_index in range(0, len(self.target_community)):
                if edge_direction == OUTGOING:
                    if self.graph.are_connected(self.target_community[
                                                    source_node_index], self.target_community[
                                                    target_node_index]):
                        result[source_node_index] = (result[
                            source_node_index]) + 1
                else:
                    if self.graph.are_connected(self.target_community[
                                                    target_node_index], self.target_community[
                                                    source_node_index]):
                        result[source_node_index] = result[
                                                        source_node_index] + 1

        for index, value in enumerate(result):
            result_d[self.target_community[index]] = value
        return result_d

    ## For each community, number of edges inside the community
    def getCommunityInternalEdges(self):
        result_d = {}
        com_index = 0
        for community in self.communities:
            result_d[com_index] = self.graph.induced_subgraph(community, implementation="auto").ecount()
            com_index = com_index + 1

        return result_d

    ## For each community, number of edges where only one node is inside the community
    def getCommunityExternalEdges(self):
        result = np.zeros(len(self.communities))
        com_index = 0
        result_d = {}
        for community in self.communities:
            # sum of the degree of all nodes minus degree internal
            for node in community:
                # find neighbours and exclude those in coms=H with set difference
                neigs = set(self.graph.neighbors(node, mode="all"))
                external_nodes = list(neigs.difference(set(self.target_community)))
                result[com_index] = (result[
                    com_index]) + len(external_nodes)
            com_index = com_index + 1
        for index, value in enumerate(result):
            result_d[index] = value
        return result_d

    ## For each i, number of INCOMING/OUTGOING edges from/to nodes outside comH
    def getEdgesExternal(self, edge_direction):
        result = []
        result_d = {}
        for node_index in range(0, len(self.target_community)):
            if edge_direction == OUTGOING:
                total_degree = self.graph.degree(self.target_community[
                                                     node_index], mode="out")
                external_degree = total_degree - self.internal_outgoing_edges[
                    self.target_community[
                        node_index]]
                result.append(external_degree)
            else:
                total_degree = self.graph.degree(self.target_community[
                                                     node_index], mode='in')
                external_degree = total_degree - self.internal_incoming_edges[
                    self.target_community[
                        node_index]]
                result.append(external_degree)

        for index, value in enumerate(result):
            result_d[self.target_community[index]] = value
        return result_d
