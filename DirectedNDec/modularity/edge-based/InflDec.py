import copy
import random
from operator import *
import timeit

import numpy as np
import pandas as pd
from igraph import *

from cdlib import algorithms

class InflDec:
    def __init__(self, dataset_path):
        self.local_path = dataset_path
        self.graph = None
        self.induced_subgraph = None
        self.total_influence = None
        self.internal_weights = None
        self.external_weights = None
        self.external_influence_per_community = None
        self.target_community = None
        self.target_community_size = None
        self.communities = None
        self.weights = None
        #
        self.communities_after = None
        self.communities_after_object = None
        self.big_i = None  # total influence in the network
        self.eta = None  # total intra_community influence in the network
        self.theta = None  # total inter_community influence in the network
        self.sigma = None

        #
        self.fairness=None
        self.goodness=None

    def initializeCommunityData(self, detection_algo, dataset, recompute_communities):
        self.big_i = self.getTotalNetworkInfluence()
        if recompute_communities:
            self.communities = self.computeCommunitiesCDLIB(detection_algo)  # returns a cluster of communities (CS)
            self.target_community_id = self.getTargetCommunityID()  # chooses the target community id based on size
            self.target_community = self.communities[self.target_community_id]

        self.eta = self.getEta()
        self.theta = self.big_i - self.eta
        self.sigma = self.getSigma()

        self.target_community_size = len(self.target_community)

        self.induced_subgraph = self.getInducedSubgraph(
            self.target_community)  # induced_subgraph of the target community
        # Notice that we might need to add variables taking in-out degrees into account
        self.total_influence = self.getNodeTotalWeightedDegrees()  # a list of each target node's degree
        self.internal_weights = self.getNodeTotalWeightedDegreesInducedSubgraph(self.target_community)
        self.external_weights = list(
            map(sub, self.total_influence, self.internal_weights))  # total weight - internal weight
        self.external_influence_per_community = self.getTotalExternalEdgeWeightPerNodeAndCommunity()
        self.node_influence_matrix = np.load(self.local_path + dataset + "/" + dataset + ".influence_matrix.npy")

    def getTotalNetworkInfluence(self):
        return np.sum(self.graph.strength(mode="out", weights="weight"))

    def getEta(self):
        et = 0
        res = []
        for com in self.communities:
            res = self.getInducedSubgraph(com).strength(range(0, len(com)), mode="out", weights="weight")
            res = [0 if math.isnan(i) else i for i in res]

            et = et + np.sum(res)

        return et

    def getSigma(self):
        sig = 0
        for com in self.communities:
            sumOutCom = np.sum(self.graph.strength(com, mode="out", weights="weight"))
            sumInCom = np.sum(self.graph.strength(com, mode="in", weights="weight"))
            sig = sig + (sumOutCom * sumInCom)
        return sig

    def getComOutDegree(self, com):
        res = self.graph.strength(com, mode="out", weights="weight")
        res = [0 if math.isnan(i) else i for i in res]
        return sum(res)

    def getComInDegree(self, com):
        res = self.graph.strength(com, mode="in", weights="weight")
        res = [0 if math.isnan(i) else i for i in res]
        return sum(res)

    def computeCommunitiesAfterUpdateCDLIB(self,detection_algo,g):
        #print("after deception=",g.ecount())
        start = timeit.default_timer()
        if detection_algo == "leiden":
            coms = algorithms.leiden(g)
        elif detection_algo == "infomap":
            coms = algorithms.infomap(g)
        elif detection_algo == "dm":
            coms=algorithms.rb_pots(g)
        elif detection_algo == "surprise":
            coms=algorithms.surprise_communities(g)
        elif detection_algo == "ds":
            coms = algorithms.gdmp2(g)
        elif detection_algo == "gemsec":
            coms = algorithms.gemsec(g)
        elif detection_algo == "walk":
            coms = algorithms.walktrap(g)
        elif detection_algo == "tc":
            coms = algorithms.threshold_clustering(g)
        end = timeit.default_timer()
        self.communities_after = coms.communities
        self.communities_after_object=coms
        return coms.communities,(end - start)

    def computeCommunitiesCDLIB(self,detection_algo):
        g = self.graph
        weights = g.es['weight']
        #print("before deception=",g.ecount())
        start = timeit.default_timer()
        if detection_algo == "leiden":
            coms = algorithms.leiden(g)
        elif detection_algo == "infomap":
            coms = algorithms.infomap(g)
        elif detection_algo == "dm":
            coms=algorithms.rb_pots(g)
        elif detection_algo == "surprise":
            coms=algorithms.surprise_communities(g)
        elif detection_algo == "walk":
            coms=algorithms.walktrap(g)
        elif detection_algo == "ds":
            coms = algorithms.gdmp2(g)
        elif detection_algo == "gemsec":
            coms = algorithms.gemsec(g)
        elif detection_algo == "tc":
            coms = algorithms.threshold_clustering(g)
        end = timeit.default_timer()
        self.communities = coms.communities
        self.communities_object=coms
        return coms.communities#,(end - start)

    # Returns bridges (edges that if deleted would disconnect the  graph induced
    ## by communty members. The method returns edges with nodeIDs from the original
    ##graph and nodeIDs from the induced subgraph (recall that the induced subgraph starts ids from 0)
    def getCommunityBridges(self):
        bridge_edges = []
        bridge_edges_original_ids = []
        indG = self.induced_subgraph
        numCompConn = len(Graph.decompose(indG))
        for e in indG.es:
            copy_induced = copy.deepcopy(indG)
            copy_induced.delete_edges(e)
            new_number = len(copy_induced.decompose())
            if (new_number) > numCompConn:
                bridge_edges.append((indG.vs[e.source].index, indG.vs[e.target].index))
                bridge_edges_original_ids.append((indG.vs[e.source]["name"], indG.vs[e.target]["name"]))
        return bridge_edges_original_ids, bridge_edges

    def getNodeWeights(nodes, graph):
        return [graph.strength(nodes, mode="out", weights="weight")]

    def getOriginalGraphNodeLabel(self, graph, node):
        return graph.vs[node]["name"]

    def getInducedGraphNodeId(self, graph, nodeLabel):
        # find the id of the node having that specfic nodeLabel
        return graph.vs[nodeLabel]["name"]

        # This has to be called once for each new dataset
    def preprocess_influence_graph(self, name):
        data = pd.read_csv(self.local_path + name + ".edgelist")
        graph = Graph.TupleList(data.itertuples(index=False), directed=False, weights=True)
        graph.vs['name'] = range(0, graph.vcount())
        edges = pd.Series(graph.get_edgelist())
        df = pd.DataFrame.from_records(edges, columns=["source", "target"])
        if "weight" in df.columns:
            ##normalize weights
            normalized_df = (data["weight"]) / (data["weight"].max())
        else:
            data = np.random.rand(graph.ecount())
            normalized_df = pd.DataFrame(data, columns=['weight'])
            normalized_df = (normalized_df["weight"]) / (normalized_df["weight"].max())
        df = df.join(pd.DataFrame(normalized_df))
        df.to_csv(self.local_path + name + ".edgelist_ordered_igraph", index_label=False, index=False)
        ##generate node similarity matrix (only for test)
        a = np.random.rand(graph.vcount(), graph.vcount())
        node_similarity_matrix = np.tril(a) + np.tril(a, -1).T
        np.save(self.local_path + name + ".similarity_matrix", node_similarity_matrix)

    # Create and initialize a Graph object from the edge list
    def read_weighted_graph(self, name):
        self.dataset_name = name

        data = pd.read_csv(self.local_path + name + "/" + name + ".edgelist")

        self.weights = np.array(data["weight"])

        # Modified
        graph = Graph.TupleList(data.itertuples(index=False), directed=True, weights=True)
        self.graph = graph
        return graph

    # This is used to find the community closest to a given size (value parameter)
    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    # Find the target community inside a community structure
    # This community is found by looking at the size of all communities
    def getTargetCommunityID(self):
        com_len = [len(self.communities[i]) for i in range(0, len(self.communities))]
        # Avoid a too small or too big size
        preferred_size = int(np.ceil(max(com_len) * 0.5)) / 2
        closest = self.find_nearest(np.array(com_len), preferred_size)
        target_community_index = com_len.index(closest)
        return target_community_index

    ##The induced subgraph changes the nodes ids
    ## the original node ids are maintained in the attribute "names" of the nodes
    def getInducedSubgraph(self, g):
        return self.graph.induced_subgraph(g, implementation="auto")

    # This finds the total degree of each target community member
    def getNodeTotalWeightedDegrees(self):
        res = self.graph.strength(self.target_community, mode="all", weights="weight")

        res = [0 if math.isnan(i) else i for i in res]

        return res

    # This finds for each node of the target community the total weight toward other communities
    def getTotalExternalEdgeWeightPerNodeAndCommunity(self):
        external_weight_per_com_matrix = np.zeros((len(self.target_community), len(self.communities)))
        for com_index in range(0, len(self.communities)):
            community = self.communities[com_index]
            # we need to exclude the target community
            for community_node_index in range(0, len(community)):
                for member_of_target_index in range(0, len(self.target_community)):
                    if self.graph.are_connected(community[community_node_index],
                                                self.target_community[member_of_target_index]):
                        if com_index != self.target_community_id:
                            w = self.graph[
                                (community[community_node_index], self.target_community[member_of_target_index])]
                            ['weight']
                            external_weight_per_com_matrix[member_of_target_index, com_index] = \
                                external_weight_per_com_matrix[member_of_target_index, com_index] + w
        return external_weight_per_com_matrix

    ##Here the nodes are considered from 0 to the size of the community
    ## If the community is [4,3,2,5] the nodes in the induced
    ## subgraph will be labeled as [0,1,2,3] with the correspondences
    ## 4-->0, 3-->1, etc.
    def getNodeTotalWeightedDegreesInducedSubgraph(self, g):
        ##Position i of res will contain the degree of the i-th member of the community
        res = self.getInducedSubgraph(g).strength(range(0, len(g)), mode="all", weights="weight")
        res = [0 if math.isnan(i) else i for i in res]
        return res

    # Similarity between nodes in terms of feature matrix
    def getNodeInfluence(self, source_node, target_node):
        return self.node_influence_matrix[source_node][target_node]

    # Check if a node belongs to a given community
    def belongToCommunity(self, node, com_id):
        community = self.communities[com_id]
        if node in community:
            return True
        else:
            return False


    def getBestOutAdditionInflDec(self, sampleInternalPercentage, sampleExternalAddPercentage):
        total_influence = self.total_influence
        external_influence_per_community = self.external_influence_per_community
        max_external_members_c = [np.max(external_influence_per_community[k]) for k in
                                  range(0, len(self.target_community))]
        # the best difference, is the largest difference between total weights
        # and the largest node-community totWeight
        best_difference_total_external = [
            0 if (total_w == 0 or best_e == 0) else
            (total_w - best_e) / total_w for (total_w, best_e) in zip(total_influence, max_external_members_c)]

        # indices of largest target node_com weights in decending order
        best_difference_total_external_indexes = [b[0] for b in
                                                  sorted(enumerate(best_difference_total_external), reverse=True,
                                                         key=lambda i: i[1])]

        # X % sample of community nodes from which an edge addition should originate
        values = int(np.ceil(len(best_difference_total_external_indexes) * sampleInternalPercentage))

        # At this stage subset_add_internal_nodes contains nodes
        # these are indexes of elements
        # values cuts the top-values elements
        subset_add_internal_nodes = best_difference_total_external_indexes[:values]

        # convert into the ids of the nodes of the community
        # This contains the original ids of the target nodes that should initiate the modification
        subset_add_internal_nodes = [self.target_community[nodeInd] for nodeInd in subset_add_internal_nodes]

        subset_add_internal_node_indexes = [self.target_community.index(nodeInd) for nodeInd in
                                            subset_add_internal_nodes]

        best_target_com = -1
        # print(len(total_weights))

        best_source_target_candidates = []
        for (community_node, community_node_index) in zip(subset_add_internal_nodes, subset_add_internal_node_indexes):

            com_weights = external_influence_per_community[community_node_index]

            maximum = np.max(com_weights)

            com_with_most_weight = np.where(com_weights == maximum)[0]

            com_with_most_weight =com_with_most_weight [0]

            #if len(com_with_most_weight) > 1:
            #    com_with_most_weight = com_with_most_weight[0]

            # X % sample of the graph to serve ad edge endpoints from the community where the node has max weight
            subset_add_external_nodes = random.sample(range(0, len(self.communities[com_with_most_weight])),
                                                      int(np.ceil(len(self.communities[
                                                                          com_with_most_weight]) * sampleExternalAddPercentage)))
            #
            max_potential_edge_weight = 0

            external_candidates = []
            for external_node in subset_add_external_nodes:

                if self.graph.are_connected(community_node, external_node) == False:
                    best_source_node = community_node
                    best_target_node = external_node
                    inf_source_target = self.getNodeInfluence(community_node, external_node)
                    maximum_weight = maximum
                    best_target_com = np.array(com_with_most_weight)

                    factor = inf_source_target / ((self.big_i ** 2) * ((self.big_i + inf_source_target) ** 2))
                    add_gain_p1 = self.eta * self.big_i * inf_source_target
                    add_gain_p2 = (self.big_i ** 2) * (self.eta + self.getComOutDegree(
                        self.communities[com_with_most_weight]) + self.getComInDegree(self.target_community))

                    add_gain_p3 = (2 * self.big_i + inf_source_target) * self.sigma

                    add_gain = (add_gain_p1 + add_gain_p2 - add_gain_p3) * factor

                    external_candidates.append([best_source_node,
                                                best_target_node,
                                                add_gain,
                                                best_target_com])

            external_candidates = np.array(external_candidates)



            best_external_idx = np.where(external_candidates[:, 2] == np.max(external_candidates[:, 2]))

            best_external_idx = int(best_external_idx[0])

            best_source_target_candidates.append(external_candidates[best_external_idx, :])

        best_source_target_candidates = np.array(best_source_target_candidates)
        best_choice_idx = np.where(best_source_target_candidates[:, 2] == np.max(best_source_target_candidates[:, 2]))
        best_choice_idx=best_choice_idx[0][0]

     

        source_add_node = int(best_source_target_candidates[best_choice_idx, 0])
        target_add_node = int(best_source_target_candidates[best_choice_idx, 1])

        gain = best_source_target_candidates[best_choice_idx, 2]
        best_target_com = int(best_source_target_candidates[best_choice_idx, 3])

        #print(gain)

        return source_add_node, target_add_node, gain, best_target_com

    def getBestInAdditionInflDec(self, sampleInternalPercentage, sampleExternalAddPercentage):
        """This method returns the best inter-community candidate for addition,
        where the destination node is in the target community, and the source
        node is in an external community."""

        total_weights = self.total_influence

        external_weights_per_community = self.external_influence_per_community

        max_external_members_c = [np.max(external_weights_per_community[k]) for k in
                                  range(0, len(self.target_community))]

        # the best difference, is the largest difference between total weights
        # and the largest node-community totWeight
        best_difference_total_external = [
            0 if (total_w == 0 or best_e == 0) else
            (total_w - best_e) / total_w for (total_w, best_e) in zip(total_weights, max_external_members_c)]

        # indices of largest target node_com weights in decending order
        best_difference_total_external_indexes = [b[0] for b in
                                                  sorted(enumerate(best_difference_total_external), reverse=True,
                                                         key=lambda i: i[1])]

        # X % sample of community nodes from which an edge addition should originate
        values = int(np.ceil(len(best_difference_total_external_indexes) * sampleInternalPercentage))

        # At this stage subset_add_internal_nodes contains nodes
        # these are indexes of elements
        # values cuts the top-values elements
        subset_add_internal_nodes = best_difference_total_external_indexes[:values]

        # convert into the ids of the nodes of the community
        # This contains the original ids of the target nodes that should initiate the modification
        subset_add_internal_nodes = [self.target_community[nodeInd] for nodeInd in subset_add_internal_nodes]

        subset_add_internal_node_indexes = [self.target_community.index(nodeInd) for nodeInd in
                                            subset_add_internal_nodes]

        best_target_com = -1
        # print(len(total_weights))

        best_source_target_candidates = []
        for (community_node, community_node_index) in zip(subset_add_internal_nodes, subset_add_internal_node_indexes):

            com_weights = external_weights_per_community[community_node_index]

            maximum = np.max(com_weights)

            com_with_most_weight = np.where(com_weights == maximum)[0]

            com_with_most_weight =com_with_most_weight [0]

           
            # X % sample of the graph to serve ad edge endpoints from the community where the node has max weight
            subset_add_external_nodes = random.sample(range(0, len(self.communities[com_with_most_weight])),
                                                      int(np.ceil(len(self.communities[
                                                                          com_with_most_weight]) * sampleExternalAddPercentage)))
            #
            max_potential_edge_weight = 0
            external_candidates = []
            for external_node in subset_add_external_nodes:
                if self.graph.are_connected(external_node, community_node) == False:
                    best_source_node = external_node
                    best_target_node = community_node
                    inf_source_target = self.getNodeInfluence(external_node, community_node)
                    maximum_weight = maximum
                    best_target_com = np.array(com_with_most_weight)

                    factor = inf_source_target / ((self.big_i ** 2) * ((self.big_i + inf_source_target) ** 2))
                    add_gain_p1 = self.eta * self.big_i * inf_source_target
                    add_gain_p2 = (self.big_i ** 2) * (
                            self.eta + self.getComOutDegree(self.target_community) + self.getComInDegree(
                        self.communities[com_with_most_weight]))
                    add_gain_p3 = (2 * self.big_i + inf_source_target) * self.sigma

                    add_gain = (add_gain_p1 + add_gain_p2 - add_gain_p3) * factor

                    external_candidates.append([best_source_node,
                                                best_target_node,
                                                add_gain,
                                                best_target_com])

            external_candidates = np.array(external_candidates)

            best_external_idx = np.where(external_candidates[:, 2] == np.max(external_candidates[:, 2]))

            best_external_idx = int(best_external_idx[0])

            best_source_target_candidates.append(external_candidates[best_external_idx, :])

        best_source_target_candidates = np.array(best_source_target_candidates)
        best_choice_idx = np.where(best_source_target_candidates[:, 2] == np.max(best_source_target_candidates[:, 2]))

        best_choice_idx = int(best_choice_idx[0])

        source_add_node = int(best_source_target_candidates[best_choice_idx, 0])
        target_add_node = int(best_source_target_candidates[best_choice_idx, 1])

        gain = best_source_target_candidates[best_choice_idx, 2]
        best_target_com = int(best_source_target_candidates[best_choice_idx, 3])

        return source_add_node, target_add_node, gain, best_target_com

    def getBestIntraDeletionInflDec(self):
        """This method finds the best intra-community edge to be deleted"""
        external_weights_per_community = self.external_influence_per_community
        ###### BRIDGE nodes: avoid to delete them and disconnect the target community
        bridge_edges_original_ids, bridge_edges = self.getCommunityBridges()
        copy_induced_excluding_bridges = self.getInducedSubgraph(self.target_community)
        copy_induced_excluding_bridges.delete_edges(bridge_edges)

        best_del_gain = -1
        # Check each edge inside the community (induced subgraph)
        # excluding bridge edges (that would disconnect the community)
        ######
        for e in copy_induced_excluding_bridges.es:
            source_node_index = e.tuple[0]  # this is the id of the source node in the induced subgraph

            target_node_index = e.tuple[1]  # this is the id of the target node in the induced subgraph\
            com_weights_source = external_weights_per_community[source_node_index]
            maximum_source = np.max(com_weights_source)
            com_weights_target = external_weights_per_community[target_node_index]
            maximum_target = np.max(com_weights_target)
            total_weights_source = self.total_influence[source_node_index]
            total_weights_target = self.total_influence[target_node_index]

            inf_source_target = self.getNodeInfluence(source_node_index, target_node_index)

            factor = inf_source_target / ((self.big_i - inf_source_target) ** 2)
            first_part_gain = (self.eta * inf_source_target) / self.big_i
            second_part_gain = self.theta
            third_part_gain = ((2 * self.big_i - inf_source_target) / self.big_i ** 2) * self.sigma
            fourth_part_gain = self.getComInDegree(self.target_community) + self.getComOutDegree(self.target_community)

            del_gain = factor * (first_part_gain + second_part_gain + third_part_gain - fourth_part_gain)

            if del_gain > best_del_gain:
                best_del_gain = del_gain
                source_del_node = self.target_community[
                    self.target_community.index(self.getOriginalGraphNodeLabel(self.induced_subgraph, e.tuple[0]))]
                target_del_node = self.target_community[
                    self.target_community.index(self.getOriginalGraphNodeLabel(self.induced_subgraph, e.tuple[1]))]

        best_del_gain = best_del_gain / (self.target_community_size)

        return source_del_node, target_del_node, best_del_gain

    def convertEdgeIdInducedToOriginal(self, e):
        converted_e = tuple((self.getOriginalGraphNodeLabel(self.induced_subgraph, e.source),
                             self.getOriginalGraphNodeLabel(self.induced_subgraph, e.target), e["weight"]))
        return converted_e

    def getDeceptionScore(self, communities):
        number_communities = len(communities)
        # number of the targetCommunity members in the various communities
        member_for_community = []
        for member in self.target_community:
            # 1 in the position i if member is in the ith community
            current_community_member = [1 if member in community else 0 for community in communities]
            member_for_community.append(current_community_member)
        # Each index of the list (representing the id of a community) reports the number of members of the target community included
        member_for_community = [sum(x) for x in zip(*member_for_community)]

        # ratio of the targetCommunity members in the various communities
        ratio_community_members = [members_for_c / len(com) for (members_for_c, com) in
                                   zip(member_for_community, communities)]

        # print("ratio_community_members=",ratio_community_members)

        ##In how many commmunities are the members of the target spread?
        spread_members = sum([1 if value_per_com > 0 else 0 for value_per_com in ratio_community_members])


        second_part = 1 / 2 * ((spread_members - 1) / number_communities) + 1 / 2 * (
                    1 - sum(ratio_community_members) / spread_members)
       
        #num_components = len(self.getInducedSubgraph().decompose(mode=STRONG))
        #first_part = 1 - ((num_components - 1) / (self.community_size - 1))
        dec_score =  second_part
        return dec_score, member_for_community


"""Edge Prediction with Fairness Goodness code begins here ..."""
"""which is a modified version of (Kumar et. al. 2016) code"""


def initiliaze_scores(self, G):
    fairness = {}
    goodness = {}

    nodes = G.nodes()
    for node in nodes:
        fairness[node] = 1
        try:
            goodness[node] = G.in_degree(node, weight='weight') * 1.0 / G.in_degree(node)
        except:
            goodness[node] = 0
    return fairness, goodness


    def compute_fairness_goodness(self, G):
        fairness, goodness = self.initiliaze_scores(G)
        nodes = G.nodes()
        iter = 0
        while (iter < 100):
            df = 0
            dg = 0

            for node in nodes:
                inedges = G.in_edges(node, data='weight')
                g = 0
                for edge in inedges:
                    g += fairness[edge[0]] * edge[2]

                try:
                    dg += abs(g / len(inedges) - goodness[node])
                    goodness[node] = g / len(inedges)
                except:
                    pass

            # print('Updating fairness')
            for node in nodes:
                outedges = G.out_edges(node, data='weight')
                f = 0
                for edge in outedges:
                    f += 1.0 - abs(edge[2] - goodness[edge[1]]) / 2.0
                try:
                    df += abs(f / len(outedges) - fairness[node])
                    fairness[node] = f / len(outedges)
                except:
                    pass

            if df < math.pow(10, -6) and dg < math.pow(10, -6):
                break
            iter += 1

        return fairness, goodness


    def getFairnessGoodness(self, name):
        G = nx.DiGraph()

        f = open(self.local_path + name + "/" + name + ".edgelist", "r")
        for l in f:

            ls = l.strip().split(",")
            try:
                G.add_edge(ls[0], ls[1], weight=float(ls[2]))  ## the weight should already be in the range of -1 to 1
            except:
                continue

        f.close()

        return self.compute_fairness_goodness(G)


    def predictEdgeWeight(self, source, target):
        """This method predicts edge weight based on
        FxG score (Kumar et. al. 2016)"""

        source = str(source)
        target = str(target)
        source_fairness = self.fairness[source]
        target_goodness = self.goodness[target]

        return source_fairness * target_goodness


"""Edge Prediction with Fairness Goodness code ends here ..."""

def main():
    d="wiki"
    dataset_path = "./datasets/"+d+"/"
    dataset = d

    path=dataset

    # creating and initializing ExperimentsInlfDec object
    infl_dec_for_data_preprocessing = InflDec(dataset_path)

    infl_dec_for_data_preprocessing.preprocess_influence_graph(path)



if __name__ == '__main__':
    main()
