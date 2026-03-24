from igraph import *
import numpy as np
import igraph
import copy
import leidenalg as la
import timeit
from cdlib import algorithms
from cdlib import viz
import networkx as nx
import pickle
import os
import random
from pathlib import Path


def _resolve_dataset_path(dataset_path: str) -> str:
    """
    Prefer a shared root datasets directory when present; fallback to provided path.
    """
    candidates = []
    env_root = os.getenv("DATASETS_ROOT")
    if env_root:
        candidates.append(Path(env_root))
    # repo root / datasets
    repo_root = Path(__file__).resolve().parents[1]
    candidates.append(repo_root / "datasets")
    # module-local and cwd relative
    candidates.append(Path(__file__).resolve().parent / dataset_path)
    candidates.append(Path.cwd() / dataset_path)
    candidates.append(Path(dataset_path))

    for candidate in candidates:
        if candidate.is_dir():
            return str(candidate.resolve()) + "/"
    # fallback to last
    return str(candidates[-1].resolve()) + "/"

#def __init__(self, graph, edges_sum, detection_func, func_args, interval, partitions, path, **kwargs):


class Utils_DIRECTED:
    def __init__(self, dataset_path):
        self.local_path = _resolve_dataset_path(dataset_path)
        self.graph=None
        self.induced_subgraph=None
        self.target_community=None
        self.target_node=None
        self.destination_community = None
        self.destination_community_id = None
        self.communities=None
        self.communities_object=None
        self.ground_truth_communities=None
        #
        self.communities_after=None
        self.communities_after_object=None
        self.target_community_id=None
        self.time_detection=None
        self.num_communities=0
        self.community_membership_dict=None
        self.post_community_membership_dict=None
        self.node_count=0

    def nodeCommunityDict(self, communities):
        communityDictionary = dict()
        for index,community in enumerate(communities):
            for node in community:
                communityDictionary.update({int(node): index})
        return communityDictionary

    def initializeCommunityData(self, detection_algo,dataset, recompute_communities):
        if recompute_communities:
            self.communities,self.time_detection = self.computeCommunitiesCDLIB(detection_algo)
            self.community_id = self.getTargetCommunityID()
            self.target_community_id =  self.community_id
            self.target_community = self.communities[self.community_id]
            #self.target_node = self.getTargetNode()
            #self.destination_community_id = self.getDestinationCommunityID()
            #self.destination_community = self.communities[self.destination_community_id]
            
        self.community_size=len(self.target_community)
        self.induced_subgraph = self.getInducedSubgraph()
        self.num_communities=len(self.communities)
        self.community_membership_dict = self.nodeCommunityDict(self.communities)
        
    
    def getPostDeceptionTargetIndividualMainAttributes(self, network, communities, node):
        
        inDegree = network.strength(node, mode='in', weights=None)
        outDegree = network.strength(node, mode='out', weights=None)
        
        self.post_community_membership_dict = self.nodeCommunityDict(self.communities_after_object.communities)
        
        currentHostCommunity = self.getPostNodeCommunity(node)
        
        return inDegree, outDegree, currentHostCommunity


    def getPreDeceptionTargetIndividualMainAttributes(self, node):
        
        inDegree = self.graph.strength(node, mode='in', weights=None)
        outDegree = self.graph.strength(node, mode='out', weights=None)
        
        currentHostCommunity = self.getNodeCommunity(node)
        
        return inDegree, outDegree, currentHostCommunity
    
    def initializeCommunityDataPassingCommunityID(self, detection_algo, dataset, recompute_communities, community_id):
        if recompute_communities:
            self.communities, self.time_detection = self.computeCommunitiesCDLIB(detection_algo)
            self.community_id = community_id
            self.target_community_id = self.community_id
            self.target_community = self.communities[self.community_id]
        self.community_size = len(self.target_community)
        self.induced_subgraph = self.getInducedSubgraph()
        self.num_communities=len(self.communities)

    def computeCommunitiesAfterUpdate(self, detection_algo,g):
        if detection_algo == "louv":
            communities = g.community_multilevel()
        elif detection_algo == "opt":
            communities = g.community_optimal_modularity()
        elif detection_algo == "infomap":
            communities = g.community_infomap()
        elif detection_algo == "walk":
            communities = g.community_walktrap().as_clustering()
        elif detection_algo == "greedy":
            communities = g.community_fastgreedy().as_clustering()
        elif detection_algo == "spin":
            communities = g.community_spinglass()
        elif detection_algo == "labelp":
            communities = g.community_label_propagation()
        elif detection_algo == "btw":
            communities = g.community_edge_betweenness().as_clustering()
        elif detection_algo == "eig":
            communities = g.community_leading_eigenvector()
        elif detection_algo == "leiden":
            communities = la.find_partition(g, la.ModularityVertexPartition)
        self.communities = communities
        return communities
    def computeCommunitiesAfterUpdateCDLIB(self,detection_algo,g):
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
            coms = self._compute_gemsec(g)
        elif detection_algo == "walk":
            coms = algorithms.walktrap(g)
        elif detection_algo == "tc":
            coms = algorithms.threshold_clustering(g)
        elif detection_algo == "rb":
            coms = algorithms.rb_pots(g)
        end = timeit.default_timer()
        self.communities_after = coms.communities
        self.communities_after_object=coms
        return coms.communities,(end - start)

    def computeCommunitiesCDLIB(self,detection_algo):
        g = self.graph
        path = self.local_path +self.dataset_name
        com_obj_path = path +"_"+detection_algo+".comms"
        time=0
        coms = None
        if os.path.exists(com_obj_path):
            try:
                print("Reading communities=",detection_algo)
                with open(com_obj_path, "rb") as file:
                    coms= pickle.load(file)
                print("Finished reading communities=",detection_algo)
                print("# of communities found by ", detection_algo, "=", len(coms.communities))
            except Exception:
                coms = None
        if coms is None:
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
            elif detection_algo == "gemsec" or detection_algo == "gem":
                coms = self._compute_gemsec(g)
            elif detection_algo == "tc":
                coms = algorithms.threshold_clustering(g)
            elif detection_algo == "rb":
                coms = algorithms.rb_pots(g)
            end = timeit.default_timer()
            time=end-start
            try:
                with open(com_obj_path, "wb") as file:
                    pickle.dump(coms, file)
            except Exception:
                pass
        print("# of communities found by ",detection_algo,"=",len(coms.communities))
        self.communities = coms.communities
        self.communities_object=coms
        return coms.communities,time

    def _compute_gemsec(self, g):
        # Prefer cdlib wrapper if available; otherwise fall back to karateclub.
        if hasattr(algorithms, "gemsec"):
            return algorithms.gemsec(g)
        try:
            from karateclub import GEMSEC
            from cdlib.classes import NodeClustering
        except Exception as exc:
            raise RuntimeError(
                "GEMSEC not available. Install karateclub or use a different detector."
            ) from exc

        # Convert igraph -> networkx
        if isinstance(g, igraph.Graph):
            H = nx.Graph()
            H.add_nodes_from(range(g.vcount()))
            H.add_edges_from(g.get_edgelist())
        else:
            H = g

        n = max(H.number_of_nodes(), 1)
        clusters = max(2, int(np.sqrt(n)))
        model = GEMSEC(dimensions=64, clusters=clusters)
        model.fit(H)
        memberships = model.get_memberships()
        comm_map = {}
        for node, cid in memberships.items():
            comm_map.setdefault(cid, []).append(node)
        communities = list(comm_map.values())
        return NodeClustering(communities, graph=g, method_name="gemsec-karateclub")

    def computeCommunities(self,detection_algo):
        g=self.graph
        start = timeit.default_timer()
        if detection_algo == "louv":
            communities = g.community_multilevel()
        elif detection_algo == "opt":
            communities = g.community_optimal_modularity()
        elif detection_algo == "infomap":
            communities = g.community_infomap()
        elif detection_algo == "walk":
            communities = g.community_walktrap().as_clustering()
        elif detection_algo == "greedy":
            communities = g.community_fastgreedy().as_clustering()
        elif detection_algo == "spin":
            communities = g.community_spinglass()
        elif detection_algo == "labelp":
            communities = g.community_label_propagation()
        elif detection_algo == "btw":
            communities = g.community_edge_betweenness().as_clustering()
        elif detection_algo == "eig":
            communities = g.community_leading_eigenvector()
        elif detection_algo == "leiden":
            communities = la.find_partition(g, la.ModularityVertexPartition)
        end = timeit.default_timer()
        self.communities=communities
        return communities,(end - start)

    #@Giuseppe: 3.1.2023
    #Bridges should be excluded when deleting!
    def getGraphBridges(self):
        bridge_edges = []
        numCompConn = len(Graph.decompose(self.graph))
        for e in self.graph.es:
            copy_graph = copy.deepcopy(self.graph)
            copy_graph.delete_edges(e)
            new_number = len(copy_graph.decompose())
            if (new_number) > numCompConn:
                bridge_edges.append((self.graph.vs[e.source].index, self.graph.vs[e.target].index))
        return bridge_edges

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
            new_number=len(copy_induced.decompose())
            if (new_number) > numCompConn:
                bridge_edges.append((indG.vs[e.source].index, indG.vs[e.target].index))
                bridge_edges_original_ids.append((indG.vs[e.source]["name"], indG.vs[e.target]["name"]))
        return bridge_edges_original_ids, bridge_edges

    def getModularity(graph, communities):
        return Graph.modularity(graph, communities)

    def getNodeDegrees(nodes, graph):
        return [graph.strength(nodes, mode="all", weights="weight")]

    def getOriginalGraphNodeLabel(self, graph, node):
        return graph.vs[node]["name"]

    def getInducedGraphNodeId(self, graph, nodeLabel):
        #find the id of the node having that specfic nodeLabel
        return graph.vs[nodeLabel]["name"]


    def read_directed_graph(self,name):
        path=self.local_path + name
        obj_path=path+ ".obj"
        file_path=path+ ".edgelist"
        self.dataset_name = name
        if os.path.exists(obj_path):
            print("Reading preprocessed network=",name)
            with open(obj_path, "rb") as file:
                self.graph = pickle.load(file)
        else:
            self.graph=Graph.Read_Ncol(file_path, directed=True)
            with open(obj_path, "wb") as file:
                pickle.dump(self.graph, file)
        self.node_count=self.graph.vcount()
        return self.graph

    def read_ground_truth_communities(self, name):
        fileDict={}

        self.dataset_name = name
        with open(self.local_path + name + ".ground_truth", 'r') as f:
            for line in f:
                first, second = line.split(' ', 1)
                if fileDict.__contains__(int(second.strip())):
                    fileDict[int(second.strip())] += [(first.strip())]
                else:
                    fileDict[int(second.strip())] = [(first.strip())]
        ground_truth_communities = fileDict.values()
        return ground_truth_communities

    def read_unweighted_graph(self, name):
        self.dataset_name=name
        graph=igraph.read(self.local_path + name + ".edgelist", directed=False)
        self.graph=graph
        return graph

    def find_nearest(self,array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def getTargetCommunityID(self):
        com_len=[len(self.communities[i])   for i in range(0,len(self.communities))]
        preferred_size=int(np.ceil(max(com_len)*0.1))/2
        closest=self.find_nearest(np.array(com_len), preferred_size)
        target_community_index=com_len.index(closest)

        return target_community_index

    def getNodeCommunity(self, node):
        #print("Node=",node, " Node_count=",self.node_count)
        if int(node) <= self.node_count:
            return self.community_membership_dict.get(int(node), -1)
        else:
            return -1
        
    def getPostNodeCommunity(self, node):
        #print("Node=",node, " Node_count=",self.node_count)
        if int(node) <= self.node_count:
            return self.post_community_membership_dict.get(int(node), -1)
        else:
            return -1
        
    ##The induced subgraph changes the nodes ids
    ## the original node ids are maintained in the attribute "names" of the nodes
    def getInducedSubgraph(self):
        return self.graph.induced_subgraph(self.target_community, implementation="auto")

    def convertEdgeIdInducedToOriginal(self,e):
        converted_e=tuple((self.getOriginalGraphNodeLabel(self.induced_subgraph, e.source), self.getOriginalGraphNodeLabel(self.induced_subgraph, e.target), e["weight"]))
        return converted_e

    def getDeceptionScore(self,communities):
        number_communities=len(communities)
        #number of the targetCommunity members in the various communities
        member_for_community=[]
        for member in self.target_community:
            #1 in the position i if member is in the ith community
            current_community_member=[1 if member in community else 0 for community in communities]
            member_for_community.append(current_community_member)
        # @TODO: This is ok
        #Each index of the list (representing the id of a community) reports the number of members of the target community included
        member_for_community = [sum(x) for x in zip(*member_for_community)]

        #print("member_for_community=",member_for_community)

        #@TODO: This is ok
        #ratio of the targetCommunity members in the various communities
        ratio_community_members=[members_for_c/len(com) for (members_for_c,com) in zip(member_for_community,communities)]

        #print("ratio_community_members=",ratio_community_members)

        # @TODO: This is ok
        ##In how many commmunities are the members of the target spread?
        spread_members=sum([1 if value_per_com>0 else 0 for value_per_com in ratio_community_members])

        if spread_members == 0:
            return 0, member_for_community
        second_part = 1 / 2 * ((spread_members - 1) / number_communities) + 1/2 * (1 - sum(ratio_community_members) / spread_members)
        num_components = len(self.getInducedSubgraph().decompose(mode=STRONG))
        first_part = 1 if self.community_size <= 1 else 1 - ((num_components - 1) / (self.community_size - 1))
        dec_score = first_part * second_part
        return dec_score, member_for_community


    #Utilities for Jaccard
    def count_common_parts(self,a_partitions, b_partitions):
        result = 0

        for a_module in a_partitions:
            for b_module in b_partitions:
                a_module = set(a_module)
                b_module = set(b_module)

                intersection_length = len(a_module.intersection(b_module))
                if intersection_length > 1:
                    result += intersection_length * (intersection_length - 1) / 2

        return result

    def count_combinations(self,modules):
        result = 0

        for module in modules:
            length = len(module)
            result += length * (length - 1) / 2

        return result

    def getJaccardAndRecallScores(self,a_partitions, b_partitions):

        r = self.count_common_parts(a_partitions, b_partitions)
        p = self.count_combinations(a_partitions)
        q = self.count_combinations(b_partitions)

        u_v_2r = p + q
        jaccard_index = r / (u_v_2r - r)
        recall_index = r / p

        return jaccard_index, recall_index
### End Jaccard
