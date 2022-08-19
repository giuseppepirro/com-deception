from igraph import *
import pandas as pd
from operator import *
import numpy as np
import copy
import random

class DIFFUSER:
    def __init__(self, dataset_path):
        self.local_path = dataset_path
        self.graph=None
        self.induced_subgraph=None
        self.total_weights=None
        self.internal_weights=None
        self.external_weights=None
        self.external_weights_per_community=None
        self.target_community=None
        self.communities=None
        self.weights =None

    def initializeCommunityData(self, detection_algo,dataset, recompute_communities):
        if recompute_communities:
            self.communities = self.computeCommunities(detection_algo)
            self.target_community_id = self.getTargetCommunityID()
            self.target_community = self.communities[self.target_community_id]
        self.community_size=len(self.target_community)
        self.induced_subgraph = self.getInducedSubgraph()
        self.total_weights = self.getNodeTotalWeightedDegrees()
        self.internal_weights = self.getNodeTotalWeightedDegreesInducedSubgraph()
        self.external_weights = list(map(sub, self.total_weights, self.internal_weights))
        self.external_weights_per_community=self.getTotalExternalEdgeWeightPerNodeAndCommunity()
        self.node_sim_matrix= np.load(self.local_path+dataset+"/"+dataset+".similarity_matrix.npy")

    def computeCommunitiesAfterUpdate(self, detection_algo, g):
        weights=g.es['weight']
        if detection_algo == "louv":
            communities = g.community_multilevel(weights)
        elif detection_algo == "opt":
            communities = g.community_optimal_modularity(weights)
        elif detection_algo == "infomap":
            communities = g.community_infomap(weights)
        elif detection_algo == "walk":
            communities = g.community_walktrap().as_clustering(weights)
        elif detection_algo == "greedy":
            communities = g.community_fastgreedy().as_clustering(weights)
        elif detection_algo == "spin":
            communities = g.community_spinglass(weights)
        elif detection_algo == "labelp":
            communities = g.community_label_propagation(weights)
        elif detection_algo == "btw":
            communities = g.community_edge_betweenness().as_clustering(weights)
        elif detection_algo == "eig":
            communities = g.community_leading_eigenvector(weights)
        self.communities = communities
        return communities

    def computeCommunities(self,detection_algo):
        g=self.graph
        weights=g.es['weight']
        if detection_algo == "louv":
            communities = g.community_multilevel(weights)
        elif detection_algo == "opt":
            communities = g.community_optimal_modularity(weights)
        elif detection_algo == "infomap":
            communities = g.community_infomap(weights)
        elif detection_algo == "walk":
            communities = g.community_walktrap().as_clustering(weights)
        elif detection_algo == "greedy":
            communities = g.community_fastgreedy().as_clustering(weights)
        elif detection_algo == "spin":
            communities = g.community_spinglass(weights)
        elif detection_algo == "labelp":
            communities = g.community_label_propagation(weights)
        elif detection_algo == "btw":
            communities = g.community_edge_betweenness().as_clustering(weights)
        elif detection_algo == "eig":
            communities = g.community_leading_eigenvector(weights)
        self.communities=communities
        return communities


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

    def getNodeWeights(nodes, graph):
        return [graph.strength(nodes, mode="all", weights="weight")]

    def getOriginalGraphNodeLabel(self, graph, node):
        return graph.vs[node]["name"]

    def getInducedGraphNodeId(self, graph, nodeLabel):
        #find the id of the node having that specfic nodeLabel
        return graph.vs[nodeLabel]["name"]

        #This has to be called once for each new dataset
    def preprocess_weighted_graph(self, name):
        data = pd.read_csv(self.local_path + name + ".edgelist")
        graph = Graph.TupleList(data.itertuples(index=False), directed=False, weights=True)
        graph.vs['name']=range(0,graph.vcount())
        edges=pd.Series(graph.get_edgelist())
        df = pd.DataFrame.from_records(edges,columns=["source","target"])
        if "weight"  in df.columns:
            ##normalize weights
            normalized_df = (data["weight"] ) / (data["weight"].max())
        else:
            data=np.random.rand(graph.ecount())
            normalized_df =pd.DataFrame(data, columns = ['weight'])
            normalized_df = (normalized_df["weight"] ) / (normalized_df["weight"].max())
        df=df.join(pd.DataFrame(normalized_df))
        df.to_csv(self.local_path+name+".edgelist_ordered_igraph",index_label=False,index=False)
        ##generate node similarity matrix (only for test)
        a = np.random.rand(graph.vcount(), graph.vcount())
        node_similarity_matrix = np.tril(a) + np.tril(a, -1).T
        np.save(self.local_path+name+".similarity_matrix",node_similarity_matrix)

    #read a new dataset
    def read_weighted_graph(self,name):
        self.dataset_name=name
        data = pd.read_csv(self.local_path + name +"/"+name+ ".edgelist")
        self.weights=np.array(data["weight"])
        graph = Graph.TupleList(data.itertuples(index=False), directed=False, weights=True)
        self.graph=graph
        return graph

    #This is used to find the community closest to a given size (value parameter)
    def find_nearest(self,array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    #Find the target community inside a community structure
    #This community is found by looking at the size of all communities
    def getTargetCommunityID(self):
        com_len=[len(self.communities[i]) for i in range(0,len(self.communities))]
        # Avoid a too small or too big size
        preferred_size=int(np.ceil(max(com_len)*0.5))/2
        closest=self.find_nearest(np.array(com_len), preferred_size)
        target_community_index=com_len.index(closest)
        return target_community_index

    ##The induced subgraph changes the nodes ids
    ## the original node ids are maintained in the attribute "names" of the nodes
    def getInducedSubgraph(self):
        return self.graph.induced_subgraph(self.target_community, implementation="auto")

    #This find the total degree of each target community member
    def getNodeTotalWeightedDegrees(self):
        res = self.graph.strength(self.target_community, mode="all", weights="weight")
        res = [0 if math.isnan(i) else i for i in res]
        return res

    #This finds for each node of the target community the total weight toward other communities
    def getTotalExternalEdgeWeightPerNodeAndCommunity(self):
        external_weight_per_com_matrix = np.zeros((len(self.target_community), len(self.communities)))
        for com_index in range(0, len(self.communities)):
            community = self.communities[com_index]
            # we need to exclude the target community
            for community_node_index in range(0, len(community)):
                for member_of_target_index in range(0, len(self.target_community)):
                    if self.graph.are_connected(community[community_node_index], self.target_community[member_of_target_index]):
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
    def getNodeTotalWeightedDegreesInducedSubgraph(self):
        ##Position i of res will contain the degree of the i-th member of the community
        res = self.getInducedSubgraph().strength(range(0, self.community_size), mode="all", weights="weight")
        res = [0 if math.isnan(i) else i for i in res]
        return res

    #Similarity between nodes in terms of feature matrix
    def getNodeSimilarity(self, source_node,target_node):
        return self.node_sim_matrix[source_node][target_node]

    #Check if a node belongs to a given community
    def belongToCommunity(self,node,com_id):
        community=self.communities[com_id]
        if node in community:
            return True
        else:
            return False

############# WEIGHTED DIFFUSION #####################
    def getBestAdditionDiffusion(self, sampleInternalPercentage, sampleExternalAddPercentage):
        total_weights = self.total_weights
        external_weights_per_community=self.external_weights_per_community
        max_external_members_c=[np.max(external_weights_per_community[k]) for k in range(0,len(self.target_community))]
        best_difference_total_external = [
            0 if (total_w == 0 or best_e == 0) else
            (total_w-best_e)/total_w for (total_w,best_e) in zip(total_weights,max_external_members_c)]
        best_difference_total_external_indexes = [b[0] for b in sorted(enumerate(best_difference_total_external), reverse=True,
                                                                  key=lambda i: i[1])]
        # X % sample of community nodes from which an edge addition should originate
        values = int(np.ceil(len(best_difference_total_external_indexes) * sampleInternalPercentage))

        #At this stage subset_add_internal_nodes contains nodes
        # these are indexes of elements
        #values cuts the top-values elements
        subset_add_internal_nodes = best_difference_total_external_indexes[:values]

        #convert into the ids of the nodes of the community
        subset_add_internal_nodes =[self.target_community[nodeInd] for nodeInd in subset_add_internal_nodes ]
        subset_add_internal_node_indexes =[self.target_community.index(nodeInd) for nodeInd in subset_add_internal_nodes ]

        best_triple_add=[]

        best_target_com=-1
        for (community_node,community_node_index) in zip(subset_add_internal_nodes,subset_add_internal_node_indexes):
            com_weights = external_weights_per_community[community_node_index ]
            maximum = np.max(com_weights)
            com_with_most_weight = np.where(com_weights == maximum)[0]
            if len(com_with_most_weight)>1:
                com_with_most_weight=com_with_most_weight[0]

            # X % sample of the graph to serve ad edge endpoints from the community where the node has max weight
            subset_add_external_nodes = random.sample(range(0, len(self.communities[com_with_most_weight])),
                                                      int(np.ceil(len(self.communities[com_with_most_weight]) * sampleExternalAddPercentage)))
            #
            max_potential_edge_weight=0
            for external_node in subset_add_external_nodes:
                if self.graph.are_connected(community_node, external_node) == False:
                    source_target_feature_similarity = self.getNodeSimilarity(community_node,external_node)
                    if source_target_feature_similarity >max_potential_edge_weight:
                        best_source_node=community_node
                        best_target_node=external_node
                        max_potential_edge_weight=source_target_feature_similarity
                        maximum_weight=maximum
                        best_target_com=np.array(com_with_most_weight)

                    best_triple_add.append([best_source_node,best_target_node,max_potential_edge_weight,maximum_weight,best_target_com])
                    source_target_gain=[[[source_target_pot_weight[0],source_target_pot_weight[1]],
                                        ((source_target_pot_weight[2]*(self.total_weights[self.target_community.index(source_target_pot_weight[0])]
                                                    -source_target_pot_weight[3]))/(self.total_weights[self.target_community.index(source_target_pot_weight[0])] *(self.total_weights[self.target_community.index(source_target_pot_weight[0])]
                                                    +source_target_pot_weight[2]))),source_target_pot_weight[4]] for source_target_pot_weight in best_triple_add]

                    source_target_gain_array=np.array(source_target_gain)
                    best_values=source_target_gain_array[:,1]
                    maximum = np.max(best_values)
                    best_pair_index = np.where(best_values == maximum)[0][0]
                    source_add_node=source_target_gain[best_pair_index][0][0]
                    target_add_node=source_target_gain[best_pair_index][0][1]
                    best_target_com=source_target_gain[best_pair_index][2]
                    add_gain=maximum/3/(self.community_size)
        return source_add_node,target_add_node,add_gain,best_target_com

    def getBestExtDeletionDiffusion(self, sampleInternalPercentage):
        total_weights = self.total_weights
        external_weights_per_community=self.external_weights_per_community
        max_external_members_c = [np.max(external_weights_per_community[k]) for k in
                                  range(0, len(self.target_community))]
        best_difference_total_external = [
            0 if (total_w == 0 or best_e == 0) else
            (best_e) / total_w for (total_w, best_e) in zip(total_weights, max_external_members_c)]
        best_difference_total_external_indexes = [b[0] for b in sorted(enumerate(best_difference_total_external), reverse=True,
                                                                  key=lambda i: i[1])]

        # X % sample of community nodes from which an edge addition should originate
        values = int(np.ceil(len(best_difference_total_external_indexes) * sampleInternalPercentage))

        #At this stage subset_add_internal_nodes contains nodes
        # identified as in the induced subgraph.
        subset_ex_del_internal_nodes = best_difference_total_external_indexes[:values]
        #convert into the ids of the nodes of the community
        subset_ex_del_internal_nodes =[self.target_community[nodeInd] for nodeInd in subset_ex_del_internal_nodes ]
        subset_ex_del_internal_node_indexes =[self.target_community.index(nodeInd) for nodeInd in subset_ex_del_internal_nodes ]

        best_triple_ex_del=[]
        for (community_node,community_node_index) in zip(subset_ex_del_internal_nodes,subset_ex_del_internal_node_indexes):
            com_weights = external_weights_per_community[community_node_index]
            max_external_value = np.max(com_weights)
            com_with_most_weight = np.where(com_weights == max_external_value)[0]

            best_edge_w=0
            source_del_node=-1
            target_del_node=-1
            edge_direction=0
            del_gain=0
            edges_id_from_community_node=self.graph.incident(community_node,mode="all")
            edges_from_community_node =self.graph.es[edges_id_from_community_node]

           ## this to only consider nodes that have values different from 0 toward external communities
            if  max_external_value > 0:
                for edge in edges_from_community_node:
                    #print(community_node,edge.source,edge.target)
                    if edge.source==community_node:
                        edge_direction=0
                    elif edge.target==community_node:
                        direction_edge=1
                    if edge_direction==0:
                        belong_to_c=self.belongToCommunity(edge.target,com_with_most_weight) or self.belongToCommunity(edge.target, self.target_community_id)
                    elif edge_direction==1:
                        belong_to_c=self.belongToCommunity(edge.source,com_with_most_weight) or self.belongToCommunity(edge.source, self.target_community_id)

                    #If one of the edge endpoints does not belong to the community then we consider it (recall we are considering INTER-EDGE deletions; one of the endpoints must not be neither in C nor in the community linked with max weight)
                    if not belong_to_c:
                        #we must update only the nodes that are in the community
                        #one of the is outside (we are considering INTER-EDGE deletions)
                        source_target_feature_similarity = edge["weight"]
                        best_source_node = edge.source
                        best_target_node = edge.target
                        if source_target_feature_similarity >best_edge_w:
                            best_edge_w=source_target_feature_similarity
                            max_edge_weight = source_target_feature_similarity
                            if direction_edge==0:
                                #print("direction_edge==0",community_node, best_source_node, best_target_node)
                                best_source_node=edge.source
                                best_target_node=edge.target
                            elif direction_edge==1:
                                #print("direction_edge==1",community_node, best_source_node, best_target_node)
                                best_source_node = edge.source
                                best_target_node = edge.target

                        #we must pay attention to the direction of the edge in order to update the node inside C
                        best_triple_ex_del.append([best_source_node,best_target_node,max_edge_weight,max_external_value,direction_edge])

                        source_target_gain_del=[[[source_target_pot_weight[0],source_target_pot_weight[1]],((source_target_pot_weight[2] *source_target_pot_weight[3])/
                                            (self.total_weights[self.target_community.index(source_target_pot_weight[0])]*(self.total_weights[self.target_community.index(source_target_pot_weight[0])]
                                                            -source_target_pot_weight[2]))),source_target_pot_weight[4]]
                                       for source_target_pot_weight in best_triple_ex_del]

                        source_target_gain_array = np.array(source_target_gain_del)
                        best_values = source_target_gain_array[:, 1]
                        max_external_value = np.max(best_values)
                        best_pair_index = np.where(best_values == max_external_value)[0][0]
                        source_del_node = source_target_gain_del[best_pair_index][0][0]
                        target_del_node = source_target_gain_del[best_pair_index][0][1]
                        edge_direction = source_target_gain_del[best_pair_index][2]
                        del_gain=max_external_value/3/(self.community_size)

        return source_del_node,target_del_node,del_gain,edge_direction

    def getBestIntDeletionDiffusion(self):
        external_weights_per_community = self.external_weights_per_community
        ###### BRIDGE nodes: avoid to delete them and disconnect the target community
        bridge_edges_original_ids, bridge_edges = self.getCommunityBridges()
        copy_induced_excluding_bridges = self.getInducedSubgraph()
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
            total_weights_source = self.total_weights[source_node_index]
            total_weights_target = self.total_weights[target_node_index]
            first_part_gain = ((2* e["weight"])/(3*(len(self.target_community)-1)))
            second_part_gain = (((e["weight"])* maximum_source)/(3*total_weights_source*(total_weights_source-e["weight"])))
            third_part_gain = (((e["weight"])* maximum_target)/(3*total_weights_target*(total_weights_target-e["weight"])))
            del_gain=first_part_gain+second_part_gain+third_part_gain
            if del_gain>best_del_gain:
                best_del_gain=del_gain
                source_del_node=self.target_community[self.target_community.index(self.getOriginalGraphNodeLabel(self.induced_subgraph, e.tuple[0]))]
                target_del_node=self.target_community[self.target_community.index(self.getOriginalGraphNodeLabel(self.induced_subgraph, e.tuple[1]))]
        best_del_gain=best_del_gain / (self.community_size)

        return source_del_node, target_del_node, best_del_gain
############# END WEIGHTED DIFFUSION

    def convertEdgeIdInducedToOriginal(self,e):
        converted_e=tuple((self.getOriginalGraphNodeLabel(self.induced_subgraph, e.source), self.getOriginalGraphNodeLabel(self.induced_subgraph, e.target), e["weight"]))
        return converted_e

    def getDeceptionScore(self):
        number_communities=len(self.communities)
        #number of the targetCommunity members in the various communities
        member_for_community=[]
        for member in self.target_community:
            #1 in the position i if member is in the ith community
            current_community_member=[1 if member in community else 0 for community in self.communities]
            member_for_community.append(current_community_member)
        #Each index of the list (representing the id of a community) reports the number of members of the target community included
        member_for_community = [sum(x) for x in zip(*member_for_community)]
        #ratio of the targetCommunity members in the various communities
        ratio_community_members=[members_for_c/len(com) for (members_for_c,com) in zip(member_for_community,self.communities)]
        ##In how many commmunities are the members of the target spread?
        spread_members=sum([1 if value_per_com>0 else 0 for value_per_com in ratio_community_members])
        second_part = 1 / 2 * ((spread_members - 1) / number_communities) + 1/2 * (1 - sum(ratio_community_members) / spread_members)
        #####
        num_components = len(self.getInducedSubgraph().decompose())
        first_part = 1 - ((num_components - 1) / (self.community_size - 1))
        dec_score =first_part * second_part
        return dec_score
