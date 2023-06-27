#Please treat the code as confidential
from cdlib import evaluation
import igraph as ig
import numpy as np
import copy
import random
import timeit
import math
import sys
import pandas as pd

from cdlib import algorithms
class DirNODEC:
    def __init__(self, dataset_path):
        self.local_path = dataset_path
        self.graph=None
        
        self.induced_subgraph=None
        
        
        #influece of each community (sum of the influences of the nodes)
        self.community_influences=None
        
        self.community_out_influences=None
        self.community_in_influences=None
        
        
        #Detection time
        self.time_detection=None


        self.degi_Ci_inf=None
        self.degi_Ci_outInf=None
        self.degi_Ci_inInf=None
        
        
        self.target_community=None
        self.communities=None
        self.communities_object = None
        self.community_membership_dict=None
        
        self.target_community_sorted_neighbors=None
        
        self.communities_after = None
        self.communities_after_object = None
        self.ratio_community_members=None
        
        
        self.total_network_influence = None #total influence in the network
        self.eta = None #total intra_community influence in the network
        self.sigma = None 
        
        
        
        
         
        

    def getModularity(self,graph, communities):
        return evaluation.newman_girvan_modularity(graph,communities).score
    
    def getModularityDIN(self, communities):
        
        I = self.total_network_influence
        
        formula_p1 = self.eta/I
        
        coms_out_inf = self.getCommunityOutInfluences(communities)
        coms_in_inf = self.getCommunityInInfluences(communities)
        
        formula_p2 = self.getUpdatedSigma(coms_out_inf, coms_in_inf)/(I**2)
        
        return formula_p1 - formula_p2
    
    def initializeCommunityData(self, detection_algo, recompute_communities,worst_case_scenario):
        self.total_network_influence = self.getTotalNetworkInfluence()
        
        if recompute_communities:
            
            self.communities,self.time_detection= self.computeCommunitiesCDLIB(detection_algo)
            
            if worst_case_scenario:
                self.target_community_id = self.getTargetCommunityID()
                self.target_community = self.communities[self.target_community_id]
            else:
                self.target_community =self.getTargetCommunityNOWorstCase()
                self.target_community_id=-1

        self.community_size=len(self.target_community)
        self.induced_subgraph = self.getInducedSubgraph()
        self.community_membership_dict=self.nodeCommunityDict()
        
        
        self.degi_Ci_outInf = self.getTargetNodeOutInfluencePerCommunity()
        self.degi_Ci_inInf = self.getTargetNodeInInfluencePerCommunity()
        self.degi_Ci_inf = self.degi_Ci_outInf+self.degi_Ci_inInf
        
        inf_matrix = np.asarray(self.degi_Ci_inf)
        
        self.target_community_sorted_neighbors = sorted(self.communities, key=lambda x: sum(inf_matrix[:, self.communities.index(x)]), reverse=True)
        
        self.target_community_sorted_neighbors.remove(self.target_community)
        
        self.eta = self.getEta()
        
        self.community_influences=self.getCommunityInfluences(self.communities)
        self.community_out_influences=self.getCommunityOutInfluences(self.communities)
        self.community_in_influences=self.getCommunityInInfluences(self.communities)
        
        self.sigma = self.getUpdatedSigma(self.community_out_influences, self.community_in_influences)

    
    def updateStructures(self):
        """This method updates dynamic structures following
        each node operation."""
        
        self.total_network_influence = self.getTotalNetworkInfluence()
        self.eta = self.getEta()
        
        
        self.degi_Ci_outInf = self.getTargetNodeOutInfluencePerCommunity()
        self.degi_Ci_inInf = self.getTargetNodeInInfluencePerCommunity()
        self.degi_Ci_inf = self.degi_Ci_outInf+self.degi_Ci_inInf
        
        
        self.community_influences = self.getCommunityInfluences(self.communities) 
        self.community_out_influences=self.getCommunityOutInfluences(self.communities)
        self.community_in_influences=self.getCommunityInInfluences(self.communities)
        
        self.sigma = self.getUpdatedSigma(self.community_out_influences, self.community_in_influences)
        
        self.community_membership_dict = self.nodeCommunityDict()
        

    
    def computeCommunitiesCDLIB(self,detection_algo):
        g = self.graph
        start = timeit.default_timer()
        if detection_algo == "leiden":
            coms = algorithms.leiden(g)
        elif detection_algo == "louv":
            coms = algorithms.louvain(g)
        elif detection_algo == "infomap":
            coms = algorithms.infomap(g)
        elif detection_algo == "combo":
            coms=algorithms.pycombo(g)
        elif detection_algo == "eig":
            coms=algorithms.eigenvector(g)
        elif detection_algo == "walk":
            coms=algorithms.walktrap(g)
        elif detection_algo == "scd":
            coms = algorithms.scd(g)
        elif detection_algo == "kcut":
            coms = algorithms.kcut(g)
        elif detection_algo == "rspec":
            coms = algorithms.r_spectral_clustering(g,method="regularized", percentile=20)
        elif detection_algo == "labelp":
            coms = algorithms.label_propagation(g)
        elif detection_algo == "greedy":
            coms = algorithms.greedy_modularity(g)
        end = timeit.default_timer()
        self.communities = coms.communities
        self.communities_object=coms
        return coms.communities,(end - start)

    def computeCommunitiesAfterUpdateCDLIB(self,detection_algo,g):
        start = timeit.default_timer()
        if detection_algo == "leiden":
            coms = algorithms.leiden(g)
        elif detection_algo == "louv":
            coms = algorithms.louvain(g)
        elif detection_algo == "infomap":
            coms = algorithms.infomap(g)
        elif detection_algo == "combo":
            coms = algorithms.pycombo(g)
        elif detection_algo == "eig":
            coms = algorithms.eigenvector(g)
        elif detection_algo == "walk":
            coms = algorithms.walktrap(g)
        elif detection_algo == "scd":
            coms = algorithms.scd(g)
        elif detection_algo == "kcut":
            coms = algorithms.kcut(g)
        elif detection_algo == "rspec":
            coms = algorithms.r_spectral_clustering(g, method="regularized", percentile=20)
        elif detection_algo == "labelp":
            coms = algorithms.label_propagation(g)
        elif detection_algo == "greedy":
            coms = algorithms.greedy_modularity(g)
        end = timeit.default_timer()
        self.communities_after = coms.communities
        self.communities_after_object = coms
        return coms.communities,(end - start)

    def getCommunityDegrees(self):
        com_degs=[]
        for com_index in range(len(self.communities)):
            community = self.communities[com_index]
            com_degs.append(sum(self.graph.degree(community)))
        return np.array(com_degs)
    
    def getCommunityOutDegrees(self):
        com_out_degs=[]
        for com_index in range(len(self.communities)):
            community = self.communities[com_index]
            com_out_degs.append(sum(self.graph.outdegree(community)))
        return np.array(com_out_degs)
    
    def getCommunityInDegrees(self):
        com_in_degs=[]
        for com_index in range(len(self.communities)):
            community = self.communities[com_index]
            com_in_degs.append(sum(self.graph.indegree(community)))
        return np.array(com_in_degs)
    
    def getCommunityInfluences(self, coms):
        com_infs=[]
        for community in coms:
            com_infs.append(sum(self.graph.strength(community, mode='all', weights='weight')))
        return np.array(com_infs)
    
    def getCommunityOutInfluences(self, coms):
        com_out_infs=[]
        for community in coms:
            com_out_infs.append(sum(self.graph.strength(community, mode='out', weights='weight')))
        return np.array(com_out_infs)
    
    def getCommunityInInfluences(self, coms):
        com_in_infs=[]
        for community in coms:
            com_in_infs.append(sum(self.graph.strength(community, mode='in', weights='weight')))
        return np.array(com_in_infs)

    def getTargetNodeDegreePerCommunity(self):
        # print(self.communities.subgraph(1))
        internal_degree_per_com_matrix = np.zeros((len(self.target_community), len(self.communities)))
        for com_index in range(0, len(self.communities)):
            community = self.communities[com_index]
            for community_node_index in range(0, len(community)):
                for member_of_target_index in range(0, len(self.target_community)):
                    if self.graph.are_connected(community[community_node_index], self.target_community[
                        member_of_target_index]) or self.graph.are_connected(
                            self.target_community[member_of_target_index], community[community_node_index]):
                            internal_degree_per_com_matrix[member_of_target_index, com_index] = \
                                internal_degree_per_com_matrix[member_of_target_index, com_index] + 1
        return internal_degree_per_com_matrix
    
    def getTargetNodeOutDegreePerCommunity(self):
        # print(self.communities.subgraph(1))
        internal_out_degree_per_com_matrix = np.zeros((len(self.target_community), len(self.communities)))
        for com_index in range(0, len(self.communities)):
            community = self.communities[com_index]
            for community_node_index in range(0, len(community)):
                for member_of_target_index in range(0, len(self.target_community)):
                    if self.graph.are_connected(
                            self.target_community[member_of_target_index], community[community_node_index]):
                            internal_out_degree_per_com_matrix[member_of_target_index, com_index] = \
                                internal_out_degree_per_com_matrix[member_of_target_index, com_index] + 1
        return internal_out_degree_per_com_matrix
    
    def getTargetNodeInDegreePerCommunity(self):
        # print(self.communities.subgraph(1))
        internal_in_degree_per_com_matrix = np.zeros((len(self.target_community), len(self.communities)))
        for com_index in range(0, len(self.communities)):
            community = self.communities[com_index]
            for community_node_index in range(0, len(community)):
                for member_of_target_index in range(0, len(self.target_community)):
                    if self.graph.are_connected(community[community_node_index], self.target_community[
                        member_of_target_index]):
                            internal_in_degree_per_com_matrix[member_of_target_index, com_index] = \
                                internal_in_degree_per_com_matrix[member_of_target_index, com_index] + 1
        return internal_in_degree_per_com_matrix
    
    
    
    def getTargetNodeOutInfluencePerCommunity(self):
        # print(self.communities.subgraph(1))
        internal_out_inf_per_com_matrix = np.zeros((len(self.target_community), len(self.communities)))
        for com_index in range(0, len(self.communities)):
            community = self.communities[com_index]
            for community_node_index, community_node in enumerate(community):
                for member_of_target_index, member_of_target in enumerate(self.target_community):
                    if self.graph.are_connected(member_of_target
                            , community_node):
                            
                            eid = self.graph.get_eid(*(member_of_target, community_node))
                            internal_out_inf_per_com_matrix[member_of_target_index, com_index] = \
                                internal_out_inf_per_com_matrix[member_of_target_index, com_index] + self.graph.es[eid]["weight"]
        return internal_out_inf_per_com_matrix
    
    def getTargetNodeInInfluencePerCommunity(self):
        # print(self.communities.subgraph(1))
        internal_in_inf_per_com_matrix = np.zeros((len(self.target_community), len(self.communities)))
        for com_index in range(0, len(self.communities)):
            community = self.communities[com_index]
            for community_node_index, community_node in enumerate(community):
                for member_of_target_index, member_of_target in enumerate(self.target_community):
                    if self.graph.are_connected(community_node
                            , member_of_target):
                            
                            eid = self.graph.get_eid(*(community_node, member_of_target))
                            internal_in_inf_per_com_matrix[member_of_target_index, com_index] = \
                                internal_in_inf_per_com_matrix[member_of_target_index, com_index] + self.graph.es[eid]["weight"]
        return internal_in_inf_per_com_matrix


    

    def getCommunityBridges(self):
        bridge_edges = []
        bridge_edges_original_ids = []
        indG = self.induced_subgraph
        numCompConn = len(ig.Graph.decompose(indG))
        for e in indG.es:
            copy_induced = copy.deepcopy(indG)
            copy_induced.delete_edges(e)
            new_number=len(copy_induced.decompose())
            if (new_number) > numCompConn:
                bridge_edges.append((indG.vs[e.source].index, indG.vs[e.target].index))
                bridge_edges_original_ids.append((indG.vs[e.source]["name"], indG.vs[e.target]["name"]))
        return bridge_edges_original_ids, bridge_edges

    def getOriginalGraphNodeLabel(self, graph, node):
        return graph.vs[node]["name"]

    def getInducedGraphNodeId(self, graph, nodeLabel):
        return graph.vs[nodeLabel]["name"]
    
   
    #Create and initialize a Graph object from the edge list
    def read_weighted_network(self,name):
        self.dataset_name=name
        
        
        try:
            data = pd.read_csv(self.local_path +name + ".csv", dtype={'source':str, 'target':str})
            
            
            self.weights=np.array(data["weight"])
            
            graph = ig.Graph.TupleList(data.itertuples(index=False), directed=True, weights=True)
            
            #print(graph)
        except:
            print("Error: can't recognise dataset format!")
            exit()
        
        self.graph=graph
        
        
        return graph
    

    def find_nearest(self,array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def getTargetCommunityID(self):
        com_len=[len(self.communities[i]) for i in range(0,len(self.communities))]
        preferred_size=int(np.ceil(max(com_len)*0.5))/2
        closest=self.find_nearest(np.array(com_len), preferred_size)
        target_community_index=com_len.index(closest)
        return target_community_index

    def getTargetCommunityNOWorstCase(self):
        com_len=[len(self.communities[i])   for i in range(0,len(self.communities))]
        preferred_size = int(np.ceil(max(com_len) * 0.5)) / 2
        number_of_nodes_in_target_com = self.find_nearest(np.array(com_len), preferred_size)
        target_com=random.sample(range(0, self.graph.vcount()), number_of_nodes_in_target_com)
        return target_com

    def getInducedSubgraph(self):
        com_excluding_deleted=self.target_community
        
        com_excluding_deleted[:] = (value for value in com_excluding_deleted if value != -1)
        
        return self.graph.induced_subgraph(com_excluding_deleted, implementation="auto")

    def getTargetNodeDegrees(self):
        degrees=np.zeros(len(self.target_community))
        for i in range(0, len(self.target_community)):
            if self.target_community[i]==-1:
                degrees[i]=0
            else:
                degrees[i] =self.graph.degree((self.target_community[i]), mode="all")
        return degrees
    
    def getTargetNodeInfluences(self):
        influences=np.zeros(len(self.target_community))
        for i in range(0, len(self.target_community)):
            if self.target_community[i]==-1:
                influences[i]=0
            else:
                influences[i] =self.graph.strength((self.target_community[i]), mode="all", weights="weight")
        return influences

    def getNodeTotalDegreesInducedSubgraph(self):
        res = self.getInducedSubgraph().degree(range(0, self.community_size), mode="all")
        res = [0 if math.isnan(i) else i for i in res]
        return res

    def belongToCommunity(self,node,com_id):
        return node in self.communities[com_id]

    def nodeCommunityDict(self):
        communityDictionary = dict()
        index = 0
        for comm in self.communities:
            for nod in comm:
                communityDictionary.update({nod: index})
            index += 1
        return communityDictionary

    def getNodeCommunity(self, node):
        
        return self.community_membership_dict[node]
        
        
    def convertEdgeIdInducedToOriginal(self, e):
        converted_e = tuple((self.getOriginalGraphNodeLabel(self.induced_subgraph, e.source),
                             self.getOriginalGraphNodeLabel(self.induced_subgraph, e.target), e["weight"]))
        return converted_e
    
    
    def getDeceptionScore(self, after_updates):
        
        if after_updates==False:
            community=self.target_community
            communities = self.communities
            number_communities = len(communities)
        else:
            community=self.target_community
            communities = self.communities_after
            number_communities = len(communities)
        member_for_community = []
        
        for member in community:
            current_community_member = [1 if member in community else 0 for community in communities]
            member_for_community.append(current_community_member)
        
        
        
        
        member_for_community = [sum(x) for x in zip(*member_for_community)]
        
        
        ratio_community_members = [members_for_c / len(com) for (members_for_c, com) in
                                   zip(member_for_community, communities)]
        
        
        
        self.ratio_community_members=ratio_community_members
        
        spread_members = sum([1 if value_per_com > 0 else 0 for value_per_com in ratio_community_members])
        
        second_part = 1 / 2 * ((spread_members - 1) / number_communities) + 1 / 2 * (
                    1 - sum(ratio_community_members) / spread_members)
        
        #print("**************")
        #print(self.target_community)
        num_components = len(self.getInducedSubgraph().decompose())
        first_part = 1 - ((num_components - 1) / (self.community_size - 1))

        dec_score = first_part * second_part

        return dec_score
   
    
    def plot_communities(self,communities):
        node_labels = range(0, self.graph.vcount())
        ig.plot(communities, mark_groups=True, vertex_size=20, vertex_label=node_labels)
        
        
    def get_intra_community_influence(self, com):
        return self.getCommunityInducedSubgraph(com).strength(range(0, len(com)), mode="out", weights="weight")
    ########################################### END UTILITIES ##############################

    
    def computeDeltaNodeDeletion(self):
        
        #Cannot delete all target community members
        if len(self.target_community) <= 1:
            deltaDel = -sys.maxsize
            return deltaDel, -1
        
        
        range_nodes= range(0,len(self.target_community))
        #print("...  ",range_nodes)
        node_internal_influence=[self.degi_Ci_inf[index][self.getNodeCommunity(node)] for (index,node) in zip(range_nodes,self.target_community)]
        
        
        influences_of_target_nodes = np.array(self.getTargetNodeInfluences())
        
        
        I = self.total_network_influence
        
        mod_before = (self.eta/I)- (self.sigma/(I**2))
        
        deltaDelDic = dict()
        
        for delCandidateIdx, delCandidate in enumerate(self.target_community):
            
            mod_after_p1 = (self.eta-node_internal_influence[delCandidateIdx])/(I-influences_of_target_nodes[delCandidateIdx]) 
            
            community_out_inf_change = self.degi_Ci_inInf[delCandidateIdx]
            community_out_inf_change[self.target_community_id] = node_internal_influence[delCandidateIdx] + (sum(self.degi_Ci_outInf[delCandidateIdx])-self.degi_Ci_outInf[delCandidateIdx][self.target_community_id])
            
            
            community_in_inf_change = self.degi_Ci_outInf[delCandidateIdx]
            community_in_inf_change[self.target_community_id] = node_internal_influence[delCandidateIdx] + (sum(self.degi_Ci_inInf[delCandidateIdx])-self.degi_Ci_inInf[delCandidateIdx][self.target_community_id])
            
            
            updatedSigma = self.getUpdatedSigma(self.community_out_influences-community_out_inf_change, self.community_in_influences-community_in_inf_change)
            
            mod_after_p2 = updatedSigma/((I-influences_of_target_nodes[delCandidateIdx])**2)
            
            mod_after = mod_after_p1 - mod_after_p2
            
            deltaDelDic[delCandidate] = mod_before - mod_after
        
        
        maxKey = max(deltaDelDic, key=deltaDelDic.get)
        
        return deltaDelDic[maxKey], maxKey 


    def computeDeltaNodeAddition(self):
        
        I = self.total_network_influence
        
        mod_before = (self.eta/I)- (self.sigma/(I**2))
        
        
        addition_list = list()
        
        for Cj in self.target_community_sorted_neighbors:
            
            Cj_index= self.communities.index(Cj)
            
            #setting internal influence of the new node to be equal to 
            #the that of the node with minimum intra_influence in destination
            #community 
            inf_i_Cj= min(self.graph.induced_subgraph(Cj, implementation='auto').strength(weights='weight'))
                          
            inf_i = inf_i_Cj*2
            
            #recording the change that affects the in/out influence of community Cj
            #as a result of adding node i
            
            coms_inf_change = np.zeros(len(self.communities))
           
            external_destination_community_idx = self.target_community_id
            
            coms_inf_change[Cj_index] = 1.5*inf_i_Cj
            coms_inf_change[external_destination_community_idx] = 0.5*inf_i_Cj
            
            #external_source_community_idx = random.choice([i for i in range(len(self.communities)) if i != Cj_index])
            #community_out_inf_change[Cj_index] = 1.5*inf_i_Cj
            #community_out_inf_change[external_destination_community_idx] = 0.5*inf_i_Cj
            
            #all inter-community edges of the new node should connect with the target
            #community to cause the maximum disrupt in its modularity
            
            #community_in_inf_change[Cj_index] = inf_i_Cj+(inf_i-inf_i_Cj)
            #community_in_inf_change[external_destination_community_idx] = (inf_i-inf_i_Cj)/2
            
            
            mod_after_p1 = (self.eta+inf_i_Cj)/(I+inf_i)
            
            updatedSigma = self.getUpdatedSigma(self.community_out_influences+coms_inf_change, self.community_in_influences+coms_inf_change)
            
            mod_after_p2 = updatedSigma/((I+inf_i)**2)
            
            mod_after = mod_after_p1 - mod_after_p2
            
            deltaAddMod = mod_before - mod_after
            
            addition_list.append([deltaAddMod, inf_i_Cj, inf_i, Cj_index, external_destination_community_idx])
            
        
        
        addition_list = np.array(addition_list)
        
        best_index = np.argmax(addition_list[:, 0])
        
        return addition_list[best_index]
    
    
    
    def computeDeltaNodeMoving(self):
        nodes_delta_mov=[]
        
        #Don't move if the target community is a singleton
        if len(self.target_community) <= 1:
            deltaMov = -sys.maxsize
            return deltaMov, [-1, -1, -1, -1]
        
        ordered_spread = self.ratio_community_members
        
        I = self.total_network_influence
        
        mod_before = (self.eta/I)- (self.sigma/(I**2))
    
        for node in self.target_community:
            
            nod_idx_in_target = self.target_community.index(node)
            
            node_current_inf = self.graph.strength(node, mode='all', weights='weight')
            
            current_community=self.getNodeCommunity(node)
            
            node_inf_current_com = self.degi_Ci_inf[nod_idx_in_target][current_community]
            
            new_community = ordered_spread.index(ordered_spread[np.argmin(ordered_spread)])
            
            if new_community != current_community:
                
                new_community_len= len(self.communities[new_community])
                
                node_inf_already_in_new_com = self.degi_Ci_inf[nod_idx_in_target][new_community]
                
                new_inf_new_com = self.community_influences[new_community]/new_community_len - node_inf_already_in_new_com
                
                node_inf_in_new_com = node_inf_already_in_new_com + new_inf_new_com
                
                new_node_inf = (node_current_inf - node_inf_current_com) + new_inf_new_com
                
                
                mod_after_p1 = (self.eta - node_inf_current_com + node_inf_in_new_com) / (I-node_current_inf+new_node_inf)
        
                community_out_inf_change_current = np.zeros(len(self.communities))
                community_out_inf_change_current[current_community] = node_inf_current_com + (sum(self.degi_Ci_outInf[nod_idx_in_target])-self.degi_Ci_outInf[nod_idx_in_target][current_community])
                
                
                community_in_inf_change_current = np.zeros(len(self.communities))
                community_in_inf_change_current[current_community] = node_inf_current_com + (sum(self.degi_Ci_inInf[nod_idx_in_target])-self.degi_Ci_inInf[nod_idx_in_target][current_community])
                
                community_out_inf_change_new = np.zeros(len(self.communities))
                community_out_inf_change_new[new_community] = new_inf_new_com + (sum(self.degi_Ci_outInf[nod_idx_in_target])-self.degi_Ci_outInf[nod_idx_in_target][current_community]-self.degi_Ci_outInf[nod_idx_in_target][new_community])
                
                
                community_in_inf_change_new = np.zeros(len(self.communities))
                community_in_inf_change_new[new_community] = new_inf_new_com + (sum(self.degi_Ci_inInf[nod_idx_in_target])-self.degi_Ci_outInf[nod_idx_in_target][current_community]-self.degi_Ci_inInf[nod_idx_in_target][new_community])
                
                
                updatedSigma = self.getUpdatedSigma(self.community_out_influences-community_out_inf_change_current+community_out_inf_change_new, self.community_in_influences-community_in_inf_change_current+community_in_inf_change_new)
                
                
                mod_after_p2 = updatedSigma/((I-node_current_inf+new_node_inf)**2)
                
                mod_after = mod_after_p1 - mod_after_p2
                
                delta_mov_mod =mod_before - mod_after
                
                nodes_delta_mov.append(np.array([node, delta_mov_mod, new_community, new_inf_new_com]))
            else:
                nodes_delta_mov.append(np.array([node, -1, -1, -1]))
        

        nodes_delta_mov=np.array(nodes_delta_mov)
        best_index=np.argmax(nodes_delta_mov[:,1])
        best_delta=nodes_delta_mov[best_index][1]
        
        
        return best_delta,nodes_delta_mov[best_index]
                

    
    def getTotalNetworkInfluence(self):
        
        return np.sum(self.graph.strength(mode="out", weights="weight"))

    
    def getEta(self):
        et = 0
        res = []
        for com in self.communities:
            
            res = self.getCommunityInducedSubgraph(com).strength(range(0, len(com)), mode="out", weights="weight")
            res = [0 if math.isnan(i) else i for i in res]
            
            et = et + np.sum(res)
        
        return et

    def getSigma(self):
        
        sig = 0.0
        for com in self.communities:
            sumOutCom = np.sum(self.graph.strength(com, mode="out", weights="weight"))
            sumInCom = np.sum(self.graph.strength(com, mode="in", weights="weight"))
            
            sig = sig + (sumOutCom * sumInCom)
        
        return sig
    
    def getUpdatedSigma(self, comOutInf, comInInf):
        sig = 0.0
        for idx in range(len(comOutInf)):
            sig = sig + (comOutInf[idx] * comInInf[idx])
        
        return sig
    
    ##The induced subgraph changes the nodes ids
    ## the original node ids are maintained in the attribute "names" of the nodes
    def getCommunityInducedSubgraph(self, g):
        return self.graph.induced_subgraph(g, implementation="auto")

    