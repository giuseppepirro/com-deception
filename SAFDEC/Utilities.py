from igraph import *
import numpy as np
import igraph
import copy
import timeit

class Utils:
    def __init__(self, dataset_path):
        self.local_path = dataset_path
        self.graph=None
        self.induced_subgraph=None
        self.target_community=None
        self.communities=None
        self.weights =None
        self.time_detection=None
        self.community_id=None
        self.total_edge_count =None
        self.internal_edge_count =None
        self.external_edge_count =None
        self.ratio_internal_external=None

        ##Compute once community-related information and general information
    def initializeCommunityData(self, detection_algo, recompute_communities):
        if recompute_communities:
            self.communities,self.time_detection = self.computeCommunities(detection_algo)
            self.community_id = self.getTargetCommunityID()
            self.target_community = self.communities[self.community_id]
        self.community_size=len(self.target_community)
        self.induced_subgraph = self.getInducedSubgraph()
        ############################################################
        self.total_edge_count = [self.graph.degree(n) for n in self.target_community]
        #for v in self.target_community:
            #print("CommunityNode=",v)
            #print("InducedNode=",self.convertCommunityNodeTOInduced(v))
            #print("OriginalGraph=",self.convertInducedIDToOriginalGraphID(self.induced_subgraph,self.convertCommunityNodeTOInduced(v)))
           # print(self.convertInducedIDToOriginalGraphID(self.induced_subgraph, v))
        self.internal_edge_count = [self.induced_subgraph.degree(self.convertCommunityNodeTOInduced(n)) for n in self.target_community]
        self.external_edge_count = [a_i - b_i for a_i, b_i in zip(self.total_edge_count, self.internal_edge_count)]
        self.ratio_internal_external = [a_i /b_i for a_i, b_i in zip(self.internal_edge_count, self.external_edge_count)]

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
        self.communities = communities
        return communities

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
        end = timeit.default_timer()
        self.communities=communities
        return communities,(end - start)

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

    def convertCommunityNodeTOInduced(self,communityNode):
        return self.target_community.index(communityNode)

    def convertInducedIDToOriginalGraphID(self, graph, nodeLabel):
        #convert an ID from the induced subgraph to the id of the original graph
        return graph.vs[nodeLabel]["name"]

    def read_network(self, name):
        self.dataset_name=name
        graph=igraph.read(self.local_path + name + ".edgelist", directed=False)
        graph.vs['name']=range(0,graph.vcount())
        #for novel graphs use the lines below for the old the line above
        #data = pd.read_csv(self.local_path + name + ".edgelist")
        #graph = Graph.TupleList(data.itertuples(index=True), directed=False)
        self.graph=graph
        return graph

    def find_nearest(self,array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def getTargetCommunityID(self):
        com_len=[len(self.communities[i])  for i in range(0,len(self.communities))]
        preferred_size=int(np.ceil(max(com_len)*0.5))/2
        closest=self.find_nearest(np.array(com_len), preferred_size)
        target_community_index=com_len.index(closest)
        return target_community_index

    ##The induced subgraph changes the nodes ids
    ## the original node ids are maintained in the attribute "names" of the nodes
    def getInducedSubgraph(self):
        return self.graph.induced_subgraph(self.target_community, implementation="auto")

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

    def getLowestRatioInternalExternalNodes(self):
        # there can be more than one
        nodes_lower_ratio_indexes = np.where(self.ratio_internal_external == np.min(self.ratio_internal_external))[0]
        nodes_lower_ratio =self.target_community[nodes_lower_ratio_indexes]
        return nodes_lower_ratio

    def getBestAdditionSafeness(self):
        candidate_source_add=self.getLowestRatioInternalExternalNodes()
        return 1

    def getBestDeletionSafeness(self):
        bridge_edges_original_ids, bridge_edges = self.getCommunityBridges()
        copy_induced_excluding_bridges = self.getInducedSubgraph()
        copy_induced_excluding_bridges.delete_edges(bridge_edges)





        return 0



def main():
    dataset="lesmis"  #geometry, pubmed, webkb, cora, citeseer,blogcatalog,facebook

    path="/Users/gpirro/Documents/PyCharm_WORKSPACE/CommunityDeception/datasets/edgelist/"
    community_detection_algo="louv"
    sampleInternalPercentage=0.1
    sampleExternalAddPercentage=0.03
    budget_updates=0.2
    ut=Utils(dataset_path=path)


    ut.read_network(dataset)
    ut.initializeCommunityData(detection_algo=community_detection_algo,recompute_communities=True)
    ##
    budget_updates=int(np.ceil(budget_updates*ut.induced_subgraph.ecount()))
    print("Number of communities=", len(ut.communities))
    print("Target community=", ut.target_community)
    print("Target community size=", len(ut.target_community))
    print("Budget=", budget_updates)

    print("Deception score before=", ut.getDeceptionScore())

    ##################@TODO this must be algo specific
    for i in range (0,budget_updates):
        #print("Safeness before=", (ut.getWeightedSafenessScore()))


        print("implement strategy to pick bet uppdate")
       ## ut.getBestUpdateSecrecy(sampleInternalPercentage, sampleExternalAddPercentage)


    ########################
    #recompute the communities and check
    ut.computeCommunities(community_detection_algo)
    print("Deception score after=", ut.getDeceptionScore()) #viene negativo

if __name__ == '__main__':
    main()


