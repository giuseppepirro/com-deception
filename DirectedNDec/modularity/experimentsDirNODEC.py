#Please treat the code as confidential

from DirNODEC import DirNODEC
from datetime import datetime
import random
import timeit
class ExperimentsDirNODEC:
    def __init__(self, dataset_path, datasets, detection_algos, budget_percentages, sampleInternalPercentage,
                 sampleExternalAddPercentage, rounds,worst_case_scenario):
        self.dataset_path = dataset_path
        self.datasets = datasets
        self.detection_algos = detection_algos
        self.budget_updates = budget_percentages
        self.sampleInternalPercentage = sampleInternalPercentage
        self.sampleExternalAddPercentage = sampleExternalAddPercentage
        self.rounds = rounds
        self.nodec = DirNODEC(dataset_path=self.dataset_path)
        exp_file_name = "./results/"+ "nodec_"+datetime.today().strftime('%Y-%m-%d-%H:%M:%S') + ".experiments"
        self.file_experiments = open(exp_file_name, "w+")
        self.current_algo = None
        self.worst_case_scenario=worst_case_scenario
        self.experiment_configuration=False


    def com_decept_nodec(self, budget):
        beta = budget

        while (True):
            ##this gives the best value the one giving hghest DELTA_del
            delta_del,best_node_to_delete=self.nodec.computeDeltaNodeDeletion()
            print("Best Del Delta=",delta_del, " , Best Node=", best_node_to_delete)
            
            #this gives a single value as we are adding a new node
            delta_add,inf_i_Cj,inf_i,Cj_index=self.nodec.computeDeltaNodeAddition()
            
            """
            #print("Addition details=> ", delta_add, inf_i_Cj, inf_i, Cj_index)
            print(self.nodec.graph.summary())
            print(self.nodec.graph)
            print(self.nodec.communities)
            self.performNewNodeAddition(inf_i_Cj, inf_i, Cj_index)
            print("&&&&&&&&&&&&&&&&&&&&&&&&&")
            print("added node: ", self.nodec.graph.vcount())
            print(self.nodec.graph.summary())
            print(self.nodec.graph)
            print(self.nodec.communities)
            exit()
            #print("Best Add Delta=",delta_add, " , Internal node influence= ", inf_i_Cj, "Total node influence= ", inf_i)
            """
            ##this gives a value for each community member
            delta_move,move_details= self.nodec.computeDeltaNodeMoving()
            #print(delta_move, move_details)
            
            #print("#################")
            #print("Best Mov Delta=",delta_move, " , Move_details= ", move_details)
            #print("#################")
            
            #deltas=[delta_add,delta_del,delta_move]
            #best_index=deltas.index(max(deltas))
            #print(deltas)
            #print(best_index)
            
            
            if float(delta_move)>0:
                
                self.performNodeMove(move_details)
                print("Move")
            else:
                if delta_del > delta_add:
                    
                    self.performNodeDeletion(best_node_to_delete)
                    
                    self.performNewNodeAddition(inf_i_Cj, inf_i, Cj_index)
                    
                else:
                    self.performNewNodeAddition(inf_i_Cj, inf_i, Cj_index)
                    
            beta = beta - 1
            if (beta == 0) or self.nodec.target_community.count(-1)==self.nodec.initial_community_size-2: ##if the community contains all -1 then stop !
                break
        
        return self.nodec.graph

    def performNodeDeletion(self, node_to_delete):
        out_neighs=self.nodec.graph.neighborhood(vertices=node_to_delete, order=1, mode="out", mindist=1)
        in_neighs=self.nodec.graph.neighborhood(vertices=node_to_delete, order=1, mode="in", mindist=1)
        
        
        for n in out_neighs:
            self.nodec.graph.delete_edges([(node_to_delete,n)])
            
        for n in in_neighs:
            self.nodec.graph.delete_edges([(n, node_to_delete)])
            
        
        nodex_index_in_com=self.nodec.target_community.index(node_to_delete)

        

        node_com=self.nodec.getNodeCommunity(node_to_delete)
        #self.nodec.target_community[nodex_index_in_com]=-1
        copy_coms=self.nodec.communities
        #print(max(self.nodec.community_membership_dict.keys()))
        
        
        if node_to_delete<max(self.nodec.community_membership_dict.keys()):
            if len(copy_coms[node_com])>1:#this avoids to have zero-nodes communities
                copy_coms[node_com].remove(node_to_delete)
                
        
                
        self.nodec.communities=copy_coms
        degree_in_com_copy=self.nodec.degi_Ci
        node_gree_in_its_com=degree_in_com_copy[nodex_index_in_com, node_com]
        
        
        self.nodec.degi_Ci = self.nodec.getTargetNodeDegreePerCommunity()
        
        self.nodec.degi_Ci_outInf = self.nodec.getTargetNodeOutInfluencePerCommunity()
        self.nodec.degi_Ci_inInf = self.nodec.getTargetNodeInInfluencePerCommunity()
        self.nodec.degi_Ci_inf = self.nodec.degi_Ci_outInf+self.nodec.degi_Ci_inInf
        
        self.nodec.community_degrees[node_com]=self.nodec.community_degrees[node_com]-node_gree_in_its_com
        
        self.nodec.E_Ci[node_com]=self.nodec.E_Ci[node_com]-node_gree_in_its_com
        
        self.nodec.graph.delete_vertices(node_to_delete)
        
        
        
        

    def performNewNodeAddition(self, degree_i_in_dest_com, total_degree_i, dest_com_index):
        new_node_id = self.nodec.graph.vcount()
        new_node_name=str(new_node_id)
        
        self.nodec.graph.add_vertices(1)
        self.nodec.graph.vs[new_node_id]["name"]= new_node_name
        
        dest_com=self.nodec.communities[dest_com_index]
        dest_nodes_indexes=random.sample(range(0, len(dest_com)), degree_i_in_dest_com)
        
        
        for i in dest_nodes_indexes:
            self.nodec.graph.add_edges([(new_node_name, dest_com[i])])
        
        #self.nodec.target_community.append(new_node_id)
        
        self.nodec.communities[dest_com_index].append(new_node_name)
        
        
        self.nodec.degi_Ci = self.nodec.getTargetNodeDegreePerCommunity()
        
        self.nodec.degi_Ci_outInf = self.nodec.getTargetNodeOutInfluencePerCommunity()
        self.nodec.degi_Ci_inInf = self.nodec.getTargetNodeInInfluencePerCommunity()
        self.nodec.degi_Ci_inf = self.nodec.degi_Ci_outInf+self.nodec.degi_Ci_inInf
        
        
        
        coms_degs_copy = self.nodec.community_degrees
        
        
        coms_degs_copy[dest_com_index] = coms_degs_copy[dest_com_index] + degree_i_in_dest_com
        self.nodec.community_degrees = coms_degs_copy
        
        
        internal_edges_copy = self.nodec.E_Ci
        internal_edges_copy[dest_com_index] = internal_edges_copy[dest_com_index] + degree_i_in_dest_com
        self.nodec.E_Ci = internal_edges_copy
        
        
        return self.nodec.graph

    def performNodeMove(self, move_details):
        node_to_move= move_details[0] #node to move from Ci in Cj
        new_community=int(move_details[2]) #  Cj
        old_community=self.nodec.getNodeCommunity(node_to_move) #Ci
        new_edges_new_com=float(move_details[3]) #total edges in Cj
        
        edges_to_be_deleted=float(move_details[4]) #edges in in Ci
        node_deg_already_new_com=float(move_details[5]) #nodes tha alredy has in Cj
        nodex_index_in_com=self.nodec.target_community.index(node_to_move)

        target_com=self.nodec.communities[new_community]
        target_nodes_indexes = range(0, len(target_com))

        ## EDGE ADDITIONS in Cj
        for i in target_nodes_indexes:
            if not (self.nodec.graph.are_connected(node_to_move, target_com[i]) and self.nodec.graph.are_connected(target_com[i],node_to_move)):
                self.nodec.graph.add_edges([(node_to_move, target_com[i])])

        ## EDGE DELETIONS in Ci
        neighs = self.nodec.graph.neighborhood(vertices=node_to_move, order=1, mode="all", mindist=1)
        for n in neighs:
            if old_community==self.nodec.getNodeCommunity(str(n)):
                if (self.nodec.graph.are_connected(node_to_move, n) ):
                    self.nodec.graph.delete_edges([(node_to_move, n)])
                else:
                    if (self.nodec.graph.are_connected(n,node_to_move)):
                        self.nodec.graph.delete_edges([(n,node_to_move)])
        degree_in_com_copy = self.nodec.degi_Ci
        node_gree_in_its_com = degree_in_com_copy[nodex_index_in_com, old_community]

        degree_in_com_copy[nodex_index_in_com,old_community] = 0##its degree becomes 0 in Ci

        degree_in_com_copy[nodex_index_in_com, new_community] = new_edges_new_com-node_deg_already_new_com#degree_in_com_copy[nodex_index_in_com, new_community]

        self.nodec.degi_Ci = degree_in_com_copy

        coms_degs_copy = self.nodec.community_degrees
        coms_degs_copy[old_community] = coms_degs_copy[old_community] - node_gree_in_its_com
        coms_degs_copy[new_community] = coms_degs_copy[new_community] +(new_edges_new_com-node_deg_already_new_com)
        self.nodec.community_degrees = coms_degs_copy

        internal_edges_copy = self.nodec.E_Ci
        internal_edges_copy[old_community] = internal_edges_copy[old_community] - node_gree_in_its_com
        internal_edges_copy[new_community] = internal_edges_copy[new_community] + (new_edges_new_com-node_deg_already_new_com)
        self.nodec.E_Ci = internal_edges_copy
        return self.nodec.graph




    def runExperiment(self, network, detection_algo,budget,internalnalP,externalP, worst_case_scenario):
        current_experiment_configuration=network+"_"+detection_algo
        #self.nodec.read_network(network)
        #self.nodec.read_directed_graph(network)
        self.nodec.read_weighted_network(network)
        if current_experiment_configuration != self.experiment_configuration:
            recompute_communities = True
            
            
            self.nodec.initializeCommunityData(detection_algo=detection_algo,
                                               recompute_communities=recompute_communities,
                                               worst_case_scenario=worst_case_scenario)
            
            self.experiment_configuration = current_experiment_configuration

        coms_before = self.nodec.communities
        num_coms_before = len(coms_before)
        print("coms_before = ", coms_before)
        deception_before = self.nodec.getDeceptionScore(after_updates=False)#after_updates=False
        #deception_before = 0 #to be REMOVED
        print("Initial ComH=", self.nodec.target_community)
        print("Initial modularity=", self.nodec.getModularity(self.nodec.graph,self.nodec.communities_object))
        print("Initial decScore=", deception_before)
        print("Initial number of nodes and edges=", self.nodec.graph.ecount(),self.nodec.graph.vcount())
        
        ##### BEGIN DirNODEC HERE
        start = timeit.default_timer()
        new_graph=self.com_decept_nodec(budget)
        # communities in the updated graph
        g = new_graph
        stop = timeit.default_timer()
        time=stop - start
        ##### END DirNODEC HERE

        coms_after = self.nodec.computeCommunitiesAfterUpdateCDLIB(detection_algo, g)
        num_coms_after = len(coms_after)
        deception_after= self.nodec.getDeceptionScore(after_updates=True)#after_updates=True
        #deception_after = 0 #to be REMOVED
        print("coms_after = ", coms_after)
        exit()
        print("Final ComH=", self.nodec.target_community)
        print("Final modularity=", self.nodec.getModularity(g,self.nodec.communities_after_object))
        print("Final number of nodes and edges=", new_graph.ecount(),new_graph.vcount())

        self.file_experiments.write(
           "\n" + network + "\t" + detection_algo + "\t" + str(budget) + "\t" + str(deception_before) + "\t" +
           str(deception_after) + "\t" +
           str(time) + str(self.nodec.target_community_id) + "\t" + str(self.nodec.community_size) + "\t" + str(
                num_coms_before) + "\t" + str(num_coms_after) )

        print(network + "\t" + detection_algo + "\t" + str(budget) + "\t" + str(deception_before) + "\t" +
              str(deception_after) + "\t" +
              str(time) + str(self.nodec.target_community_id) + "\t" + str(self.nodec.community_size) + "\t" + str(
                num_coms_before) + "\t" + str(num_coms_after)  )

    def runAllExperiments(self):
        self.file_experiments.write(
            "Network" + "\t" + "DetectionAlgorithm" + "\t" + "BudgetUpdates" + "\t" + "InitialDeceptionValue" + "\t" +
            "FinalDeceptionValue" + "\t" +
            "Time(s)" + "\t" + "SizeTargeetCom" + "\t" + "#CommunitiesBefore" + "\t" + "#CommunitiesAfter" )
        print(
            "Network" + "\t" + "DetectionAlgorithm" + "\t" + "BudgetUpdates" + "\t" + "InitialDeceptionValue" + "\t" +
            "FinalDeceptionValue" + "\t" +
            "Time(s)" + "\t" + "SizeTargetCom" + "\t" + "#CommunitiesBefore" + "\t" + "#CommunitiesAfter" )

        for budget in self.budget_updates:
            for network in self.datasets:
                for com_algo in self.detection_algos:
                    for internalP,externalP in zip(self.sampleInternalPercentage, self.sampleExternalAddPercentage):
                        for i in range(0, self.rounds):
                            self.runExperiment(network, com_algo,budget,internalP,externalP,self.worst_case_scenario)
        self.file_experiments.close()


def main():
    dataset_path = "./datasets/"
    datasets = ["rete"]  # "terrorist2","terrorist1","blogcatalog" #"pubmed","webkb", "cora", "citeseer","facebook","netscience"
    detection_algorithms = ["leiden"]  #louv # "walk" "infomap" "labelp" "greedy" # not in the paper "btw","eig","opt"
    budget_updates = [2]
    internal_edge_sample_sizes = [0.5]
    external_edge_sample_sizes = [0.5]
    worst_case_scenario=True

    num_rounds = 1
    experiments = ExperimentsDirNODEC(dataset_path, datasets, detection_algorithms, budget_updates, internal_edge_sample_sizes,
                                   external_edge_sample_sizes, num_rounds,worst_case_scenario)

    experiments.runAllExperiments()

if __name__ == '__main__':
    main()