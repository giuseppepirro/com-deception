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
            #print("*****", beta, "*****")
            ##this gives the best value the one giving hghest DELTA_del
            delta_del,best_node_to_delete=self.nodec.computeDeltaNodeDeletion()
            #print("Best Del Delta=",delta_del, " , Best Node=", best_node_to_delete)
            
            #this gives a single value as we are adding a new node
            delta_add,inf_i_Cj,inf_i,Cj_index, external_community_idx=self.nodec.computeDeltaNodeAddition()
            
            
            #print("Best Add Delta=",delta_add, " , Internal node influence= ", inf_i_Cj, "Total node influence= ", inf_i)
            
            ##this gives a value for each community member
            delta_move,move_details= self.nodec.computeDeltaNodeMoving()
            #print(delta_move, move_details)
            
            #print("#################")
            #print("Best Mov Delta=",delta_move, " , Move_details= ", move_details)
            #print("#################")
            #exit()
            #deltas=[delta_add,delta_del,delta_move]
            #best_index=deltas.index(max(deltas))
            #print(deltas)
            #print(best_index)
            
            
            if float(delta_move)>0:
                
                self.performNodeMove(move_details)
                print("Moved node ", move_details[0])
            else:
                if delta_del > delta_add:
                    
                    self.performNodeDeletion(best_node_to_delete)
                    print("deleted one node (", best_node_to_delete,")")
                    
                else:
                    self.performNewNodeAddition(inf_i_Cj, inf_i, Cj_index, external_community_idx)
                    print("added one node")
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
        
        
        
        

    def performNewNodeAddition(self, inf_i_in_dest_com, total_inf_i, dest_com_index, ext_com_idx):
        
        new_node_id = self.nodec.graph.vcount()
        new_node_name=str(new_node_id)
        
        #check if new name is already reserved
        while True:
            if new_node_name not in self.nodec.graph.vs["name"]:
                break
            new_node_id += 1
            new_node_name=str(new_node_id)
        
        
        new_node_id = self.nodec.graph.vcount()
        
        self.nodec.graph.add_vertices(1)
        self.nodec.graph.vs[new_node_id]["name"]= new_node_name
        
        dest_com=self.nodec.communities[dest_com_index]
        
        
        #setting up new edges for the new node
        #dest_nodes_indexes=random.sample(range(0, len(dest_com)), len(dest_com))
        dest_nodes_indexes = range(len(dest_com))
        
        edge_weight = inf_i_in_dest_com/len(dest_com)
        
        for i in dest_nodes_indexes:
            self.nodec.graph.add_edges([(new_node_name, dest_com[i])])
            self.nodec.graph.es[self.nodec.graph.get_eid(new_node_name, dest_com[i])]["weight"] = edge_weight
        
        
        #setting up inter-community edges adjacent to new node
        external_community = self.nodec.communities[ext_com_idx]
        external_dest_nodes=random.sample(external_community, len(external_community))
        inter_com_inf_budget = total_inf_i-inf_i_in_dest_com
        
        for i in external_dest_nodes:
            if inter_com_inf_budget >= edge_weight:
                self.nodec.graph.add_edges([(new_node_name, i)])
                self.nodec.graph.es[self.nodec.graph.get_eid(new_node_name, i)]["weight"] = edge_weight
                inter_com_inf_budget -= edge_weight
            
            else: 
                break
            
            
            
        #self.nodec.target_community.append(new_node_id)
        
        self.nodec.communities[dest_com_index].append(new_node_name)
        
        
        self.nodec.degi_Ci = self.nodec.getTargetNodeDegreePerCommunity()
        
        self.nodec.degi_Ci_outInf = self.nodec.getTargetNodeOutInfluencePerCommunity()
        self.nodec.degi_Ci_inInf = self.nodec.getTargetNodeInInfluencePerCommunity()
        self.nodec.degi_Ci_inf = self.nodec.degi_Ci_outInf+self.nodec.degi_Ci_inInf
        
        
        
        coms_degs_copy = self.nodec.community_degrees
        
        
        coms_degs_copy[dest_com_index] = coms_degs_copy[dest_com_index] + inf_i_in_dest_com
        self.nodec.community_degrees = coms_degs_copy
        
        
        internal_edges_copy = self.nodec.E_Ci
        internal_edges_copy[dest_com_index] = internal_edges_copy[dest_com_index] + inf_i_in_dest_com
        self.nodec.E_Ci = internal_edges_copy
        
        
        return self.nodec.graph



    def performNodeMove(self, move_details):
        node_to_move= move_details[0] #node to move from Ci in Cj
        new_community=int(move_details[2]) #  Cj
        old_community=self.nodec.communities[self.nodec.getNodeCommunity(node_to_move)] #Ci
        new_inf_new_com=float(move_details[3]) #total edges in Cj
        
        #the old community is the target community
        #node_index_in_old_com=self.nodec.target_community.index(node_to_move)

        destination_com=self.nodec.communities[new_community]
        #destination_nodes_indexes = range(0, len(destination_com))


        
        ## EDGE DELETIONS in Ci
        out_neighs = [self.nodec.graph.vs[i]["name"] for i in self.nodec.graph.neighbors(node_to_move, mode="out") if i in old_community]
        in_neighs = [self.nodec.graph.vs[i]["name"] for i in self.nodec.graph.neighbors(node_to_move, mode="in") if i in old_community]
        
        
        
        for n in range(1, len(out_neighs)): #keep the first out-edge to preserve communication
            self.nodec.graph.delete_edges( [(node_to_move, out_neighs[n])] )
        
        for n in range(1, len(in_neighs)):  #keep the first in-edge to preserve communication
            self.nodec.graph.delete_edges( [(in_neighs[n], node_to_move)] )
        
        node_out_neighs = [self.nodec.graph.vs[i]["name"] for i in self.nodec.graph.neighbors(node_to_move, mode="out")]
        
        new_neighbors_in_new_com = [i for i in destination_com if i not in node_out_neighs]
        print("node to move = ",node_to_move)
        #print(node_out_neighs)
        #print(destination_com)
        #print(new_neighbors_in_new_com)
        edge_weight = new_inf_new_com/len(new_neighbors_in_new_com)
        
        ## EDGE ADDITIONS in Cj
        for i in new_neighbors_in_new_com:
            self.nodec.graph.add_edges([(node_to_move, i)])
            self.nodec.graph.es[self.nodec.graph.get_eid(node_to_move, i)]["weight"] = edge_weight
     
        
        degree_in_com_copy = self.nodec.degi_Ci
        
        self.nodec.degi_Ci = degree_in_com_copy
        
        self.nodec.target_community.remove(node_to_move)

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
        #print(self.nodec.graph)
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
        print("Final ComH=", self.nodec.target_community)
        print("Final modularity=", self.nodec.getModularity(g,self.nodec.communities_after_object))
        print("Final number of nodes and edges=", new_graph.ecount(),new_graph.vcount())
        #print(self.nodec.graph)
        
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
    datasets = ["freeman"]  # "terrorist2","terrorist1","blogcatalog" #"pubmed","webkb", "cora", "citeseer","facebook","netscience"
    detection_algorithms = ["leiden"]  #louv # "walk" "infomap" "labelp" "greedy" # not in the paper "btw","eig","opt"
    budget_updates = [10]
    internal_edge_sample_sizes = [0.5]
    external_edge_sample_sizes = [0.5]
    worst_case_scenario=True

    num_rounds = 1
    experiments = ExperimentsDirNODEC(dataset_path, datasets, detection_algorithms, budget_updates, internal_edge_sample_sizes,
                                   external_edge_sample_sizes, num_rounds,worst_case_scenario)

    experiments.runAllExperiments()

if __name__ == '__main__':
    main()