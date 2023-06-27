
from InflDec import *
from datetime import datetime
import timeit
from cdlib import evaluation

class ExperimentsInflDec:
    def __init__(self, dataset_path, datasets, detection_algos, budget_percentages, 
                 sampleInternalPercentage, sampleExternalAddPercentage, rounds):
        """This method initializes ExperimentsInflDec attributes,
        which are essecially the parameters used for all experiments."""
        self.dataset_path = dataset_path
        self.datasets = datasets
        self.detection_algos = detection_algos
        self.budget_updates = budget_percentages
        self.sampleInternalPercentage = sampleInternalPercentage
        self.sampleExternalAddPercentage = sampleExternalAddPercentage
        self.rounds = rounds
        self.infldec = InflDec(dataset_path=self.dataset_path)

        exp_file_name = self.dataset_path +"results/infldec/"+ "InflDec_"+datetime.today().strftime('%Y-%m-%d-%H:%M:%S') + ".experiments.txt"
        self.file_experiments = open(exp_file_name, "w+")
        self.current_algo = None


    
    def com_deception_infldec(self, budget, sampleInternalPercentage, sampleExternalAddPercentage):
        add_gain = 0
        del_gain = 0
        beta = budget

        while (True):
            #get best OUT-GOING inter-community addition candidate
            best_source_add_out, best_target_add_out, add_gain_out, best_target_com_out = self.infldec.getBestOutAdditionInflDec(
                sampleInternalPercentage, sampleExternalAddPercentage)

            #get best IN-GOING inter-community addition candidate
            best_source_add_in, best_target_add_in, add_gain_in, best_target_com_in = self.infldec.getBestInAdditionInflDec(
                sampleInternalPercentage, sampleExternalAddPercentage)


            #choose best among addition candidates
            is_out_inter_add = True
            if add_gain_out >= add_gain_in:
                best_source_add = best_source_add_out 
                best_target_add = best_target_add_out
                best_target_com = best_target_com_out
                add_gain = add_gain_out
            else: 
                is_out_inter_add = False
                best_source_add = best_source_add_in
                best_target_add = best_target_add_in
                best_target_com = best_target_com_in
                add_gain = add_gain_in
           
            
            # add
            #best_source_add, best_target_add, best_target_com = self.diffuser.getBestInterAdditionDiffusion()
            
            # del
            #gets best intra-community edge candidate for deletion
            best_source_del, best_target_del, del_gain = self.infldec.getBestIntraDeletionInflDec()

            if add_gain >= del_gain and add_gain > 0:
                if is_out_inter_add:
                    #print("addition_out")
                    edge_weight = self.infldec.getNodeInfluence(best_source_add, best_target_add)
                    weight = {"weight": edge_weight}
                    self.infldec.graph.add_edges([(best_source_add, best_target_add)], attributes=weight)
                    self.infldec.external_weights[self.infldec.target_community.index(best_source_add)] = \
                    self.infldec.external_weights[
                        self.infldec.target_community.index(
                            best_source_add)] + edge_weight
                    self.infldec.total_influence[self.infldec.target_community.index(best_source_add)] = \
                    self.infldec.total_influence[
                        self.infldec.target_community.index(
                            best_source_add)] + edge_weight
                    self.infldec.external_influence_per_community[
                        self.infldec.target_community.index(best_source_add)][best_target_com] = \
                        self.infldec.external_influence_per_community[
                            self.infldec.target_community.index(best_source_add)][best_target_com] + edge_weight
                else:
                    #print("addition_in")
                    edge_weight = self.infldec.getNodeInfluence(best_source_add, best_target_add)
                    weight = {"weight": edge_weight}
                    self.infldec.graph.add_edges([(best_source_add, best_target_add)], attributes=weight)
                    self.infldec.external_weights[self.infldec.target_community.index(best_target_add)] = \
                    self.infldec.external_weights[
                        self.infldec.target_community.index(
                            best_target_add)] + edge_weight
                    self.infldec.total_influence[self.infldec.target_community.index(best_target_add)] = \
                    self.infldec.total_influence[
                        self.infldec.target_community.index(
                            best_target_add)] + edge_weight
                    self.infldec.external_influence_per_community[
                        self.infldec.target_community.index(best_target_add)][best_target_com] = \
                        self.infldec.external_influence_per_community[
                            self.infldec.target_community.index(best_target_add)][best_target_com] + edge_weight
                    
            # return 1
            elif del_gain > 0:
                print("intra_del")
                del_weight = self.infldec.graph.es[(best_source_del, best_target_del)]["weight"][0]
                if self.infldec.graph.are_connected(best_source_del, best_target_del):
                    self.infldec.graph.delete_edges([(best_source_del, best_target_del)])
                    # Update the weights of the involved nodes
                    self.infldec.internal_weights[self.infldec.target_community.index(best_source_del)] = \
                    self.infldec.internal_weights[
                        self.infldec.target_community.index(
                            best_source_del)] - del_weight
                    self.infldec.total_influence[self.infldec.target_community.index(best_source_del)] = \
                        self.infldec.total_influence[
                        self.infldec.target_community.index(best_source_del)] - del_weight
                    try:
                        self.infldec.internal_weights[self.infldec.target_community.index(best_target_del)] = \
                        self.infldec.internal_weights[
                            self.infldec.target_community.index(
                                best_target_del)] - del_weight
                    except ValueError as e:
                        self.infldec.internal_weights[self.infldec.target_community.index(best_source_del)] = \
                            self.infldec.internal_weights[
                                self.infldec.target_community.index(
                                    best_source_del)] - del_weight
                    try:
                        self.infldec.total_influence[self.infldec.target_community.index(best_target_del)] = \
                        self.infldec.total_influence[
                            self.infldec.target_community.index(
                                best_target_del)] - del_weight
                    except ValueError as e:
                        self.infldec.total_influence[self.infldec.target_community.index(best_source_del)] = \
                            self.infldec.total_influence[
                                self.infldec.target_community.index(
                                    best_source_del)] - del_weight   
            beta = beta - 1
            if (beta > 0 and (add_gain > 0 or del_gain > 0)):
                continue
            else:
                break

        return self.infldec.graph


    def runExperiment(self, network, detection_algo,budget, internalP,externalP):
        self.infldec.read_weighted_graph(network)
        recompute_communities = True
        self.infldec.initializeCommunityData(detection_algo=detection_algo, dataset=network,
                                             recompute_communities=recompute_communities)

        print("Number edges before=",self.infldec.graph.ecount())
        coms_before = self.infldec.communities_object
        #print(coms_before.communities)
        num_coms_before = len(coms_before.communities)
        # print("Number of communities before=", (num_coms_before))
        deception_before, member_for_community = self.infldec.getDeceptionScore(coms_before.communities)
        ##### BEGIN InflDec HERE
        start = timeit.default_timer()
        #
        sampleInternalPercentage=internalP
        sampleExternalAddPercentage=externalP
        updated_graph=self.com_deception_infldec(budget, sampleInternalPercentage,
                                             sampleExternalAddPercentage)

        #
        # communities in the updated graph
        stop = timeit.default_timer()
        time=stop - start
        ##### END InflDec HERE

        coms_after = self.infldec.computeCommunitiesAfterUpdateCDLIB(detection_algo, updated_graph)
        coms_after = self.infldec.communities_after_object
        num_coms_after = len(coms_after.communities)

        print(coms_after.communities)
        # print("Number of communities after=", (num_coms_after))
        deception_after, member_for_community = self.infldec.getDeceptionScore(coms_after.communities)
        member_for_community = np.count_nonzero(np.array(member_for_community) > 0)
        nmi=evaluation.normalized_mutual_information(coms_before, coms_after)[0]

        self.file_experiments.write(
            "\n" + network + "\t" + detection_algo + "\t" + str(budget) + "\t" + str(deception_before) + "\t" +
            str(deception_after) + "\t" +
            str(time) + str(self.infldec.target_community_id) + "\t" + str(self.infldec.target_community_size) + "\t" + str(
                num_coms_before) + "\t" + str(num_coms_after) + "\t" + str(nmi) + "\t" + str(
                self.infldec.target_community_id) + "\t" + str(member_for_community))

        print(network + "\t" + detection_algo + "\t" + str(budget) + "\t" + str(deception_before) + "\t" +
              str(deception_after) + "\t" +
              str(time) + str(self.infldec.target_community_id) + "\t" + str(
            self.infldec.target_community_size) + "\t" + str(
            num_coms_before) + "\t" + str(num_coms_after) + "\t" + str(nmi) + "\t" + str(
            self.infldec.target_community_id) + "\t" + str(member_for_community))


    def runAllExperiments(self):
        self.file_experiments.write(
            "Network" + "\t" + "DetectionAlgorithm" + "\t" + "BudgetUpdates" + "\t" + "InitialDeceptionValue" + "\t" +
            "FinalDeceptionValue" + "\t" +
            "Time(s)" + "\t" + "SizeTargetCom" + "\t" + "#CommunitiesBefore" + "\t" + "#CommunitiesAfter" + "\t" + "NMI" + "\t" + "TargetComID" + "\t" + "CommunitySplit")
        print(
            "Network" + "\t" + "DetectionAlgorithm" + "\t" + "BudgetUpdates" + "\t" + "InitialDeceptionValue" + "\t" +
            "FinalDeceptionValue" + "\t" +
            "Time(s)" + "\t" + "SizeTargetCom" + "\t" + "#CommunitiesBefore" + "\t" + "#CommunitiesAfter" + "\t" + "NMI" + "\t" + "TargetComID" + "\t" + "CommunitySplit")


        for budget in self.budget_updates:
                    for network in self.datasets:
                        for com_algo in self.detection_algos:
                            for internalP,externalP in zip(self.sampleInternalPercentage, self.sampleExternalAddPercentage):
                                for i in range(0, self.rounds):
                                    self.runExperiment(network, com_algo, budget, internalP,externalP)
        self.file_experiments.close()


def main():
    dataset_path = "./datasets/"
    datasets = ["wiki"]  # 
    detection_algorithms = ["leiden"]  
    budget_updates = [20]
    internal_edge_sample_sizes = [0.5]
    external_edge_sample_sizes = [0.5]
    num_rounds = 1
    
    #creating and initializing ExperimentsInlfDec object
    experiments = ExperimentsInflDec(dataset_path, datasets, detection_algorithms,
                                     budget_updates, internal_edge_sample_sizes,
                                     external_edge_sample_sizes, num_rounds)
    #call method to runAllExperiments with all selected parameter variations
    experiments.runAllExperiments()

if __name__ == '__main__':
    main()
