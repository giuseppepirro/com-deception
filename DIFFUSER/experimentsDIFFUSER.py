
from DIFFUSER import *
from datetime import datetime
import timeit
class ExperimentsDIFFUSER:
    def __init__(self, dataset_path, datasets, detection_algos, budget_percentages, sampleInternalPercentage,
                 sampleExternalAddPercentage, rounds):
        self.dataset_path = dataset_path
        self.datasets = datasets
        self.detection_algos = detection_algos
        self.budget_updates = budget_percentages
        self.sampleInternalPercentage = sampleInternalPercentage
        self.sampleExternalAddPercentage = sampleExternalAddPercentage
        self.rounds = rounds
        self.diffuser = DIFFUSER(dataset_path=self.dataset_path)
        exp_file_name = self.dataset_path +"results/diffuser-weighted/"+ "diffuser_weighted_"+datetime.today().strftime('%Y-%m-%d-%H:%M:%S') + ".experiments"
        self.file_experiments = open(exp_file_name, "w+")
        self.current_algo = None

    def com_decept_diffusion(self, budget, sampleInternalPercentage, sampleExternalAddPercentage):
        add_gain = 0
        del_gain = 0
        beta = budget
        while (True):
            # add
            best_source_add, best_target_add, add_gain, best_target_com = self.diffuser.getBestAdditionDiffusion(
                sampleInternalPercentage, sampleExternalAddPercentage)

            # del
            best_inter_source_del, best_target_inter_del, inter_del_gain, edge_direction = self.diffuser.getBestExtDeletionDiffusion(
              sampleInternalPercentage)
            #inter_del_gain=0

            best_intra_source_del, best_intra_target_del, intra_del_gain = self.diffuser.getBestIntDeletionDiffusion()
            if inter_del_gain >= intra_del_gain and inter_del_gain > 0:
                del_gain = inter_del_gain
            else:
                best_source_del = best_intra_source_del
                best_target_del = best_intra_target_del
                del_gain = intra_del_gain

            if add_gain >= del_gain and add_gain > 0:
                edge_weight = self.diffuser.getNodeSimilarity(best_source_add, best_target_add)
                weight = {"weight": edge_weight}
                self.diffuser.graph.add_edges([(best_source_add, best_target_add)], attributes=weight)
                self.diffuser.external_weights[self.diffuser.target_community.index(best_source_add)] = \
                self.diffuser.external_weights[
                    self.diffuser.target_community.index(
                        best_source_add)] + edge_weight
                self.diffuser.total_weights[self.diffuser.target_community.index(best_source_add)] = \
                self.diffuser.total_weights[
                    self.diffuser.target_community.index(
                        best_source_add)] + edge_weight

                self. diffuser.external_weights_per_community[
                    self.diffuser.target_community.index(best_source_add)][best_target_com] = \
                    self.diffuser.external_weights_per_community[
                        self.diffuser.target_community.index(best_source_add)][best_target_com] + edge_weight
            # return 1
            elif del_gain > 0:
                del_weight = self.diffuser.graph.es[(best_source_del, best_target_del)]["weight"][0]
               #
                if self.diffuser.graph.are_connected(best_source_del,best_target_del):
                    self.diffuser.graph.delete_edges([(best_source_del, best_target_del)])
                    # Update the weights of the involved nodes
                    self.diffuser.internal_weights[self.diffuser.target_community.index(best_source_del)] = \
                    self.diffuser.internal_weights[
                        self.diffuser.target_community.index(
                            best_source_del)] - del_weight
                    self.diffuser.total_weights[self.diffuser.target_community.index(best_source_del)] = \
                        self.diffuser.total_weights[
                        self.diffuser.target_community.index(best_source_del)] - del_weight
                    try:
                        self.diffuser.internal_weights[self.diffuser.target_community.index(best_target_del)] = \
                        self.diffuser.internal_weights[
                            self.diffuser.target_community.index(
                                best_target_del)] - del_weight
                    except ValueError as e:
                        self.diffuser.internal_weights[self.diffuser.target_community.index(best_source_del)] = \
                            self.diffuser.internal_weights[
                                self.diffuser.target_community.index(
                                    best_source_del)] - del_weight
                    try:
                        self.diffuser.total_weights[self.diffuser.target_community.index(best_target_del)] = \
                        self.diffuser.total_weights[
                            self.diffuser.target_community.index(
                                best_target_del)] - del_weight
                    except ValueError as e:
                        self.diffuser.total_weights[self.diffuser.target_community.index(best_source_del)] = \
                            self.diffuser.total_weights[
                                self.diffuser.target_community.index(
                                    best_source_del)] - del_weight

            beta = beta - 1
            if (beta > 0 and (add_gain > 0 or del_gain > 0)):
                continue
            else:
                break
        return self.diffuser.graph


    def runExperiment(self, network, detection_algo,budget,internalnalP,externalP):
        self.diffuser.read_weighted_graph(network)
        recompute_communities = True
        self.diffuser.initializeCommunityData(detection_algo=detection_algo, dataset=network,
                                              recompute_communities=recompute_communities)
        coms_before = self.diffuser.communities
        num_coms_before = len(coms_before)

        deception_before = self.diffuser.getDeceptionScore()

        ##### BEGIN DIFFUSER HERE
        start = timeit.default_timer()

        sampleInternalPercentage=internalnalP
        sampleExternalAddPercentage=externalP
        new_graph=self.com_decept_diffusion(budget, sampleInternalPercentage,
                                            sampleExternalAddPercentage)
        # communities in the updated graph
        g = new_graph
        stop = timeit.default_timer()
        time=stop - start
        ##### END DIFFUSER HERE

        coms_after = self.diffuser.computeCommunitiesAfterUpdate(detection_algo, g)
        num_coms_after = len(coms_after)
        deception_after = self.diffuser.getDeceptionScore()
        nmi = compare_communities(coms_before, coms_after, method="nmi")

        self.file_experiments.write(
           "\n" + network + "\t" + detection_algo + "\t" + str(budget) + "\t" + str(deception_before) + "\t" +
           str(deception_after) + "\t" +
           str(time) + str(self.diffuser.target_community_id) + "\t" + str(self.diffuser.community_size) + "\t" + str(
                num_coms_before) + "\t" + str(num_coms_after) + "\t" + str(nmi))

        print(network + "\t" + detection_algo + "\t" + str(budget) + "\t" + str(deception_before) + "\t" +
              str(deception_after) + "\t" +
              str(time) + str(self.diffuser.target_community_id) + "\t" + str(self.diffuser.community_size) + "\t" + str(
                num_coms_before) + "\t" + str(num_coms_after) + "\t" + str(nmi))

    def runAllExperiments(self):
        self.file_experiments.write(
            "Network" + "\t" + "DetectionAlgorithm" + "\t" + "BudgetUpdates" + "\t" + "InitialDeceptionValue" + "\t" +
            "FinalDeceptionValue" + "\t" +
            "Time(s)" + "\t" + "SizeTargeetCom" + "\t" + "#CommunitiesBefore" + "\t" + "#CommunitiesAfter" + "\t" + "NMI")
        print(
            "Network" + "\t" + "DetectionAlgorithm" + "\t" + "BudgetUpdates" + "\t" + "InitialDeceptionValue" + "\t" +
            "FinalDeceptionValue" + "\t" +
            "Time(s)" + "\t" + "SizeTargetCom" + "\t" + "#CommunitiesBefore" + "\t" + "#CommunitiesAfter" + "\t" + "NMI")

        for budget in self.budget_updates:
            for network in self.datasets:
                for com_algo in self.detection_algos:
                    for internalP,externalP in zip(self.sampleInternalPercentage, self.sampleExternalAddPercentage):
                        for i in range(0, self.rounds):
                            self.runExperiment(network, com_algo, budget,internalP,externalP)
        self.file_experiments.close()


def main():
    dataset_path = "./datasets/"
    datasets = ["terrorist2"]  # "terrorist2","terrorist1","blogcatalog" #"pubmed","webkb", "cora", "citeseer","facebook"
    detection_algorithms = ["louv"]  #louv # "walk" "infomap" "labelp" "greedy" # not in the paper "btw","eig","opt"
    budget_updates = [1,2,3,4,5,6,7,8,9,10]
    internal_edge_sample_sizes = [0.5]
    external_edge_sample_sizes = [0.5]

    num_rounds = 1
    experiments = ExperimentsDIFFUSER(dataset_path, datasets, detection_algorithms, budget_updates, internal_edge_sample_sizes,
                                      external_edge_sample_sizes, num_rounds)
    experiments.runAllExperiments()

if __name__ == '__main__':
    main()
