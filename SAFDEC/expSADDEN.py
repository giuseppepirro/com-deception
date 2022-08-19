from Utilities import *
from datetime import datetime
import SADDEN
import timeit
class ExperimentsSADDEN:
    def __init__(self, dataset_path, datasets, detection_algos, budget_percentages,rounds):
        self.dataset_path = dataset_path
        self.datasets = datasets
        self.detection_algos = detection_algos
        self.budget_updates = budget_percentages
        self.rounds = rounds
        self.runner = Utils(dataset_path=self.dataset_path)
        exp_file_name = self.dataset_path +"results/sadden/"+ "sadden_"+datetime.today().strftime('%Y-%m-%d-%H:%M:%S') + ".experiments"
        self.file_experiments = open(exp_file_name, "w+")
        self.current_algo = None

    def runExperiment(self, network, detection_algo,budget):
        self.runner.read_network(network)
        recompute_communities = True
        self.runner.initializeCommunityData(detection_algo=detection_algo,
                                            recompute_communities=recompute_communities)
        coms_before = self.runner.communities
        num_coms_before = len(coms_before)
        print("Number of communities before=", (num_coms_before))
        deception_before = self.runner.getDeceptionScore()
        print("Deception score before=", deception_before)

        ##### BEGIN SADDEN HERE
        start = timeit.default_timer()
        IG_edgeList_=SADDEN.runSADDENSingleCom(self.runner.graph,self.runner.target_community,self.runner.communities,budget)
        # communities in the updated graph
        g = igraph.Graph(directed=False)
        g.add_vertices(self.runner.graph.vcount())
        g.add_edges(IG_edgeList_)
        stop = timeit.default_timer()
        time=stop - start
        ##### END SADDEN HERE

        coms_after = self.runner.computeCommunitiesAfterUpdate(detection_algo,g)
        num_coms_after = len(coms_after)
        print("Number of communities after=", (num_coms_after))
        deception_after = self.runner.getDeceptionScore()
        print("Deception score after=", deception_after)
        nmi = compare_communities(coms_before, coms_after, method="nmi")

        self.file_experiments.write(
            "\n" + network + "\t" + detection_algo + "\t" + str(budget) + "\t" + str(deception_before) + "\t" +
            str(deception_after) + "\t" +
            str(time) + str(self.runner.target_community_id) + "\t" + str(self.runner.community_size) + "\t" + str(
                num_coms_before) + "\t" + str(num_coms_after) + "\t" + str(nmi))

        print(network + "\t" + detection_algo + "\t" + str(budget) + "\t" + str(deception_before) + "\t" +
              str(deception_after) + "\t" +
              str(time) + str(self.runner.target_community_id) + "\t" + str(
            self.runner.community_size) + "\t" + str(
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
                        for i in range(0, self.rounds):
                             self.runExperiment(network, com_algo, budget)
                self.file_experiments.close()


def main():
    dataset_path = "../datasets/edgelist/"
    datasets = ["football"]  # ,"terrorist2","blogcatalog" #"geometry" #"pubmed","webkb", "cora", "citeseer","facebook"
    com_algos = ["louv"]  # btw very slow #,"infomap","walk","greedy","spin","labelp","btw","eig","opt"
    budget_updates = [2,3,5,10]
    num_rounds = 1

    experiments = ExperimentsSADDEN(dataset_path, datasets, com_algos, budget_updates, num_rounds)
    experiments.runAllExperiments()

if __name__ == '__main__':
    main()
