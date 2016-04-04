# Community Deception
Community Deception Algorithms and Datasets

This project implemenents algorithms for the novel Community Deception problem. Simply stated this problem is about how to hide a community (i.e., a subset of nodes in a network) from community detection algorithms.

We  have implemented two community deception algorithms. The first based on Modularity and the second one based on the novel notion of Community Safeness.

We provide instructions about how to replicate experiments reported in the paper under review at PKDD2016.

* Move all the network datasets in a subfolder named datasets
* Set the parameters (e.g., community detection algorithm, budget, deception algorithm) inside the communityDeceptionRunner.R file
* Lanuch the R script. The execution will output the results in a file and also display a subset of them on the terminal
