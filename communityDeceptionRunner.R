library(igraph)
source("utils.R")
source("algoDeceptionMod.R")
source("algoDeceptionSaf.R")
##### ALGO ACRONYMS ###############
# louv: cluster_louvain
# infomap: cluster_infomap
# greedy: cluster_fast_greedy
# walk: cluster_walktrap
# spin: cluster_spinglass
# labelp: cluster_label_prop
# opt: cluster_optimal
# betw: cluster_edge_betweenness
# eig: cluster_leading_eigen
##################################
##### DATASETS ACRONYMS ###############
##### SYNTHETIC DATASETS
# Barabasi: Synthetic network generated according to the Barabasi-Albert model
# Erdos: Synthetic network generated according to the Erdos-Ranyi model
# Watts: Synthetic network generated according to the Watts-Strogatz model
##### REAL-WORLD DATASETS
# karate: Zachary Karate's Club - A karate club at a US university, as described by Wayne Zachary
# lesmis : Les Miserables - Coappearances of characters in Victor Hugo's novel "Les Miserables" compiled by Valdis Krebs
# dolphins: Associations between 62 dolphins in a community living off Doubtful Sound,New Zealand, as compiled by Lusseau et al. (2003)
# football: American football games between Division IA colleges during regular season Fall 2000, as compiled by M. Girvan and M. Newman
# madrid: Madrid Terrorist Network compiled by Jose A. Rodriguez of the University of Barcelona
# polbooks: Books about US Politics
# power: Representation of the topology of the Western States Power Grid of the United States, compiled by Duncan Watts and Steven Strogatz# jazz com troppo piccole: List of edges of the network of Jazz musicians compiled by compiled by A. Arenas and colleagues
#generated: Networks Generated with the Lancichinetti and Fortunato's benchmark. Rename the file in the dataset folder to be used to generated.edgelist
##################################

current_algo<-"louv"
dataset<-"karate"

##Community Deception Algorithm
#algo_deception<-"modularity"
algo_deception<-"both-algos"

#chose whether to show the communities
make_plot<-FALSE
#useful to print the header only once
done=FALSE

#########################################
#Budget in terms of edge changes
max_budget<-4

##Type of community: worst-case (it belongs to the output of the community detection algorithm)
types_of_experiment<-"worst_case"
#types_of_experiment<-"NO_worst_case"
## Random nodes are generated; Small=5% of total nodes; Medium=10% of total nodes; Big=15% of total nodes
size_community_NO_worst_case<-"Small"
# 

 ########################
initialNetwork <- initializeNetwork(dataset)
initialNetwork<-simplify(initialNetwork,TRUE,TRUE)
  
initialClustering <-computeCommunities(initialNetwork,current_algo)
#############
  initialM<-membership(initialClustering)
  total_num_communities<-length(initialClustering)
  total_num_nodes<-vcount(initialNetwork)
  total_num_edges<-ecount(initialNetwork)
  
  network_to_work_modularity <- initializeNetwork(dataset)
  network_to_work_modularity<-simplify(network_to_work_modularity,TRUE,TRUE)
  
  network_to_work_safeness <- initializeNetwork(dataset)
  network_to_work_safeness<-simplify(network_to_work_safeness,TRUE,TRUE)
  
  
  
  ##THIS PART CHOOSES a random \ComH among the existing one (according to the settings above)
  if(identical(types_of_experiment,"worst_case"))
  {
    community<-getCommunityToProtectWorstCase(initialM,total_num_communities)
    
  }
  if(identical(types_of_experiment,"NO_worst_case"))
  {
    community<- getCommunityToProtect(initialNetwork,total_num_nodes,size_community_NO_worst_case)
  }
  print("===================Selected Community to deceive===================")
  print(community)
  print("===================================================================")
  ## modularity and deception values for the initial network
  
  detectionTime<-system.time({
    wc <-computeCommunities(network_to_work_modularity,current_algo) 
    myM<-membership(wc)
  })
  
  #deception score
  initialdeceptionScore<-0
  initialdeceptionScore<-getdeceptionScore(myM,network_to_work_modularity,community)
  print(paste("initial_deception_score=",initialdeceptionScore))
  
  #modularity
  initial_modularity<-getModularity(myM,network_to_work_modularity)
  print(paste("initial_modularity_value=",initial_modularity))
  
  #safeness  
  initial_safeness<-getSafenessScore(myM,network_to_work_modularity,community)
  print(paste("initial_safeness_value=",initial_safeness))
  
  
  ################## INITIAL VALUES
#########################################

##
for(budget_e in 1:max_budget)
{
  budget<-budget_e
  

  file_name<-paste(current_algo,algo_deception,budget,dataset,types_of_experiment,sep="_")
  file_name<-paste("./",file_name,sep = "")
  
  file_name_for_EXCEL<-paste(current_algo,algo_deception,budget,dataset,types_of_experiment,"TAB-SEPARATED",sep="_")
  file_name_for_EXCEL<-paste("./",file_name_for_EXCEL,sep = "")
  
  ###HEADERS
  if(!done)
  {
  write(paste("NETWORK","DETECTION-ALGO","BUDGET","INITAL-MODULARITY",
              "FINAL-MODULARITY","INITAL-SAFENESS","FINAL-SAFENESS","INITAL-DECEPTION",
              "FINAL-DECEPTION-MODULARITY","FINAL-DECEPTION-SAFENESS","TIME-MODULARITY",
              "TIME-SAFENESS","TIME-DETECTION",sep="\t"),file=paste(file_name_for_EXCEL,".dat"),append=TRUE)
  done=TRUE
    }
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 
#================ BEGIN NETWORK CHANGES ====================================================
  ######################## ############ ############ ############  BEGIN changes modularity 
  modTime<-system.time({
  
        upd_modularity<-getBestUpdateModularity(myM,network_to_work_modularity,community)
    
    if(is.null(upd_modularity)){
      print(paste("NO UPDATES FOR MODULARITY ARE POSSIBLE"))
      break
    }
 
  else if(upd_modularity[3]==666)
  {
    write(paste("INTER_ADD_MODULARITY=","(",upd_modularity[1],",",upd_modularity[2],")"),file=paste(file_name,".dat"),append=TRUE)
    print(paste("MODULARITY: INTER_ADD=","(",upd_modularity[1],",",upd_modularity[2],")"))
    
    network_to_work_modularity<- add_edges(network_to_work_modularity,c(upd_modularity[1],upd_modularity[2]))
     
    }
  else  if(upd_modularity[3]==999) 
    {
    write(paste("INTRA_DEL_MODULARITY=","(",upd_modularity[1],",",upd_modularity[2],")"),file=paste(file_name,".dat"),append=TRUE)
    print(paste("MODULARITY: INTRA_DEL=","(",upd_modularity[1],",",upd_modularity[2],")"))
        network_to_work_modularity<- delete_edges(network_to_work_modularity,paste(upd_modularity[1],"|",upd_modularity[2]))
        }
  })
  

  safTime<-system.time({
  upd_safeness<-getBestUpdateSafeness(myM,network_to_work_safeness,community)

  if(is.null(upd_safeness)){
    print(paste("NO UPDATES FOR SAFENESS ARE POSSIBLE"))
    break
    
  }
  if(upd_safeness[3]==666)
  {
    write(paste("INTER_ADD_SAFENESS=","(",upd_safeness[1],",",upd_safeness[2],")"),file=paste(file_name,".dat"),append=TRUE)
    print(paste("SAFENESS: INTER_ADD=","(",upd_safeness[1],",",upd_safeness[2],")"))
    network_to_work_safeness<- add_edges(network_to_work_safeness,c(upd_safeness[1],upd_safeness[2]))
  }
  else if(upd_safeness[3]==999) {
    write(paste("INTRA_DEL_SAFENESS=","(",upd_safeness[1],",",upd_safeness[2],")"),file=paste(file_name,".dat"),append=TRUE)
    print(paste("SAFENESS: INTRA_DEL=","(",upd_safeness[1],",",upd_safeness[2],")"))
    network_to_work_safeness<- delete_edges(network_to_work_safeness,paste(upd_safeness[1],"|",upd_safeness[2]))
  }
  })

###### Compute again modularity, safeness and deception after applying all changes on the two copies of the network
##MODULARITY
  wc_modularity <-computeCommunities(network_to_work_modularity,current_algo) 
  myM_modularity<-membership(wc_modularity)
  final_deception_score_modularity<-getdeceptionScore(myM_modularity,network_to_work_modularity,community)
  print(paste("final_deception_score_modularity=",final_deception_score_modularity))
  final_modularity<-getModularity(myM_modularity,network_to_work_modularity)
  print(paste("final_modularity_value=",final_modularity))  
  
## SAFENESS:
  wc_safeness <-computeCommunities(network_to_work_safeness,current_algo) 
  myM_safeness<-membership(wc_safeness)
  final_deception_score_safeness<-getdeceptionScore(myM_safeness,network_to_work_safeness,community)
  print(paste("final_deception_score_safeness=",final_deception_score_safeness))
  final_safeness<-getSafenessScore(myM_safeness,network_to_work_safeness,community)
  print(paste("final_safeness_value_safeness=",final_safeness))

  ###TIME
  print(paste("TIME_MODULARITY=",modTime[3]))
  print(paste("TIME_SAFENESS=",safTime[3]))

################# Creating a Plot with the results
  
  #modularity
  mscore<-round(final_deception_score_modularity,4)
  num_communities_mod<-length(wc_modularity)
  
  #safeness
  sscore<-round(final_deception_score_safeness,4)
  num_communities_safeness<-length(wc_safeness)
  
  ###string representation of the community chosen
  com_to_protect_string<-capture.output(cat(community))
  
  
if(make_plot)
{
attach(mtcars)
par(mfrow=c(1,3))
ipcore<-round(initialdeceptionScore,4)

plot(initialClustering, initialNetwork,
     main=paste("Det=",current_algo, "scoreI=",ipcore, "com=", com_to_protect_string))#vertex.label=NA
plot(wc_modularity,network_to_work_modularity, main=paste(" Prot=Modularity", "scoreF=",mscore))#,vertex.label=NA
plot(wc_safeness,network_to_work_safeness, main=paste(" Prot=Safeness", "scoreF=",sscore))#,vertex.label=NA

}
 

### Write on file the results
write(paste("DeTection_Algorithm=",current_algo),file=paste(file_name,".dat"),append=TRUE)
write(paste("DeCeption_Algorithm=",algo_deception),file=paste(file_name,".dat"),append=TRUE)
write(paste("Network=",dataset),file=paste(file_name,".dat"),append=TRUE)
write(paste("#Initial_Communities=",total_num_communities),file=paste(file_name,".dat"),append=TRUE)
write(paste("#Communities_after_MODULARITY=",num_communities_mod),file=paste(file_name,".dat"),append=TRUE)
write(paste("#Communities_after_SAFENESS=",num_communities_safeness),file=paste(file_name,".dat"),append=TRUE)
write(paste("#Nodes=",total_num_nodes),file=paste(file_name,".dat"),append=TRUE)
write(paste("#Edges=",total_num_edges),file=paste(file_name,".dat"),append=TRUE)

write(paste("Community_to_Deceive=",com_to_protect_string),file=paste(file_name,".dat"),append=TRUE)
write(paste("TypeOfExperiment=",types_of_experiment),file=paste(file_name,".dat"),append=TRUE)

#modularity
write(paste("Initial_Modularity_Value=",initial_modularity),file=paste(file_name,".dat"),append=TRUE)
write(paste("Final_Modularity_Value=",final_modularity),file=paste(file_name,".dat"),append=TRUE)
#safeness
write(paste("Initial_Safeness_Value=",initial_safeness),file=paste(file_name,".dat"),append=TRUE)
write(paste("Final_Safeness_Value=",final_safeness),file=paste(file_name,".dat"),append=TRUE)

##initial deception
write(paste("Initial_Deception_Score=",initialdeceptionScore),file=paste(file_name,".dat"),append=TRUE)
#modularity
write(paste("Final_Deception_Score_MODULARITY=",round(final_deception_score_modularity,3)),file=paste(file_name,".dat"),append=TRUE)
#safeness
write(paste("Final_Deception_Score_SAFENESS=",round(final_deception_score_safeness,3)),file=paste(file_name,".dat"),append=TRUE)
###TIME
write(paste("TIME_MODULARITY_ALGO=",modTime[3]),file=paste(file_name,".dat"),append=TRUE)
write(paste("TIME_SAFENESS_ALGO=",safTime[3]),file=paste(file_name,".dat"),append=TRUE)
#
write("\n",file=paste(file_name,".dat"),append=TRUE)


##write results tab-separated for automatic import
write(paste(dataset,current_algo,budget,initial_modularity,final_modularity,initial_safeness, final_safeness,initialdeceptionScore,
            final_deception_score_modularity,final_deception_score_safeness,modTime[3],safTime[3],detectionTime[3],sep="\t"),
            file=paste(file_name_for_EXCEL,".dat"),append=TRUE)
}
