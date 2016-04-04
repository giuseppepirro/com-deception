library(igraph)
##BEST UPDATE 
getBestUpdateModularity<-function(members, network,nodes_in_my_community)
{
  
  ##IDs of communities sorted by degree
  degree_coms<-getDegreeSortedCommunities(members,network)
  
  

inter_add_delta<-NULL
intra_del_delta<-NULL

be<-finNodeHighestDegCom(members,network,nodes_in_my_community,degree_coms)

if(!is.null(be[1])){
  node1<-be[1]
  nodes_in_com2<-getNodesinCommunity(members,be[2])
  idxn2<-ceiling(runif(1, 1, length(nodes_in_com2)))
  node2<- nodes_in_com2[idxn2]
  ##########if the two nodes are already connected pick another node2
  while(are.connected(network,node1,node2)==TRUE)
  {
    idxn2<-ceiling(runif(1, 1, length(nodes_in_com2)))
    node2<- nodes_in_com2[idxn2]
  }
  
  inter_add_delta<- getDeltaModularityAfterAddition(members, network, node1, node2) 
}else inter_add_delta<-0


best_com<-computeBestCommunityForDeletionForModularity(members,network, nodes_in_my_community)

if(!is.null(best_com)){
  node_lowest<-getTargetNodeInComWithIntraEdges(members,network,best_com,nodes_in_my_community)

##get a random node (neighbour of node_lowest) IN THE SAME COMMUNITY so we are sure that the edge exists
nodeL2<-getAnotherSameCom(members,network,node_lowest)
}
else node_lowest<-NULL

#####We need to be sure that the edge exists!!!
if(!(is.null(node_lowest))) intra_del_delta<- getDeltaModularityAfterDeletion(members, network, node_lowest, nodeL2) 
else intra_del_delta<-0

best<-NULL
best[1]<-NULL
best[2]<-NULL

if(intra_del_delta<0)
{
  paste("ERROR NEGATIVE DELTA INTRA DEL EDGE MODULARITY")
}

if(inter_add_delta<0)
{
  paste("ERROR NEGATIVE DELTA INTER ADD EDGE MODULARITY")
}

if(inter_add_delta==0 && intra_del_delta==0)
  return (NULL)
  
if(inter_add_delta>=intra_del_delta)
{
  best[1]<-node1
  best[2]<-node2
  best[3]<-666 #corresponds to add
  # print(paste("BETTER INTER ADD WITH DELTA",inter_add_delta,"n1=",node1,"n2=",node2))
  
}
else 
  {
    best[1]<-node_lowest
    best[2]<-nodeL2
    best[3]<-999 #corresponds to del
    
    # print(paste("BETTER INTRA DEL WITH DELTA",intra_del_delta,"n1=",node_lowest,"n2=",nodeL2))
    
    
  }
#}

#print("Riga 139 - ritorno il cambio modMin")


return (best)
}
