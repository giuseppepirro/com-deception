library(igraph)

##BEST UPDATE 
getBestUpdateSafeness<-function(members, network, nodes_in_my_community)
{
  
  ##IDs of communities sorted by size
  size_coms<-getComSortedBySize(members,network)
  
  
  inter_add_delta<-NULL
  intra_del_delta<-NULL
  
  ## node in \comH with the maximum ratio |E(u,Ci)|/deg(u) having a possible inter-edge addition
  nodeA1<-getBestMemberByExternalEdges(members,network,nodes_in_my_community)

  
  ##highest and lowest size communities
  #nodes_in_lowest_size<-getNodesinCommunity(members,size_coms[1])
  #nodes_in_highest_size<-getNodesinCommunity(members,size_coms[length(size_coms)])
  
  ##IDs of communities sorted by degree
  degree_coms<-getDegreeSortedCommunities(members,network)

  if(!is.null(nodeA1)){
    com2<-length(degree_coms)
    nodeA2<-NULL
    while(is.null(nodeA2) && com2>0){
      nodes_in_com2<-getNodesinCommunity(members,degree_coms[com2])
      if(!(nodeA1 %in% nodes_in_com2)){
         neighs<-neighbors(network,nodeA1)
         nodes_in_com2<-setdiff(nodes_in_com2,neighs)
         if(length(nodes_in_com2)>0){
            idxn2<-ceiling(runif(1, 1, length(nodes_in_com2)))
            nodeA2<- nodes_in_com2[idxn2]
         }
      }else
        com2<-com2-1
    }
    inter_add_delta<- getDeltaSafenessScoreAfterAddition(members, network, nodes_in_my_community,nodeA1, nodeA2) 
  }else inter_add_delta<-0
  
#  print(paste("BEST ADD SAFE ",nodeA1, " ",nodeA2, " DELTA ", inter_add_delta))

  intra_community_edges_count<-getIntraEdgesCountSetOfNodes(network,nodes_in_my_community)
  
  nodeD1<-NULL
  nodeD2<-NULL
  
#  print(paste("intra_community_edges_count ",intra_community_edges_count, " NODES=",(length(nodes_in_my_community)-1)))
  
  if(intra_community_edges_count>length(nodes_in_my_community)-1){
    
    #The intra-edge with the minimal ratio on the endpoint 
    be<-computeBestDeletionForSafeness(members,network, nodes_in_my_community)
    nodeD1<-be[1]
    nodeD2<-be[2]
    intra_del_delta<- getDeltaSafenessScoreAfterDeletion(members, network, nodes_in_my_community, nodeD1, nodeD2) 
  #  print(paste("BEST DEL SAFE ",nodeD1, " ",nodeD2, " DELTA ", intra_del_delta))
  }
  
  best<-NULL
  best[1]<-NULL
  best[2]<-NULL
  
  if(is.null(nodeD1) && is.null(nodeA1)){
    print("NO EDIT POSSIBLE")
    return(NULL)
  }
  
  if(is.null(nodeD1)){
    best[1]<-nodeA1
    best[2]<-nodeA2
    best[3]<-666 #corresponds to add
    return(best)
  }
  
  if(is.null(nodeA1)){
    best[1]<-nodeD1
    best[2]<-nodeD2
    best[3]<-999 #corresponds to add
    return(best)
  }
  
  if(intra_del_delta<0)
  {
    paste("ERROR NEGATIVE DELTA INTRA DEL EDGE MODULARITY")
  }
  
  if(inter_add_delta<0)
  {
    paste("ERROR NEGATIVE DELTA INTER ADD EDGE MODULARITY")
  }
  
  
  if(inter_add_delta>=intra_del_delta)
  {
    best[1]<-nodeA1
    best[2]<-nodeA2
    best[3]<-666 #corresponds to add
    return(best)
    
  #  print(paste("BETTER INTER ADD WITH DELTA",inter_add_delta,"n1=",node1,"n2=",node2))
  }
  else 
  {
    best[1]<-nodeD1
    best[2]<-nodeD2
    best[3]<-999 #corresponds to del
    return(best)
    
  }
  
  
}
