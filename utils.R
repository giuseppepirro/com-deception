getSquaredDegreePerCom<-function(memberships,network)
  {
output<-0
for (i in 1:max(memberships)) {
  output[i]=0;
  
}

for (i in 1:length(memberships)) 
{
  output[memberships[i]]<-output[memberships[i]]+(degree(network,i))   
  
}
for (i in 1:max(memberships)) {
  output[i]<-output[i]^2;
}

return(output)
}

##get the community (integer) to which a node belongs
getCommunity<-function(memberships,node)
  {
  return(memberships[node])
  }


##check if two nodes are in the same community
isSameCommunity<-function(memberships,node1,node2)
{
  if(memberships[node1]==memberships[node2])
  return(TRUE)
  else
    return(FALSE)
}


## get the number of intra edges for a node
getIntraEdgesCount<-function(memberships,network,node)
  {
  #get all edges of a node
  #check among the nodes reached by an edge is these are in the same community
  
  count<-0
  neighs<-neighbors(network,node)
  
  for(n in neighs)
  {
    if(isSameCommunity(memberships,n,node))
       count<-count+1
  }
    
  return(count)
  #print(neighs)
}

##get the number of intra edges for a node and a given community
getIntraEdgesCountC<-function(network,nodes,node)
{
  #get all edges of a node
  #check among the nodes reached by an edge is these are in the same community
  
  count<-0
  neighs<-neighbors(network,node)
  
  for(n in neighs)
  {
    if(n %in% nodes)
      count<-count+1
  }
  
  return(count)
  #print(neighs)
}


##get the number of inter edges for a node and a given community
getInterEdgesCountC<-function(memberships,network,nodes,node)
{
  count<-0
  neighs<-neighbors(network,node)
  
  for(n in neighs)
  {
    if(!(n %in% nodes))
      count<-count+1
  }
  
  return(count)
  #print(neighs)
}

##get the number of inter edges for a node
getInterEdgesCount<-function(memberships,network,node)
{
  #get all edges of a node
  #check among the nodes reached by an edge is these are in the same community
  
  count<-0
  neighs<-neighbors(network,node)
  
  for(n in neighs)
  {
    if(!isSameCommunity(memberships,n,node))
      count<-count+1
  }
  
  return(count)
  #print(neighs)
}

##get the number of inter edges for the entire community
getInterEdgesCountCom<-function(memberships,network,comm)
{
  #get all edges of a node
  #check among the nodes reached by an edge is these are in the same community
  
  count<-0
  
  for (i in 1:length(memberships)) 
  {
    if(memberships[i]==comm)
      count<-count+getInterEdgesCount(memberships,network,i)
  }
  return(count)
  #print(neighs)
}

##get the number of inter edges for the entire community
getIntraEdgesCountCom<-function(memberships,network,comm)
{
  #get all edges of a node
  #check among the nodes reached by an edge is these are in the same community
  
  count<-0
  
  for (i in 1:length(memberships)) 
  {
    if(memberships[i]==comm)
      count<-count+getIntraEdgesCount(memberships,network,i)
  }
  count<-count/2
  return(count)
}

##return the number of nodes in a community
 getNodesinCommunityCount <- function(memberships, com_number)
{
   count <- 0
   for (i in 1:length(memberships))
   {
     if (identical(memberships[i],com_number))
       count <- count + 1
   }
   return(count)
   
 }

##return the number of nodes in a community
getMaxCommunitySize<-function(memberships)
{
  
  tatal_num_com<-max(memberships)
  
  membershipCount<-vector(mode="numeric",length=tatal_num_com)

  
  for (i in 1:tatal_num_com) 
  {
    membershipCount[i]<-0
  }
  
  for (j in 1:length(memberships)) 
  {
    membershipCount[memberships[j]]<-membershipCount[memberships[j]]+1
 
  }

  max<- membershipCount[1]
  
  for (k in 1:tatal_num_com) 
  {
    
    if(membershipCount[k]>max)
    {
     # print(membershipCount[i])
      max<-membershipCount[k]
    }
    
    
  }
  
  return(max)
  
}

##return the number of nodes in a community
getNodesinComOrderedByDegree<-function(memberships,network, com_number)
{
  nodesInCom<-getNodesinCommunity(memberships, com_number)
  
  degsNodes=c()
  
  for (i in 1:length(nodesInCom)) 
  {

    degsNodes[i]=degree(network,nodesInCom[i])
  }
  
  
  degsNodes<-sort(degsNodes, decreasing = FALSE)
  
  
  return(degsNodes)
  
}

##return the nodes in a community
getNodesinCommunity<-function(memberships, com_number)
{
  nodesinCom=c()
  
  for (i in 1:length(memberships)) 
  {
    if(memberships[i]==com_number)
    {
      nodesinCom<-c(nodesinCom,i)
    }
    
  }
  return(nodesinCom)
  
}



##return the degree of a community
getCommunityDegree<-function(memberships,network, com_number)
{
  degCom<-0
  
  nodesinC=getNodesinCommunity(memberships, com_number)
  
  for (n in nodesinC) 
  {
    degCom<-degCom+degree(network,n)
    
  }
  
  return(degCom)
  
}


#return the community with the highest degree
##return the degree of a community and not the COMMUNITY ID!!
getHighestDegreeCommunity<-function(memberships,network)
{
  coms<-getDegreeSortedCommunities(memberships,network)
  return(length(coms))
 
}


##get communities sorted by degree AT THE MOMENT DOES NOT WORK! returns the degree but not the community id
getDegreeSortedCommunities<-function(memberships,network)
{
  degrees_coms<-0
  for (k in 1:max(memberships)) {
    degrees_coms[k]=0;
  }
  
  for (k in 1:length(degrees_coms)) 
  {
    degrees_coms[k]<-getCommunityDegree(memberships,network,k)
  }
  
  
   ordered_vector<- unlist(as.matrix(sort.int(decreasing = TRUE,degrees_coms, index.return=TRUE))[1])
   ordered_vector_indexes<- unlist(as.matrix(sort.int(decreasing = TRUE,degrees_coms, index.return=TRUE))[2])
   
  
 ## print(ordered_vector)
  ##print(ordered_vector_indexes)
  
  
  outputS<-ordered_vector_indexes
  return(outputS)
}

##
getModularity<-function(memberships,network)
{
  numIntraEdges<-0
  degree<-getSquaredDegreePerCom(memberships,network)
  commDegree<-0
  for (i in 1:max(memberships)){
    numIntraEdges<-numIntraEdges+getIntraEdgesCountCom(memberships,network,i)
    commDegree<-commDegree+degree[i]
  }
  modularity<-(numIntraEdges/ecount(network))-(commDegree/(4*ecount(network)^2))
  return(modularity)
}


getDeltaModularityAfterDeletion <- function(memberships, network, n1, n2)
{
  modularity<-getModularity (memberships, network)
  
  network1<-network
  
  network1<- delete_edges(network1,paste(n1,"|",n2))
 
  
  modularity2<-getModularity (memberships, network1)
  diff<-modularity-modularity2
  return(diff)
}

getDeltaModularityAfterAddition <- function(memberships, network, n1, n2)
{
  modularity<-getModularity (memberships, network)
  network1<-network
  network1<- add_edges(network1,c(n1,n2))

  modularity2<-getModularity (memberships, network1)
  diff<-modularity-modularity2
  return(diff)
}

computeBestDeletionForModularity <- function(memberships, network, nodes)
{
  best<-NULL
  best[1]<-NULL
  best[2]<-NULL
  bestDelta<-NULL
  for(i in 1:(length(nodes)-1)){
    for(j in (i+1):length(nodes))
      {
      res<-are.connected(network,nodes[i], nodes[j])
      if(res)
        {
           delta<-getDeltaModularityAfterDeletion(memberships,network,nodes[i],nodes[j])
        
        if((is.null(bestDelta) || delta>bestDelta) && delta>0)
          {
          best[1]<-nodes[i]
          best[2]<-nodes[j]
          bestDelta<-delta
          
        }
      }
    }
  }
  return(best)
}

computeBestAdditionForModularity <- function(memberships, network, nodes)
{
  best<-NULL
  best[1]<-NULL
  best[2]<-NULL
  bestDelta<-NULL
  otherNodes<-NULL
  j<-1
  for (i in 1:length(memberships)){
    if(! i %in% nodes){
      otherNodes[[j]]<-i
      j<-j+1
    }
  }
  
  for(i in 1:length(nodes)){
    for(j in 1:length(otherNodes)){
      if(are.connected(network,nodes[i],otherNodes[j])==FALSE){
        delta<-getDeltaModularityAfterAddition(memberships,network,nodes[i],otherNodes[j])
        if((is.null(bestDelta) || delta>bestDelta) && delta>0){
          best[1]<-nodes[i]
          best[2]<-otherNodes[j]
          bestDelta<-delta
        }
      }
    }
  }
  return(best)
}


#get number of reachable nodes in the same community by passing through nodes only in the community
getReachableInCom<-function(memberships,network,node)
{
  reachable<-NULL
  num_reachable<-1
  for (i in 1:length(memberships)) 
  {
    reachable[i]<-FALSE   
  }
  queue<-new.queue()
  visited<-NULL
  for (i in 1:length(memberships)) 
  {
    visited[i]<-FALSE   
  }
  
  enqueue(queue,node)
  visited[node]<- TRUE
  
  while (! is.empty(queue)) {
    current<-dequeue(queue)
    for (n in neighbors(network,current)){
      if ((visited[n] == FALSE) && isSameCommunity(memberships,node,n)){
        num_reachable<-num_reachable+1
        reachable[n]<-TRUE
        visited[n]<- TRUE
        enqueue(queue,n)
      }
    }
  }
  
  return(num_reachable)
}


getReachable<-function(network,nodes,node)
{
  reachable<-NULL
  num_reachable<-1
  for (i in 1:length(nodes)) 
  {
    reachable[i]<-FALSE   
  }
  queue<-new.queue()
  visited<-NULL
  for (i in 1:length(nodes)) 
  {
    visited[i]<-FALSE   
  }
  
  enqueue(queue,node)
  pos<-match(node,nodes)
  visited[pos]<- TRUE
  
  while (! is.empty(queue)) {
    current<-dequeue(queue)
   
    for (n in neighbors(network,current))
      if(n %in% nodes){
        pos<-match(n,nodes)
        if (visited[pos] == FALSE){
          num_reachable<-num_reachable+1
          reachable[pos]=TRUE
          visited[pos]<- TRUE
          enqueue(queue,n)
        }
      }
  }
  
  return(num_reachable)
}

#implementation of a queue
{
  new.queue <- function() {
    ret <- new.env()
    ret$front <- new.env()
    ret$front$q <- NULL
    ret$front$prev <- NULL
    ret$last <- ret$front
    return(ret)
  }
  
  ## add to end of queue
  enqueue <- function(queue, add){
    queue$last$q <- new.env()
    queue$last$q$prev <- queue$last
    queue$last <- queue$last$q
    queue$last$val <- add
    queue$last$q <- NULL
  }
  
  ## return front of queue and remove it
  dequeue <- function(queue){
    if (is.empty(queue)) {
      stop("Attempting to take element from empty queue")
    }
    value <- queue$front$q$val
    queue$front <- queue$front$q
    queue$front$q$prev <- NULL
    return(value)
  }
  
  is.empty <- function(queue){
    return(is.null(queue$front$q))
  }
}

#get number of communities with nodes

getCommunitySpread <- function(memberships,nodes)
{
  presence<-FALSE
  for (i in 1:max(memberships)) {
    presence[i]=FALSE;
  }
  
  num_comm<-0
  
  for (n in nodes) 
  {
    if(presence[memberships[n]]==FALSE){
      num_comm<-num_comm+1
      presence[memberships[n]]<-TRUE
    }
  }
  
  return(num_comm)
}

#get average portions of nodes per community 

getCommunityExposition <- function(memberships,network,nodes)
{
  presence<-NULL
  size_comm<-0
  for (i in 1:max(memberships)) {
    presence[i]<-0
    size_comm[i]<-0
  }
  
  
  
  for (n in V(network)) 
  {
    size_comm[memberships[n]]<-size_comm[memberships[n]]+1
    if(n %in% nodes){
      presence[memberships[n]]<-presence[memberships[n]]+1
    }
  }
  
  average<-0
  count<-0
  
  for (i in 1:max(memberships)) {
    if(presence[i]>0){
      average<-average+(presence[i]/size_comm[i])
      count<-count+1
    }
  }
    
    average<-average/count
  
  return(average)
}

getSafenessScore <- function(memberships, network, nodes)
{
  safe_score<-0
  max_size<-getMaxCommunitySize(memberships)
  reachable<-c()
  result<-connectedComponents(network, nodes)

  for(k in 1:(result[[1]])){
    for(i in 1:(length((result[[k+1]]))))
      reachable[result[[k+1]][[i]]]<-length((result[[k+1]]))
    
  }

  for(n in nodes){
    safe_score<-safe_score+getSafenessScoreNode(memberships, network, nodes, n,reachable[n],max_size)
  }
  
  safe_score<-(safe_score/length(nodes))

  return(safe_score)
}

getSafenessScoreWithReachability <- function(memberships, network, nodes,reachable)
{
  safe_score<-0
  max_size<-getMaxCommunitySize(memberships)
  
  for(n in nodes){
    safe_score<-safe_score+getSafenessScoreNode(memberships, network, nodes, n,reachable[n],max_size)
  }
  
  safe_score<-(safe_score/length(nodes))

  return(safe_score)
}

getdeceptionScore <- function(memberships, network, nodes)
{
  deception_score<-0
  
  
  deception_score<-((1/2)*((getCommunitySpread(memberships, nodes)-1)/length(nodes))+(1/2)*(1-getCommunityExposition(memberships,network,nodes)))*(1-((connectedComponents(network,nodes)[[1]]-1)/(length(nodes)-1)))

  return(deception_score)
}

getSafenessScoreNode <- function(memberships, network, nodes, node,reachable,max_size)
{
  
  safe_score<-0.5*((reachable-getIntraEdgesCountC(network,nodes,node))/length(nodes))+0.5*(getInterEdgesCountC(memberships,network,nodes,node)/(degree(network,node)))
  
  
  if(is.na(safe_score))
    safe_score=0
  return(safe_score)
}



getDeltaSafenessScoreAfterDeletion <- function(memberships, network, nodes, n1, n2)
{
  safe_score<-getSafenessScore (memberships, network, nodes)
  
  network1<-NULL
  
  network1<-network
  
  network1<- delete_edges(network1,paste(n1,"|",n2))

  
  safe_score2<-getSafenessScore (memberships, network1, nodes)
  diff<-safe_score2-safe_score
  
  return(diff)
}


getDeltaSafenessScoreAfterDeletionWithReachability <- function(memberships, network, nodes, n1, n2,reachable)
{
  safe_score<-getSafenessScoreWithReachability (memberships, network, nodes,reachable)
  
  network1<-NULL
  
  network1<-network
  
  network1<- delete_edges(network1,paste(n1,"|",n2))
  
  reachableN<-c()
  
  result<-connectedComponents(network1, nodes)
  
  for(k in 1:(result[[1]])){
    for(i in 1:(length((result[[k+1]]))))
      reachableN[result[[k+1]][[i]]]<-length((result[[k+1]]))
    
  }
  
  safe_score2<-getSafenessScoreWithReachability (memberships, network1, nodes, reachableN)
  diff<-safe_score2-safe_score
  return(diff)
}
        
        

getDeltaSafenessScoreAfterAddition <- function(memberships, network, nodes, n1, n2)
{
  safe_score<-getSafenessScore (memberships, network, nodes)
  
  network1<-network
  
  network1<- add_edges(network1,c(n1,n2))
  
  safe_score2<-getSafenessScore (memberships, network1, nodes)
  
  diff<-safe_score2-safe_score
  
  return(diff)
}

computeBestCommunityForDeletionForModularity <- function(memberships, network, nodes)
{
  bestDegree<-NULL
  bestCom<-NULL
  
  for(i in 1:(length(nodes))){
        neighs<-neighbors(network,nodes[i])
        
        my_com<-getCommunity(memberships,nodes[i])
        
        nodes_in_com<-getNodesinCommunity(memberships,my_com)  
        
        neigh_same_com<-intersection(neighs,nodes_in_com)
        
        if(length(neigh_same_com)>0){
           degree<- getCommunityDegree(memberships,network, memberships[nodes[i]])
           if(is.null(bestDegree) || degree<bestDegree){
              bestDegree<-degree
              bestCom<-memberships[nodes[i]]
           }
        }
    }
  return(bestCom)
}

computeBestDeletionForSafeness <- function(memberships, network, nodes)
{
  best<-NULL
  best[1]<-NULL
  best[2]<-NULL
  bestDelta<-NULL
  #print(nodes)
  reachable<-c()
  
  result<-connectedComponents(network, nodes)
  
  for(k in 1:(result[[1]])){
    for(i in 1:(length((result[[k+1]]))))
      reachable[result[[k+1]][[i]]]<-length((result[[k+1]]))
    
  }
  
   for(i in 1:(length(nodes)-1)){
    for(j in (i+1):length(nodes)){
      if(are.connected(network,nodes[i], nodes[j])==TRUE)
        {
        delta<-getDeltaSafenessScoreAfterDeletionWithReachability(memberships,network,nodes,nodes[i],nodes[j],reachable)

        if((is.null(bestDelta) || delta>bestDelta) && delta>0){
          best[1]<-nodes[i]
          best[2]<-nodes[j]
          bestDelta<-delta
        }
      }
     
      
    }
   
  }
  return(best)
}



computeBestAdditionForSafeness <- function(memberships, network, nodes)
{
  best<-NULL
  best[1]<-NULL
  best[2]<-NULL
  bestDelta<-NULL
  otherNodes<-NULL
  j<-1
  for (i in 1:length(memberships)){
    if(! i %in% nodes){
      otherNodes[[j]]<-i
      j<-j+1
    }
  }
  
  for(i in 1:length(nodes)){
    for(j in 1:length(otherNodes)){
      if(are.connected(network,nodes[i],otherNodes[j])==FALSE){
        delta<-getDeltaSafenessScoreAfterAddition(memberships,network,nodes,nodes[i],otherNodes[j])
        if((is.null(bestDelta) || delta>bestDelta) && delta>0){
          best[1]<-nodes[i]
          best[2]<-otherNodes[j]
          bestDelta<-delta
        }
      }
    }
  }
  return(best)
}

#####
computeCommunities <- function(network, algo_name)
{
  communities<-0;
  
  if(algo_name=="louv")
  {
    communities<-cluster_louvain(network)
  }
  else if(algo_name=="opt")
    {    
    communities<-cluster_optimal(network)

    }
  else if(algo_name=="infomap")
  {
    communities<-cluster_infomap(network)
    
  }
  else if(algo_name=="walk")
  {
    communities<-cluster_walktrap(network)
    
  }
  else if(algo_name=="greedy")
  {
    communities<-cluster_fast_greedy(network)
    
  }
  else if(algo_name=="spin")
  {
    communities<-cluster_spinglass(network)
    
  }
  else if(algo_name=="labelp")
  {
    communities<-cluster_label_prop(network)
    
  }
  else if(algo_name=="betw")
  {
    communities<-cluster_edge_betweenness(network)
    
  }
  else if(algo_name=="eig")
  {
    communities<-cluster_leading_eigen(network)
    
  }
    
  return(communities)
  
}

initializeNetwork<- function(type_of_net)
  {
  
  net<-0
  
if(identical(type_of_net,"Zachary"))
{
  net<- graph.famous("Zachary")
} 
  else
    if(identical(type_of_net,"Kite"))
    {
      net<- graph.famous("Krackhardt_Kite")
    }
  else
    if(identical(type_of_net,"Erdos"))
       {
      #n: The number of vertices in the graph.
      #p.or.m Either the probability for drawing an edge between two arbitrary vertices (G(n,p) graph), or the number of edges in the graph (for G(n,m) graphs).
      #type: The type of the random graph to create, either gnp (G(n,p) graph) or gnm (G(n,m) graph).
      #directed: Logical, whether the graph will be directed, defaults to FALSE.
      #loops: Logical, whether to add loop edges, defaults to FALSE.
      n<-500
      p<-1/1000
      net<-erdos.renyi.game(n, p,directed = FALSE)
    }
  else
    if(identical(type_of_net,"Watts"))
    {
      #dim: Integer constant, the dimension of the starting lattice.
      #size: Integer constant, the size of the lattice along each dimension.
      #nei: Integer constant, the neighborhood within which the vertices of the lattice will be connected.
      #p: Real constant between zero and one, the rewiring probability.
      #loops: Logical scalar, whether loops edges are allowed in the generated graph.
      #multiple: Logical scalar, whether multiple edges are allowed int the generated graph.
      dim<-1
      size<-500
      nei<-5
      p<-0.05
      net<-watts.strogatz.game(dim, size, nei, p,directed = FALSE)
    }
  else
    if(identical(type_of_net,"Barabasi"))
    {
      #n: Number of vertices.
      #power: The power of the preferential attachment, the default is one, ie. linear preferential attachment.
      #m: Numeric constant, the number of edges to add in each time step This argument is only used if both out.dist and out.seq are omitted or NULL.
      n<-500
      power<-.3
      net<-barabasi.game(n,power,directed = FALSE)
    }
  else
    if(identical(type_of_net,"madrid"))
    {
    
      file_name<-paste("./datasets/",type_of_net,".edgelist",sep = "")
      net<- as.undirected(read.graph(file_name,format=c("edgelist")) )
      }
  else
    if(identical(type_of_net,"oclinks"))
    {
      
      file_name<-paste("./datasets/",type_of_net,".edgelist",sep = "")
      print(file_name)
      net<- as.undirected(read.graph(file_name,format=c("edgelist")) )
    }
  else
    if(identical(type_of_net,"facebook-75"))
    {
      
      file_name<-paste("./datasets/",type_of_net,".edgelist",sep = "")
      print(file_name)
      net<- as.undirected(read.graph(file_name,format=c("edgelist")) )
    }
  else
    if(identical(type_of_net,"ERDOS992"))
    {
      
      file_name<-paste("./datasets/",type_of_net,".net",sep = "")
      print(file_name)
      net<- as.undirected(read.graph(file_name,format=c("pajek")) )
    }
  else
    if(identical(type_of_net,"jazz"))
    {
      
      file_name<-paste("./datasets/",type_of_net,".net",sep = "")
      print(file_name)
      net<- as.undirected(read.graph(file_name,format=c("pajek")) )
    }
  else
    if(identical(type_of_net,"USAir97"))
    {
      
      file_name<-paste("./datasets/",type_of_net,".net",sep = "")
      net<- as.undirected(read.graph(file_name,format=c("pajek")) )
    }
  else
    if(identical(type_of_net,"email"))
    {
      
      file_name<-paste("./datasets/",type_of_net,".edgelist",sep = "")
      print(file_name)
      net<- as.undirected(read.graph(file_name,format=c("edgelist")) )
    }
  else
    if(identical(type_of_net,"ia-infect-dublin"))
    {
      
      file_name<-paste("./datasets/",type_of_net,".edgelist",sep = "")
      net<- as.undirected(read.graph(file_name,format=c("edgelist")) )
    }
  else
    if(identical(type_of_net,"generated"))
  {
      file_name<-paste("./datasets/",type_of_net,".edgelist",sep = "")
      print(file_name)
      net<- as.undirected(read.graph(file_name,format=c("edgelist")) )
  }
    else
    {
    ##READ FROM FILE AND RETURN THE NET!
file_name<-paste("./datasets/",type_of_net,".gml",sep = "")
    net<-read.graph(file_name,format=c("gml"))
      }
  
  return(net)
  
}

##Method to chose the community that we want to protect in the WORST-CASE
getCommunityToProtectWorstCase<-function(initialM,num_communties)
{
  members<-c()
  
 while(length(members)<4) 
  {
  com_number<- ceiling(runif(1, 1, total_num_communities))
  members<-getNodesinCommunity(initialM,com_number)
  }
  
  return(members)
  }
  
##Method to chose the community that we want to Protect
getCommunityToProtect<-function(initialNetwork,total_num_nodes,type_of_size)
{
    
  size_com<-0
  
  if(identical(type_of_size,"Small"))
  {
    size_com<-ceiling(0.05*total_num_nodes)
    print(size_com)
  }
  else
    if(identical(type_of_size,"Medium"))
    {
      size_com<-ceiling(0.1*total_num_nodes)
      
      
    }
  else
    if(identical(type_of_size,"Big"))
    {
      size_com<-ceiling(0.15*total_num_nodes)
      
    }
  members<-vector(mode="numeric",length=size_com)
  
  for(i in 1: size_com)
  {
    candidate<-ceiling(runif(1, i, total_num_nodes))
    print(append("Candidate",candidate))
    members<-c(members,candidate)
  }
  
  members<-members[-which(members %in% c(0))]
  
  return(members)
}


##Find the node in \comH also belonging to the community with the highest degree (useful for mod min)
finNodeHighestDegCom<-function(members,network,my_nodes,degree_coms)
{    
  my_n<-NULL
  my_n[1]<-NULL
  my_n[2]<-NULL
  
  com1<-1
  while(com1<=length(degree_coms)){
     nodes_in_com<-getNodesinCommunity(members,degree_coms[com1])
  
     my_nodes_in_com<-intersect(my_nodes,nodes_in_com)
  
     while(length(my_nodes_in_com)>0){
       node1<-getTargetNodeInCom(members,degree_coms[com1],my_nodes)
       com2<-finHighestDegComForAddition(members,network,node1,degree_coms)
       if(!is.null(com2)){
         my_n[1]<-node1
         my_n[2]<-com2
         return (my_n)
       }
   #    print(my_nodes_in_com)
       my_nodes_in_com<-setdiff(my_nodes_in_com,node1)
   #    print(paste("AFTER DELETION OF ",node1," ",my_nodes_in_com))
     }
     com1<-com1+1
   }
   
  return(my_n)
}

finHighestDegComForAddition<-function(members,network,node,degree_coms)
{
  
  res<-NULL
  i<-0
  cont<-TRUE
  while(cont && i<length(degree_coms))
  {    
    i<-i+1
    if(members[node]!=degree_coms[i]){
      nodes_in_Bcom<-getNodesinCommunity(members,degree_coms[i])
    
      for(j in 1:length(nodes_in_Bcom)){
        if(are.connected(network,node,nodes_in_Bcom[j])==FALSE){
          cont=FALSE
          break
        }
      }
    }
  }
  if(i<length(degree_coms)) return(degree_coms[i])
  else if(i==length(degree_coms) && !cont) return (degree_coms[i])
  else return (NULL)
}
###################################

##Find the nodein \comH belonging to the community with the highest degree (useful for mod min)
finNodeLowestDegCom<-function(members,network,my_nodes)
{
  degree_coms<-getDegreeSortedCommunities(members,network)
  res<-NULL
  i<-length(degree_coms)
  cont<-TRUE
  
  while(cont)
  {    
    my_n<-getTargetNodeInCom(members,degree_coms[i],my_nodes)

    
    if(!is.na(my_n))
    {
      res<-my_n
      cont<-FALSE
    }
    i<-i-1
  }
  return(my_n)
}

##return another node in the same community
getAnotherSameCom<-function(members,network,node)
{
  neighs<-neighbors(network,node)
  
  my_com<-getCommunity(members,node)
  
  nodes_in_com<-getNodesinCommunity(members,my_com)  
  
  neigh_same_com<-intersection(neighs,nodes_in_com)
  
  notMyID<-FALSE

  while(!notMyID)
  { 
    idx<-runif(1, 1, length(neigh_same_com))
   
     my_n<-neigh_same_com[idx]
   
     if(!is.na(my_n))
    {
      if(!identical(node,my_n))
      {
        res<-my_n
        notMyID<-TRUE
      }
    }

  }
  
  return(res)
}

##Given a community and the list of members of \comH's get a random node if contained
getTargetNodeInCom<-function(members,com_number,my_nodes)
{
  nodes_in_com<-getNodesinCommunity(members,com_number)
  
  my_nodes_in_com<-intersect(my_nodes,nodes_in_com)
  
  size<-length(my_nodes_in_com)
  
  idxn<-suppressWarnings(ceiling(runif(1, 1, size)))
  
  ##THE INTERSECTION CAN BE EMPTY: N/A IS CAPTURED IN THE METHODS THAT CALL THIS ONE
  node_res<- nodes_in_com[idxn]
  
return(node_res)
  
}

getTargetNodeInComWithIntraEdges<-function(members,network,com_number,my_nodes)
{
  nodes_in_com<-getNodesinCommunity(members,com_number)
  
  my_nodes_in_com<-intersect(my_nodes,nodes_in_com)
  
  size<-length(my_nodes_in_com)
  
  repeat{
  idxn<-suppressWarnings(ceiling(runif(1, 1, size)))
  
  ##THE INTERSECTION CAN BE EMPTY: N/A IS CAPTURED IN THE METHODS THAT CALL THIS ONE
  node_res<- nodes_in_com[idxn]
  
  neighs<-neighbors(network,node_res)
  
  neigh_same_com<-intersection(neighs,nodes_in_com)
  
  if(length(neigh_same_com)>0) break
  
 } 
    
  
  return(node_res)
  
}


##return the number of nodes in a community
getBestMemberByExternalEdges<-function(members,network, nodes)
{
  
  degsNodes=c()
  
  for (i in 1:length(nodes)) 
  {
    
    degsNodes[i]=getInterEdgesCountMyCommunity(network,nodes[i],nodes)/degree(network,nodes[i])
  }
  

  ordered_vector<- unlist(as.matrix(sort.int(decreasing = FALSE,degsNodes, index.return=TRUE))[1])
  ordered_vector_indexes<- unlist(as.matrix(sort.int(decreasing = FALSE,degsNodes, index.return=TRUE))[2])

  
  my_n<-NULL
  my_n[1]<-NULL
  my_n[2]<-NULL
  
  com1<-1

  while(com1<=length(ordered_vector_indexes)){

    neighs<-neighbors(network,nodes[ordered_vector_indexes[com1]])
    
    all_nodes<-V(network)
    
    all_nodes<-setdiff(all_nodes,nodes)
    
    nodes_not_in_com<-setdiff(all_nodes,neighs)
    
    if(length(nodes_not_in_com)>0){
      return (nodes[ordered_vector_indexes[com1]])
    }
    com1<-com1+1
  }
  
  return(NULL)
  
}

##get the number of inter edges for a node
getInterEdgesCountMyCommunity<-function(network,node,nodes)
{
  #get all edges of a node
  #check among the nodes reached by an edge is these are in my community
  
  count<-0
  neighs<-neighbors(network,node)
  
  for(n in neighs)
  {
    if(!(n %in% nodes))
      count<-count+1
  }
  
  return(count)
}


getComSortedBySize<-function(members,network)
{
  sizeCom=c()
  
  for (i in 1:max(members)) 
  {
    sizeCom[i]=0
  }
  
  for (i in 1:length(members)) 
  {
    sizeCom[members[i]]=sizeCom[members[i]]+1
  }
  
  
  ordered_vector<- unlist(as.matrix(sort.int(decreasing = FALSE,sizeCom, index.return=TRUE))[1])
  ordered_vector_indexes<- unlist(as.matrix(sort.int(decreasing = FALSE,sizeCom, index.return=TRUE))[2])
  
  outputS<-ordered_vector_indexes
  return(outputS)
}

##get the number of intra edges for a given community
getIntraEdgesCountSetOfNodes<-function(network,nodes)
{
  
  count<-0
  
  for(n1 in nodes){
    neighs<-neighbors(network,n1)
    
    for(n2 in neighs){
      if(n2 %in% nodes){
        count<-count+1
      }
    }
  }
  return(count/2)
}

getEdgeMinimalExternalRatio<-function(network,nodes)
{
  result<-NULL
  resultE<-vector(mode="numeric",length=2)
  
  for(n1 in nodes){
    neighs<-neighbors(network,n1)
    
    for(n2 in neighs){
      if(n2 %in% nodes){
        ratio=getInterEdgesCountMyCommunity(network,n1,nodes)/degree(network,n1)+getInterEdgesCountMyCommunity(network,n2,nodes)/degree(network,n2)
        if(is.null(result) || (ratio<result)){
          result<-ratio
          resultE[1]<-n1
          resultE[2]<-n2
        }
      }
    }
  }
  return(resultE)
}

connectedComponents<-function(network,nodes)
{
  visited<-c()
  num_components<-0
  components<-NULL
  
  for (i in 1:length(nodes)) 
    visited[i]<-FALSE
  
  for (i in 1:length(nodes)) {
    if(visited[i]==FALSE){
      component<-c()
      component[1]<-nodes[i]
      visited[i]<- TRUE
      num_components<-num_components+1
      queue<-new.queue()
      enqueue(queue,nodes[i])
      node_in_com<-2
      
      while (! is.empty(queue)) {
        current<-dequeue(queue)
        for (n in neighbors(network,current))
          if(n %in% nodes){
            pos<-match(n,nodes)
            if (visited[pos] == FALSE){
              visited[pos]<- TRUE
              component[node_in_com]<-n
              node_in_com<-node_in_com+1
              enqueue(queue,n)
            }
          }
      }
      components[num_components] <- list(component)
    }
  }
  

  result<-c()
  result[1]<-num_components
  for(i in 1:num_components)
    result[i+1]<-components[i]
  return(result)
  
}
