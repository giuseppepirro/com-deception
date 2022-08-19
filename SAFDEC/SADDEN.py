import igraph
from igraph import *
import copy
from scipy.stats import entropy
import numpy as np

def num_comm(target_comm, communities):
  uni_comm = []
  comm_list = []
  for node in target_comm:
    for c in communities:
      if node in c:
        comm_list.append(c)
        if c not in uni_comm:
          uni_comm.append(c)
          break
  return len(uni_comm), comm_list

def get_targetComm_Neighbours(target_comm, communities, Adjacency_List):
	List = []
	marked = dict()
	for i in target_comm:
		for j in Adjacency_List[i]:
			if j not in marked:
				for k in range(len(communities)):
					if j in communities[k]:
						List.append(k)
						marked[j] = j
	return List, marked

def check_neighbours(neighbours, communities):
	ctr = 0
	List = []
	for i in range(len(communities)):
		for j in communities[i]:
			if j in neighbours:
				List.append(i)
				ctr += 1
			if ctr == len(neighbours):
				return List
	return List

def get_entropy(labels, base = None):
  values, counts = np.unique(labels, return_counts = True)
  return entropy(counts, base = base)

def get_adj_list(E):
	Adjacency_List = {}
	for i in range (0, len(E)):
		e = E[i]
		s = e[0]
		t = e[1]
		if (s in Adjacency_List.keys()):
			Adjacency_List[s].append(t)
		else:
			Adjacency_List[s] = []
			Adjacency_List[s].append(t)
		if (t in Adjacency_List.keys()):
			Adjacency_List[t].append(s)
		else:
			Adjacency_List[t] = []
			Adjacency_List[t].append(s)
	return Adjacency_List

def com_decept_safeness(target_comm, communities, IG_edgeList, beta, deg, out_deg, out_ratio, num_vertices, new_adj, new_edge_list, intra_considered):
	add_gain = 0
	del_gain = 0
	while(True):
		node_list = get_min_NodeRatio_index(out_ratio)

		(add_gain, add_node_ind) = min_index_edge_addition(node_list, deg, out_deg)
		add_node = target_comm[add_node_ind]
		add_node_2 = findExternalNode(add_node, target_comm, communities, IG_edgeList)

		li = getBestDelExclBridges(target_comm, new_edge_list, new_adj, num_vertices)

		((del_node, del_node_2), max_gain) = deletion_Gain(li, intra_considered, deg, out_deg, target_comm)
		del_gain = max_gain

		if add_gain >= del_gain and add_gain > 0:
			IG_edgeList.append((add_node, add_node_2))
			for i in target_comm:
				deg_ = 0
				out_deg_ = 0
				for j in IG_edgeList:
					if i == j[0] or i == j[1]:
						deg_ = deg_ + 1
						if (i == j[0] and j[1] not in target_comm) or (i == j[1] and j[0] not in target_comm):
							out_deg_ = out_deg_ + 1
				deg[target_comm.index(i)] = deg_
				out_deg[target_comm.index(i)] = out_deg_

			for i in range (0, len(out_ratio)):
				out_ratio[i] = out_deg[i]/deg[i]
		elif del_gain > 0:
			IG_edgeList.remove((del_node,del_node_2))
			intra_considered.append((del_node,del_node_2))
			for i in target_comm:
				deg_ = 0
				out_deg_ = 0
				for j in IG_edgeList:
					if i == j[0] or i == j[1]:
						deg_ = deg_ + 1
						if (i == j[0] and j[1] not in target_comm) or (i == j[1] and j[0] not in target_comm):
							out_deg_ = out_deg_ + 1
				deg[target_comm.index(i)] = deg_
				out_deg[target_comm.index(i)] = out_deg_

			for i in range (0, len(out_ratio)):
				out_ratio[i] = out_deg[i]/deg[i]
			new_edge_list.remove((del_node, del_node_2))
			new_adj[del_node].remove(del_node_2)
			new_adj[del_node_2].remove(del_node)
		beta = beta - 1
		if (beta > 0 and (add_gain > 0 or del_gain > 0)):
			continue
		else:
			break
	return IG_edgeList

def get_min_NodeRatio_index(out_ratio):
	min_val = min(out_ratio)
	node = []
	for i in range (0, len(out_ratio)):
		if out_ratio[i] == min_val:
			node.append(i)
	return node

def min_index_edge_addition(node_list, deg, out_deg):
	node_ind = 0
	max_gain = 0
	for i in node_list:
		gain = 0.5*((out_deg[i]+1)/(deg[i]+1)-out_deg[i]/deg[i]) 
		if gain > max_gain:
			max_gain = gain
			node_ind = i
	return (max_gain, node_ind)

def findExternalNode(com_node, com, graph, edges):
	for i in graph:
		if i != com:
			for j in i:
				if ((com_node, j) or (j, com_node)) not in edges:
					return j
def deletion_Gain(li, intra_considered, deg, out_deg, target_comm):
	max_gain = 0
	node_u = 0
	node_v = 0
	for i in li:
		if i not in intra_considered:
			u = i[0]
			v = i[1]
			gain = (out_deg[target_comm.index(u)]/(2*deg[target_comm.index(u)]*(deg[target_comm.index(u)]-1)))+(out_deg[target_comm.index(v)]/(2*deg[target_comm.index(v)]*(deg[target_comm.index(v)]-1))) + (1/(len(target_comm) - 1))
			if(gain > max_gain):
				max_gain = gain
				node_u = u
				node_v = v
	return ((node_u, node_v), max_gain)

def getBestDelExclBridges(target_comm, edges, Adjacency_List, num_vertices):
	best_edges = []
	for i in edges:
		Cpy_Adj_List = copy.deepcopy(Adjacency_List)
		Cpy_Adj_List[i[0]].remove(i[1])
		Cpy_Adj_List[i[1]].remove(i[0])
		try:
			if(connectedComponents(target_comm, num_vertices, Cpy_Adj_List)) == 1:
				best_edges.append(i)
		except:
			continue
	return best_edges

def DFSUtil(target_comm, temp, v, visited, Adjacency_List):
	visited[v] = True
	temp.append(v)
	for i in Adjacency_List[target_comm[v]]:
		if visited[target_comm.index(i)] == False:
			temp = DFSUtil(target_comm, temp, target_comm.index(i), visited, Adjacency_List)
	return temp

def connectedComponents(target_comm, num_vertices, Adjacency_List):
	visited = [] 
	cc = [] 
	for i in range(num_vertices):
		visited.append(False)
	for v in range(num_vertices):
		if visited[v] == False: 
			temp = [] 
			cc.append(DFSUtil(target_comm, temp, v, visited, Adjacency_List))
	return len(cc)

def vertices_in_connectedComponents(target_comm, num_vertices, Adjacency_List, node):
	visited = [] 
	cc = [] 
	for i in range(num_vertices):
		visited.append(False)
	for v in range(num_vertices):
		if visited[v] == False: 
			temp = [] 
			cc.append(DFSUtil(target_comm, temp, v, visited, Adjacency_List))
	cc_node_list=[]
	for i in cc:
		tmp=[]
		for j in i:
			tmp.append(target_comm[j])
		cc_node_list.append(tmp)

	for i in cc_node_list:
		if (node in i):
			return len(i)
	return 0


def runSADDENSingleCom(graph,target_comm, communities,budget):
	edgesL = graph.es
	edges = []
	for e in edgesL:
		edges.append(e.tuple)
	e_ = edges
	Adjacency_List = get_adj_list(e_)
	num_vertices = graph.vcount()
	IG_edgeList = []
	for j in e_:
		IG_edgeList.append((j[0], j[1]))
	g = igraph.Graph(directed=False)
	g.add_vertices(num_vertices)
	g.add_edges(IG_edgeList)
	deg = []
	for j in target_comm:
		deg.append(g.vs[j].degree())
	out_deg = []
	for j in target_comm:
		out_ = 0
		for k in Adjacency_List[j]:
			if (k) not in target_comm:
				out_ = out_ + 1
		out_deg.append(out_)

	out_ratio = []
	for j in range(0, len(out_deg)):
		out_ratio.append(out_deg[j] / deg[j])
	new_adj = {}
	for j in Adjacency_List.keys():
		if j in target_comm:
			new_adj[j] = []
			for k in Adjacency_List[j]:
				if k in target_comm:
					new_adj[j].append(k)
	new_edge_list = []
	for j in IG_edgeList:
		if j[0] in target_comm and j[1] in target_comm:
			new_edge_list.append(j)
	if budget is 0:
		beta = int(0.3 * len(target_comm))
	else:
		beta=budget

	intra_considered = []
	IG_edgeList=list(set(IG_edgeList))
	new_edge_list=list(set(new_edge_list))
	IG_edgeList_ = com_decept_safeness(target_comm, communities, IG_edgeList, beta, deg, out_deg, out_ratio,
									   len(target_comm), new_adj, new_edge_list, intra_considered)

	return IG_edgeList_
##################################################




def SADDENALL():
	path="./dataset/"
	dataset="karate"
	graph = Graph.Read_Edgelist(path+dataset+".edgelist")

	Adjacency_List = graph.get_adjacency()

	g = igraph.Graph(directed = False)
	num_vertices = graph.vcount()
	g.add_vertices(num_vertices)

	graph = copy.deepcopy(g)
	graph2 = copy.deepcopy(g)
	e_ = list(graph.edges)
	e2_ = list(graph2.edges)

	IG_edgeList = []

	for i in e_:
		IG_edgeList.append((i[0], i[1]))

	IG_edgeList2 = IG_edgeList[:]

	g.add_edges(IG_edgeList)

	communities = g.community_walktrap().as_clustering()
	communities2 = copy.deepcopy(communities)
	comm_1 = copy.deepcopy(communities)
	safe_copy_comm = copy.deepcopy(communities)

	comm_length = len(communities)

	NMI_List = []
	Neighbourhood_NMI_List = []
	sum_comm = 0
	sum_entropy = 0

	for i in range(0, len(communities)):
		e_ = list(graph.edges)
		e2_ = list(graph2.edges)

		Adjacency_List = get_adj_list(e_)
		Adjacency_List2 = get_adj_list(e2_)

		g = igraph.Graph(directed = False)
		g2 = igraph.Graph(directed = False)

		num_vertices = 62

		g.add_vertices(num_vertices)
		g2.add_vertices(num_vertices)

		IG_edgeList = []

		for j in e_:
			IG_edgeList.append((j[0], j[1]))

		IG_edgeList2 = IG_edgeList[:]

		g.add_edges(IG_edgeList)
		g2.add_edges(IG_edgeList2)

		target_comm = communities[i]
		target_comm2 = communities2[i]
		pre_neighbours, neighbours = get_targetComm_Neighbours(target_comm, comm_1, Adjacency_List)

		deg = []
		for j in target_comm:
			deg.append(g.vs[j].degree())

		deg2 = copy.deepcopy(deg)

		out_deg = []

		for j in target_comm:
			out_ = 0
			for k in Adjacency_List[j]:
				if (k) not in target_comm:
					out_ = out_ + 1
			out_deg.append(out_)

		out_deg2 = copy.deepcopy(out_deg)

		out_ratio = []
		for j in range (0, len(out_deg)):
			out_ratio.append(out_deg[j]/deg[j])

		out_ratio2 = out_ratio[:]

		new_adj = {}
		for j in Adjacency_List.keys():
			if j in target_comm:
				new_adj[j] = []
				for k in Adjacency_List[j]:
					if k in target_comm:
						new_adj[j].append(k)

		new_adj2 = copy.deepcopy(new_adj)

		new_edge_list = []
		for j in IG_edgeList:
			if j[0] in target_comm and j[1] in target_comm:
				new_edge_list.append(j)

		beta = int(0.3*len(target_comm))
		intra_considered = []

		IG_edgeList_ = com_decept_safeness(target_comm, communities, IG_edgeList, beta, deg, out_deg, out_ratio, len(target_comm), new_adj, new_edge_list, intra_considered)

		g = igraph.Graph(directed = False)
		num_vertices = 62
		g.add_vertices(num_vertices)
		g.add_edges(IG_edgeList_)

		communities = g.community_walktrap().as_clustering()
		post_neighbours = check_neighbours(neighbours, communities)

		num_splits, comm_list = num_comm(target_comm, communities)
		sum_comm = sum_comm + num_splits

		nmi = igraph.compare_communities(comm_1, communities, method = "nmi")

		nmi_neighbourhood = igraph.compare_communities(pre_neighbours, post_neighbours, method = "nmi")

		entropy_val = get_entropy(comm_list)
		sum_entropy = sum_entropy + entropy_val

		NMI_List.append(nmi)
		Neighbourhood_NMI_List.append(nmi_neighbourhood)
		communities = safe_copy_comm