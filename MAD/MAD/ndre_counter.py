from math import log


def __inner_count(value):
    if not value:
        return 0
    return value * log(value, 2)


# --- pre-edit terms ---

def __count_pre_position_entropy(graph, add, out_distribute):
    total_degree = graph.ecount() + (1 if add else -1)
    pre_position_entropy = 0
    for node in graph.vs:
        var = out_distribute[node["name"]] / total_degree
        if var > 0:
            pre_position_entropy -= var * log(var, 2)
    return pre_position_entropy


def __count_pre_resistance(graph, partitions, partitions_inner_edges, partitions_out_volume, add):
    total_degree = graph.ecount() + (1 if add else -1)
    resistance = 0
    for index, part in enumerate(partitions):
        part_volume = partitions_out_volume[index]
        part_inner_edges = partitions_inner_edges[index]
        resistance -= (part_inner_edges / total_degree) * log(part_volume / total_degree, 2)
    return resistance


# --- deletion ---

def __count_position_entropy_by_pre_deletion(edge, pre_position_entropy, total_degree, outdegree_distribute):
    src = edge[0]
    src_now = (outdegree_distribute[src] - 1) / total_degree
    src_now = __inner_count(src_now)
    src_pre = outdegree_distribute[src] / total_degree
    src_pre = __inner_count(src_pre)
    return pre_position_entropy + src_pre - src_now


def __count_resistance_by_pre_deletion(pre_resistance, src_des, total_degree, partitions_inner_edges, partitions_out_volume):
    src, des = src_des
    if src == des:
        vs = partitions_out_volume[src]
        gs = vs - partitions_inner_edges[src]
        differ = 0 - (vs - gs - 1) / total_degree * log(((vs - 1) / total_degree), 2) + (vs - gs) / total_degree * log(
            vs / total_degree, 2
        )
    else:
        src_volume = partitions_out_volume[src]
        src_degree = partitions_inner_edges[src]
        src_pre = src_degree / total_degree * log(src_volume / total_degree, 2)
        src_now = src_degree / total_degree * log((src_volume - 1) / total_degree, 2)
        differ = src_pre - src_now
    resistance = pre_resistance + differ
    return resistance


# --- addition ---

def __count_position_entropy_by_pre_addition(edge, pre_position_entropy, total_degree, outdegree_distribute):
    src, des = edge
    src_now = (outdegree_distribute[src] + 1) / total_degree
    src_now = __inner_count(src_now)
    src_pre = outdegree_distribute[src] / total_degree
    src_pre = __inner_count(src_pre)
    return pre_position_entropy + src_pre - src_now


def __count_resistance_by_pre_addition(pre_resistance, src_des, total_degree, partitions_inner_edges, partitions_out_volume):
    src, des = src_des
    if src == des:
        vs = partitions_out_volume[src]
        gs = vs - partitions_inner_edges[src]
        differ = 0 - (vs - gs + 1) / total_degree * log(((vs + 1) / total_degree), 2) + (vs - gs) / total_degree * log(
            vs / total_degree, 2
        )
    else:
        src_volume = partitions_out_volume[src]
        src_degree = partitions_inner_edges[src]
        src_pre = src_degree / total_degree * log(src_volume / total_degree, 2)
        src_now = src_degree / total_degree * log((src_volume + 1) / total_degree, 2)
        differ = src_pre - src_now
    resistance = pre_resistance + differ
    return resistance


# --- public API ---

def count_pre_security_index(graph, partitions, partitions_inner_edges, partitions_out_volume, out_distribute, addition):
    pre_position_entropy = __count_pre_position_entropy(graph, addition, out_distribute)
    pre_resistance = __count_pre_resistance(graph, partitions, partitions_inner_edges, partitions_out_volume, addition)
    return pre_position_entropy, pre_resistance


def count_security_index_by_pre_deletion(
    pre_count, edge, src_des, total_degree, partitions_inner_edges, partitions_out_volume, outdegree_distribute
):
    pre_position_entropy, pre_resistance = pre_count
    position_entropy = __count_position_entropy_by_pre_deletion(edge, pre_position_entropy, total_degree, outdegree_distribute)
    resistance = __count_resistance_by_pre_deletion(pre_resistance, src_des, total_degree, partitions_inner_edges, partitions_out_volume)
    return resistance / position_entropy


def count_security_index_by_pre_addition(
    pre_count, edge, src_des, total_degree, partitions_inner_edges, partitions_out_volume, outdegree_distribute
):
    pre_position_entropy, pre_resistance = pre_count
    position_entropy = __count_position_entropy_by_pre_addition(edge, pre_position_entropy, total_degree, outdegree_distribute)
    resistance = __count_resistance_by_pre_addition(pre_resistance, src_des, total_degree, partitions_inner_edges, partitions_out_volume)
    return resistance / position_entropy
