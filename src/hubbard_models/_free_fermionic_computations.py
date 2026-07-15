"""
Contains primitives for computing free-fermionic quantities using the fact that free fermionic systems are efficiently simulable.
"""

from copy import deepcopy

from functools import singledispatch
from typing import TypeVar
import warnings

import networkx as nx
import numpy as np
from openfermion.ops import FermionOperator
from openfermion.utils import commutator
from openfermion.transforms import normal_ordered

T = TypeVar('T', np.ndarray, FermionOperator, nx.Graph, nx.DiGraph)

# Warning messages

node_reordering_warning = "When computing free fermionic commutator performed between graphs with different nodes the ordering of nodes may change in the result."

# Utility functions

def spectral_norm_of_free_fermionic_operator(
        data: T
) -> float:
    casted_data = cast_data_to_array(data)
    return np.linalg.norm(casted_data, ord = 'nuc') / 2

@singledispatch
def ff_commutator(op1: T, op2: T) -> T:
    raise TypeError(f"Unsupported type: {type(op1)}")

@singledispatch
def cast_data_to_array(data: T) -> np.ndarray:
    """
    The base function. This executes if the input type
    doesn't match any of the registered types below.
    """
    raise TypeError(f"Unsupported type: {type(data)}")

def multiply_graph_by(G: nx.Graph, w: float) -> nx.Graph:
    res = deepcopy(G)
    for edge in res.edges:
        try:
            res[edge[0]][edge[1]]["weight"] *= w
        except KeyError:
            res[edge[0]][edge[1]]["weight"] = w
    return res

def add_graphs(G1: nx.Graph, G2: nx.Graph, scale_2: float = 1.0) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(G1.nodes)
    G.add_weighted_edges_from([(*e, G1[e[0]][e[1]]["weight"]) for e in G1.edges])
    for edge in G2.edges:
        if edge in G.edges:
            G[edge[0]][edge[1]]["weight"] += scale_2 * G2[edge[0]][edge[1]]["weight"]
        else:
            G.add_edge(*edge)
            G[edge[0]][edge[1]]["weight"] = scale_2 * G2[edge[0]][edge[1]]["weight"]
    return G


# Workhorse code

@cast_data_to_array.register(np.ndarray)
def _(data: np.ndarray) -> np.ndarray:
    # basic checks:

    if not (data.ndim == 2 and data.shape[0] == data.shape[1]):
        raise ValueError(f"np.ndarray has invalid shape: {data.shape}")

    return data

@cast_data_to_array.register(FermionOperator)
def _(data: FermionOperator) -> np.ndarray:

    # Find the maximum orbital index
    max_index = 0
    for term in data.terms:
        if len(term) != 2 or (not (term[0][1] ^ term[1][1])):
            raise ValueError(f"Operator is not free fermionic. Encountered non free-fermionic term: {term}")
        for ladder_operator in term:
            max_index = max(max_index, ladder_operator[0])

    num_orbitals = max_index + 1

    res = np.zeros((num_orbitals, num_orbitals), dtype = np.complex128)

    for term, coefficient in data.terms.items():
        i,j = (term[0][0], term[1][0]) if term[0][1] else (term[1][0], term[0][0])
        res[i,j] = coefficient

    return res

@cast_data_to_array.register(nx.Graph)
def _(data: nx.Graph) -> np.ndarray:
    return nx.adjacency_matrix(data, dtype = np.complex128).toarray()

@cast_data_to_array.register(nx.DiGraph)
def _(data: nx.DiGraph) -> np.ndarray:
    ad_mat = nx.adjacency_matrix(data, dtype = np.complex128).toarray()
    return ad_mat - np.conj(np.transpose(ad_mat))

@ff_commutator.register(np.ndarray)
def _(op1: np.ndarray, op2: np.ndarray) -> np.ndarray:
    if not isinstance(op2,  np.ndarray):
        raise TypeError(f"Both arguments must be np.ndarray.")
    return op1 @ op2 - op2 @ op1

@ff_commutator.register(FermionOperator)
def _(op1: FermionOperator, op2: FermionOperator) -> FermionOperator:
    if not isinstance(op2,  FermionOperator):
        raise TypeError(f"Both arguments must be FermionOperator.")
    return normal_ordered(commutator(op1, op2))

@ff_commutator.register(nx.Graph)
def _(op1: nx.Graph, op2: nx.Graph | nx. DiGraph) -> nx.Graph:

    if isinstance(op2, nx.DiGraph):
        return _commutator_of_graph_and_digraph(op1,op2)
    elif isinstance(op2, nx.Graph):
        return _commutator_of_graphs(op1, op2)
    raise TypeError(f"Both arguments must be networkx Graph or DiGraph.")

@ff_commutator.register(nx.DiGraph)
def _(op1: nx.DiGraph, op2: nx.Graph | nx. DiGraph) -> nx.DiGraph | nx.Graph:
    if isinstance(op2, nx.DiGraph):
        return _commutator_of_digraphs(op1, op2)
    elif isinstance(op2, nx.Graph):
        res = _commutator_of_graph_and_digraph(op2, op1)
        _flip_edge_weights(res)
        return res
    raise TypeError(f"Both arguments must be networkx Graph or DiGraph.")

def _commutator_of_graphs(G1: nx.Graph, G2: nx.Graph) -> nx.DiGraph:
    res = nx.DiGraph()
    if G1.nodes == G2.nodes:
        nodes = G1.nodes # This keeps the node ordering of G1
    else:
        nodes = set(G1.nodes).union(G2.nodes)
        warnings.warn(node_reordering_warning)
    res.add_nodes_from(nodes)
    if nx.is_weighted(G1):
        edge_weights_1 = {edge: G1[edge[0]][edge[1]]["weight"] for edge in G1.edges}
    else:
        edge_weights_1 = {edge: 1 for edge in G1.edges}

    if nx.is_weighted(G2):
        edge_weights_2 = {edge: G2[edge[0]][edge[1]]["weight"] for edge in G2.edges}
    else:
        edge_weights_2 = {edge: 1 for edge in G2.edges}

    edge_weights_1.update({x[::-1]: y for x,y in edge_weights_1.items()})
    edge_weights_2.update({x[::-1]: y for x,y in edge_weights_2.items()})

    for edge in G1.edges:
        n1 = edge[0]
        n2 = edge[1]
        ad_nodes1 = G2.neighbors(n1)
        for other_node in ad_nodes1:
            w = edge_weights_1[edge] * edge_weights_2[(n1, other_node)]
            if res.has_edge(n2, other_node):
                res[n2][other_node]["weight"] += w
            elif res.has_edge(other_node, n2):
                res[other_node][n2]["weight"] -= w
            else:
                res.add_edge(n2, other_node, weight = w)

        ad_nodes2 = G2.neighbors(n2)
        for other_node in ad_nodes2:
            w = edge_weights_1[edge] * edge_weights_2[(n2, other_node)]
            if res.has_edge(n1, other_node):
                res[n1][other_node]["weight"] += w
            elif res.has_edge(other_node, n1):
                res[other_node][n1]["weight"] -= w
            else:
                res.add_edge(n1, other_node, weight = w)

    _remove_zero_weight_edges(res, copy = False)
    return res


def _commutator_of_graph_and_digraph(G1: nx.graph, G2: nx.DiGraph) -> nx.Graph:
    res = nx.Graph()
    if G1.nodes == G2.nodes:
        nodes = G1.nodes # This keeps the node ordering of G1
    else:
        nodes = set(G1.nodes).union(G2.nodes)
        warnings.warn(node_reordering_warning)
    res.add_nodes_from(nodes)
    G2_rev = G2.reverse()

    if nx.is_weighted(G1):
        edge_weights_1 = {edge: G1[edge[0]][edge[1]]["weight"] for edge in G1.edges}
    else:
        edge_weights_1 = {edge: 1 for edge in G1.edges}

    if nx.is_weighted(G2):
        edge_weights_2 = {edge: G2[edge[0]][edge[1]]["weight"] for edge in G2.edges}
    else:
        edge_weights_2 = {edge: 1 for edge in G2.edges}

    edge_weights_1.update({x[::-1]: y for x, y in edge_weights_1.items()})
    edge_weights_2.update({x[::-1]: -y for x, y in edge_weights_2.items()})

    for edge in G1.edges:
        n1 = edge[0]
        n2 = edge[1]

        ad_nodes1 = G2.neighbors(n1)
        for other_node in ad_nodes1:
            w = edge_weights_1[edge] * edge_weights_2[(n1, other_node)]
            if res.has_edge(n2, other_node):
                res[n2][other_node]['weight'] += w

            else:
                res.add_edge(n2, other_node, weight=w)

        ad_nodes2 = G2.neighbors(n2)
        for other_node in ad_nodes2:
            w = edge_weights_1[edge] * edge_weights_2[(n2, other_node)]
            if res.has_edge(n1, other_node):
                res[n1][other_node]['weight'] += w

            else:
                res.add_edge(n1, other_node, weight=w)

        ad_nodes3 = G2_rev.neighbors(n1)
        for other_node in ad_nodes3:
            w = edge_weights_1[edge] * edge_weights_2[(other_node, n1)]
            if res.has_edge(n2, other_node):
                res[n2][other_node]['weight'] += -w

            else:
                res.add_edge(n2, other_node, weight=-w)

        ad_nodes4 = G2_rev.neighbors(n2)
        for other_node in ad_nodes4:
            w = edge_weights_1[edge] * edge_weights_2[(other_node, n2)]
            if res.has_edge(n1, other_node):
                res[n1][other_node]['weight'] += -w

            else:
                res.add_edge(n1, other_node, weight=-w)

    _remove_zero_weight_edges(res, copy = False)

    return res

def _commutator_of_digraphs(G1: nx.DiGraph, G2: nx.DiGraph) -> nx.DiGraph:
    res = nx.DiGraph()
    if G1.nodes == G2.nodes:
        nodes = G1.nodes # This keeps the node ordering of G1
    else:
        nodes = set(G1.nodes).union(G2.nodes)
        warnings.warn(node_reordering_warning)
    res.add_nodes_from(nodes)
    G2_rev = G2.reverse()

    if nx.is_weighted(G1):
        edge_weights_1 = {edge: G1[edge[0]][edge[1]][["weight"]] for edge in G1.edges}
    else:
        edge_weights_1 = {edge: 1 for edge in G1.edges}

    if nx.is_weighted(G2):
        edge_weights_2 = {edge: G2[edge[0]][edge[1]][["weight"]] for edge in G2.edges}
    else:
        edge_weights_2 = {edge: 1 for edge in G2.edges}

    edge_weights_1.update({x[::-1]: y for x, y in edge_weights_1.items()})
    edge_weights_2.update({x[::-1]: -y for x, y in edge_weights_2.items()})

    for edge in G1.edges:
        n1 = edge[0]
        n2 = edge[1]

        ad_nodes1 = G2.neighbors(n1)
        for other_node in ad_nodes1:
            w = edge_weights_1[edge] * edge_weights_2[(n1, other_node)]
            if res.has_edge(n2, other_node):
                res[n2][other_node]['weight'] += -w

            else:
                res.add_edge(n2, other_node, weight=-w)

        ad_nodes2 = G2.neighbors(n2)
        for other_node in ad_nodes2:
            w = edge_weights_1[edge] * edge_weights_2[(n2, other_node)]
            if res.has_edge(n1, other_node):
                res[n1][other_node]['weight'] += w

            else:
                res.add_edge(n1, other_node, weight=w)

        ad_nodes3 = G2_rev.neighbors(n1)
        for other_node in ad_nodes3:
            w = edge_weights_1[edge] * edge_weights_2[(other_node, n1)]
            if res.has_edge(n2, other_node):
                res[n2][other_node]['weight'] += w

            else:
                res.add_edge(n2, other_node, weight=w)

        ad_nodes4 = G2_rev.neighbors(n2)
        for other_node in ad_nodes4:
            w = edge_weights_1[edge] * edge_weights_2[(other_node, n2)]
            if res.has_edge(n1, other_node):
                res[n1][other_node]['weight'] += -w

            else:
                res.add_edge(n1, other_node, weight=-w)

    _remove_zero_weight_edges(res, copy=False)

    return res

def _remove_zero_weight_edges(G: nx.Graph, copy: bool = False) -> nx.Graph:
    graph_to_modify = G.copy() if copy else G

    zero_weight_edges = [
        (u, v) for u, v, weight in graph_to_modify.edges(data='weight')
        if weight == 0
    ]

    graph_to_modify.remove_edges_from(zero_weight_edges)
    if copy:
        return graph_to_modify

def _flip_edge_weights(G: nx.Graph, copy: bool = False) -> nx.Graph:
    graph_to_modify = G.copy() if copy else G
    if not nx.is_weighted(graph_to_modify):
        for edge in graph_to_modify.edges:
            graph_to_modify[edge[0]][edge[1]]["weight"] = -1
    else:
        for edge in graph_to_modify.edges:
            graph_to_modify[edge[0]][edge[1]]["weight"] *= -1

    if copy:
        return graph_to_modify
