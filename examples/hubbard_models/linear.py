from openfermion import FermionOperator
from openfermion.transforms import normal_ordered

from examples.hubbard_models._primitives import hubbard_from_nx

import networkx as nx

def linear_graph(size: int) -> nx.Graph:
    edges = [(i,i+1) for i in range(size - 1)]
    g = nx.Graph()
    g.add_edges_from(edges)
    return g

def linear_hubbard_model(size: int, h: float = 1.0, u: float = 0.5, normal: bool = True) -> FermionOperator:
    graph = linear_graph(size)
    H = hubbard_from_nx(h, graph, [(graph, u)])
    if normal:
        H = normal_ordered(H)
    return H

def linear_double(size, h: float = 1.0, u1: float = 0.5, u2: float = 0.25, normal: bool = True):
    graph1 = linear_graph(size)
    graph2 = nx.Graph()
    graph2.add_edges_from([(i,i+2) for i in range(size - 2)])
    H = hubbard_from_nx(h, graph1, [(graph1, u1),(graph2, u2)])
    if normal:
        H = normal_ordered(H)
    return H
