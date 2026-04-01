from openfermion import FermionOperator
from openfermion.transforms import normal_ordered

from examples.hubbard_models._primitives import hubbard_from_nx, spinless_PPP_model, ohno_potential_2d

import networkx as nx

def linear_graph(size: int, with_positions: bool = False) -> nx.Graph:
    edges = [(i,i+1) for i in range(size - 1)]
    g = nx.Graph()
    g.add_edges_from(edges)
    if with_positions:
        pos = {}
        for n in g.nodes:
            pos[n] = (float(n),0.0)
        nx.set_node_attributes(g, pos, name="pos")
    return g

def linear_hubbard_model(size: int, h: float = 1.0, u: float = 0.5, normal: bool = True) -> FermionOperator:
    graph = linear_graph(size)
    H = hubbard_from_nx(h, graph, [(graph, u)])
    if normal:
        H = normal_ordered(H)
    return H

def linear_double(size: int, h: float = 1.0, u1: float = 0.5, u2: float = 0.25, normal: bool = True):
    graph1 = linear_graph(size)
    graph2 = nx.Graph()
    graph2.add_edges_from([(i,i+2) for i in range(size - 2)])
    H = hubbard_from_nx(h, graph1, [(graph1, u1),(graph2, u2)])
    if normal:
        H = normal_ordered(H)
    return H

def linear_spinless_PPP(size: int, h: float = 1.0, U = 1.5, spacing: float = 2.0) -> FermionOperator:
    graph = linear_graph(size, with_positions=True)
    vij = ohno_potential_2d(graph, U, scale=spacing)
    return spinless_PPP_model(h, graph, vij)