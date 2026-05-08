from openfermion import FermionOperator
from examples._primitives import hubbard_from_nx, spinless_PPP_model, ohno_potential_2d, coulomb_spinless_hamiltonian

import networkx as nx

def linear_graph(size: int, with_positions: bool = False, periodic: bool = True) -> nx.Graph:
    edges = [(i,i+1) for i in range(size - 1)]
    if periodic:
        edges += [(0, size-1)]
    g = nx.Graph()
    g.add_edges_from(edges)
    if with_positions:
        pos = {}
        for n in g.nodes:
            pos[n] = (float(n),0.0)
        nx.set_node_attributes(g, pos, name="pos")
    return g

def linear_hubbard_model(size: int, h: float = 1.0, U: float = 0.5) -> FermionOperator:
    graph = linear_graph(size)
    H = hubbard_from_nx(h, U, graph)
    return H

def linear_double(size: int, h: float = 1.0, u1: float = 0.5, u2: float = 0.25):
    graph1 = linear_graph(size)
    graph2 = nx.Graph()
    graph2.add_edges_from([(i,i+2) for i in range(size - 2)])
    H = coulomb_spinless_hamiltonian(h, graph1, [(graph1, u1),(graph2, u2)])
    return H

def linear_spinless_PPP(size: int, h: float = 1.0, U: float = 1.5, spacing: float = 2.0) -> FermionOperator:
    graph = linear_graph(size, with_positions=True)
    vij = ohno_potential_2d(graph, U, scale=spacing)
    return spinless_PPP_model(h, graph, vij)