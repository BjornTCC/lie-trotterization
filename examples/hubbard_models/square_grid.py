import networkx as nx
from openfermion import FermionOperator

from examples.hubbard_models._primitives import hubbard_from_nx, spinless_PPP_model, ohno_potential_2d

def square_grid_hubbard_model(size: int | tuple[int], h: float = 1.0, u: float = 0.5) -> FermionOperator:
    if isinstance(size, int):
        size = (size, size)
    graph = nx.convert_node_labels_to_integers(nx.grid_2d_graph(*size))

    return hubbard_from_nx(h,graph, [(graph, u)])

def square_grid_spinless_PPP(size: int | tuple[int], h: float=1.0, U: float=1.5, spacing: float = 2.68) -> FermionOperator:
    if isinstance(size, int):
        size = (size, size)
    graph = nx.grid_2d_graph(*size)

    pos = {}
    for n in graph.nodes:
        pos[n] = (float(n[0]), float(n[1]))
    nx.set_node_attributes(graph, pos, name="pos")

    vij = ohno_potential_2d(graph, U, scale=spacing)
    return spinless_PPP_model(h, graph, vij)