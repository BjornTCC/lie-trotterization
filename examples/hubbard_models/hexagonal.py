import networkx as nx
from openfermion import FermionOperator

from examples._primitives import hubbard_from_nx, spinless_PPP_model, ohno_potential_2d

def hexagonal_grid_hubbard_model(size: int | tuple[int], h: float = 1.0, U: float = 0.5) -> FermionOperator:
    if isinstance(size, int):
        size = (size, size)
    graph = nx.convert_node_labels_to_integers(nx.hexagonal_lattice_graph(*size))

    return hubbard_from_nx(h,U,graph)

def hexagonal_grid_spinless_PPP(size: int | tuple[int], h: float=1.0, U: float=1.5, spacing: float = 2.68) -> FermionOperator:
    if isinstance(size, int):
        size = (size, size)
    graph = nx.hexagonal_lattice_graph(*size, with_positions = True)
    vij = ohno_potential_2d(graph, U, scale=spacing)
    return spinless_PPP_model(h, graph, vij)