import networkx as nx
from openfermion import FermionOperator

from examples.hubbard_models._primitives import hubbard_from_nx

def hexagonal_grid_hubbard_model(size: int | tuple[int], h: float = 1.0, u: float = 0.5) -> FermionOperator:
    if isinstance(size, int):
        size = (size, size)
    graph = nx.convert_node_labels_to_integers(nx.hexagonal_lattice_graph(*size))

    return hubbard_from_nx(h,graph, [(graph, u)])