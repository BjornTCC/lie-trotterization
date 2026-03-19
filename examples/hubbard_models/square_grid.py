import networkx as nx
from openfermion import FermionOperator

from examples.hubbard_models._primitives import hubbard_from_nx

def square_hubbard_model(size: int, h: float, u: float) -> FermionOperator
    graph = nx.grid_2d_graph