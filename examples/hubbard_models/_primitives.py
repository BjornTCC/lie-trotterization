from openfermion import FermionOperator

import networkx as nx

def site_hopping(i: int, j: int) -> FermionOperator:
    return FermionOperator(((i,1),(j,0)), 1j) + FermionOperator(((j,1),(i,0)), 1j)

def coulomb_interaction(i: int, j: int) -> FermionOperator:
    return FermionOperator(((i,1),(i,0), (j,1),(j,0),), 1j)

def hubbard_from_nx(
        h: float,
        hopping_graph: nx.Graph,
        coulomb_graphs: list[tuple[nx.Graph, float]] = [],
) -> FermionOperator:
    terms = []
    for edge in hopping_graph.edges:
        terms.append(
            h*site_hopping(*edge)
        )

    for coulomb_graph, coeff in coulomb_graphs:
        for edge in coulomb_graph.edges:
            terms.append(
                coeff*coulomb_interaction(*edge)
            )
    return sum(terms)
