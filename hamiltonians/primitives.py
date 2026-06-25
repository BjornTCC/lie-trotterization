import numpy as np

from openfermion import FermionOperator
from openfermion.transforms import normal_ordered

import networkx as nx

def site_hopping(i: int, j: int) -> FermionOperator:
    return FermionOperator(((i,1),(j,0)), 1j) + FermionOperator(((j,1),(i,0)), 1j)

def hubbard_from_nx(
    h: float,
    u: float,
    hopping_graph: nx.Graph,
) -> FermionOperator:

    N = len(hopping_graph.nodes)

    spin_orbital_to_index = {
        (i,s): s*N + i for i in hopping_graph.nodes for s in [0,1]
    }

    terms = []

    for edge in hopping_graph.edges:
        for s in [0,1]:
            terms.append(
                h*site_hopping(*[spin_orbital_to_index[(e,s)] for e in edge])
            )

    for n in hopping_graph.nodes:
        terms.append(u*coulomb_interaction(
            spin_orbital_to_index[(n,0)], spin_orbital_to_index[(n,1)]
        ))

    return normal_ordered(sum(terms))

def coulomb_spinless_hamiltonian(
        h: float | None,
        hopping_graph: nx.Graph,
        coulomb_graphs: list[tuple[nx.Graph, float]] = [],
) -> FermionOperator:
    terms = []
    for edge in hopping_graph.edges:
        if h is None:
            terms.append(
                np.random.uniform()*site_hopping(*edge)
            )
        else:
            terms.append(
                h * site_hopping(*edge)
            )

    for coulomb_graph, coeff in coulomb_graphs:
        for edge in coulomb_graph.edges:
            terms.append(
                coeff*coulomb_interaction(*edge)
            )
    return normal_ordered(sum(terms))

def coulomb_interaction(i: int, j: int) -> FermionOperator:
    return FermionOperator(((i,1),(i,0), (j,1),(j,0),), 1j)

def spinless_PPP_model(
        h: float,
        hopping_graph: nx.Graph,
        vij: callable,
) -> FermionOperator:
    terms = []
    graph_index_map = {n: ind for ind,n in enumerate(hopping_graph.nodes)}
    for edge in hopping_graph.edges:
        terms.append(
            h * site_hopping(graph_index_map[edge[0]],graph_index_map[edge[1]])
        )

    N = len(hopping_graph.nodes)
    for i in range(N):
        for j in range(i+1, N):
            nodei, nodej = list(hopping_graph.nodes)[i], list(hopping_graph.nodes)[j]

            coeff = vij(nodei, nodej)

            terms.append(
                coeff * coulomb_interaction(j,i)
            )

    return normal_ordered(sum(terms))

def ohno_potential_2d(
        graph: nx.Graph,
        U: float,
        scale: float = 2.68 # Distance to scale the interatomic distances by
) -> callable:
    positions = nx.get_node_attributes(graph, 'pos')
    def res(node1, node2) -> float:
        pos1, pos2 = positions[node1], positions[node2]
        rij = scale*np.sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
        return U /(np.sqrt(1 + (U*rij/14.397)**2))
    return res