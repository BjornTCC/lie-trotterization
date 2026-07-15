from copy import deepcopy
import networkx as nx

from src.hubbard_models._free_fermionic_computations import (
    ff_commutator,
    add_graphs,
    multiply_graph_by,
)

def _color_edge_decomp(G: nx.Graph, strategy: str = "DSATUR") -> list[nx.Graph]:
    if len(G.edges) == 0:
        return []
    LG = nx.line_graph(G)
    coloring = nx.greedy_color(LG, strategy = strategy)
    d = max(coloring.values()) + 1
    res = [nx.Graph() for _ in range(d)]
    for g in res:
        g.add_nodes_from(G.nodes)
    for edge, c in coloring.items():
        try:
            w = G[edge[0]][edge[1]]["weight"]
        except KeyError:
            w = 1.0
        res[c].add_edge(*edge)

        res[c][edge[0]][edge[1]]["weight"] = w
    return res

def corrector(terms: list[nx.Graph]) -> nx.Graph:
    res = nx.Graph()
    res.add_nodes_from(terms[0])
    N = len(terms)
    for j in range(N):
        for k in range(j+1, N):
            new_term = multiply_graph_by(
                ff_commutator(terms[j], ff_commutator(terms[k], terms[j])),
            1 / 24
            ) # the sign corrects for the 1j factor missing from the terms
            res = add_graphs(res, new_term)

    for k in range(N):
        for j in range(k+1, N):
            for i in range(k+1, N):
                new_term = multiply_graph_by(
                    ff_commutator(terms[i], ff_commutator(terms[j], terms[k])),
                    1 / 12
                )# the sign corrects for the 1j factor missing from the terms
                res = add_graphs(res, new_term)
    return res

def _recombine_corrector_into_terms(F: nx.Graph, terms: list[nx.Graph]) -> tuple[list[nx.Graph], nx.Graph]:
    new_terms = [deepcopy(g) for g in terms]
    new_F = deepcopy(F)
    for edge in F.edges:
        for T in new_terms:
            if edge in T.edges:
                T[edge[0]][edge[1]]["weight"] += F[edge[0]][edge[1]]["weight"]
                new_F.remove_edge(edge[0], edge[1])
                continue
    return new_terms, new_F

def second_order_s1_tiling(G: nx.Graph, decomp: list[nx.Graph] = None, strategy: str = "DSATUR", time: float = 1.0) -> list[nx.Graph]:
    if decomp is None:
        decomp = _color_edge_decomp(G, strategy)
    
    Left = [multiply_graph_by(g, 0.5 * time) for g in decomp[:-1]]
    Right = Left[::-1]
    Mid = [multiply_graph_by(decomp[-1], time)]
    return Left + Mid + Right
    
def augmented_s1_tiling(G: nx.Graph, decomp: list[nx.Graph] = None, strategy: str = "DSATUR", time: float = 1.0) -> list[nx.Graph]:
    if decomp is None:
        decomp = _color_edge_decomp(G, strategy)

    F = multiply_graph_by(corrector(decomp), time**3)

    scaled_decomp = [multiply_graph_by(g, time) for g in decomp]

    scaled_decomp, F = _recombine_corrector_into_terms(F, scaled_decomp)
    F_decomp = _color_edge_decomp(F, strategy)

    Left = [multiply_graph_by(g, 0.5) for g in scaled_decomp]
    Right = Left[::-1]
    return Left + F_decomp + Right