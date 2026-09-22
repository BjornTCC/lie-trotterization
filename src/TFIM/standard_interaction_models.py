import networkx as nx
import numpy as np

def nearest_neighbor_from_graph(G: nx.Graph, interaction_strength: float) -> tuple[callable, np.ndarray]:
    def f(n1, n2):
        if ((n1, n2) in G.edges) or ((n2, n1) in G.edges):
            return interaction_strength/2
        return 0.0

    return f, np.array(G.nodes)

def hexagonal_power_law(Lx: int, Ly: int, interaction_strength: float, alpha: float, max_dist: float = None) -> callable:
    G = nx.hexagonal_lattice_graph(m=Ly, n=2*Lx, periodic=True, with_positions=True)
    pos = nx.get_node_attributes(G,'pos')
    positions = [pos[(Lx, Ly)]]
    for n in G.nodes:
        if n != (Lx, Ly):
            positions.append(pos[n])
    positions = 2*np.array(positions) # Factor 2 ensures that all distances >= 1

    def poly(n1, n2):
        dist = np.sqrt(np.sum((n1 - n2) ** 2))
        if max_dist is not None and dist > max_dist:
            return 0.0

        return 0.5*interaction_strength / dist ** alpha

    return poly, positions

def square_power_law(Lx: int, Ly: int, interaction_strength: float, alpha: float, max_dist: float = None) -> callable:
    assert (Lx % 2) and (Ly % 2)
    positions = [[0,0]]
    for i in range(-(Lx // 2), Lx // 2 + 1):
        for j in range(-(Lx // 2), Ly // 2 + 1):
            if [i,j] != [0,0]:
                positions.append([i,j])

    def poly(n1, n2):
        dist = np.sqrt(np.sum((n1 - n2) ** 2))
        if max_dist is not None and dist > max_dist:
            return 0.0

        return 0.5*interaction_strength / dist ** alpha

    positions = np.array(positions)
    return poly, positions

def cubic_power_law(Lx: int, Ly: int, Lz: int, interaction_strength: float, alpha: float, max_dist: float = None) -> callable:
    assert (Lx % 2) and (Ly % 2) and (Lz % 2)
    positions = [[0, 0, 0]]

    for i in range(-(Lx // 2), Lx // 2 + 1):
        for j in range(-(Ly // 2), Ly // 2 + 1):
            for k in range(-(Lz // 2), Lz // 2 + 1):
                if [i,j,k] != [0,0,0]:
                    positions.append([i,j,k])

    def poly(n1, n2):
        dist = np.sqrt(np.sum((n1 - n2) ** 2))
        if max_dist is not None and dist > max_dist:
            return 0.0

        return 0.5*interaction_strength / dist ** alpha

    positions = np.array(positions)
    return poly, positions