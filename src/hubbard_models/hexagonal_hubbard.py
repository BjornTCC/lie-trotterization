import math
import warnings

import numpy as np
import scipy as sp
import networkx as nx

from scipy.optimize import fsolve, minimize_scalar
from src.hubbard_models.free_fermionic_errors import error_between_exp_of_free_fermionic
from src.hubbard_models._free_fermionic_computations import cast_data_to_array

def s3_tile_trotterization_graphs(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    if (Ly % 2):
        raise ValueError(f"Only even values of Ly are allowed, got: {Ly}")
    whole = nx.hexagonal_lattice_graph(Lx, Ly, periodic=True)
    pos = nx.get_node_attributes(whole, 'pos')

    Gr, Gb, Gy = nx.Graph(), nx.Graph(), nx.Graph()
    for g in [Gr, Gb, Gy]:
        g.add_nodes_from(whole.nodes)
        for node, p in pos.items():
            g.nodes[node]['pos'] = p

    edge_groups = [[], [], []]

    for i in range(Ly + 1):
        d = i % 2
        for j in range(Lx + 2):
            x = i
            y = 2 * j + (i + 1) % 2
            edge_groups[d].extend([
                ((x, y), (x, (y + 1) % (2 * Lx))),
                ((x, y), (x, (y - 1) % (2 * Lx))),
                ((x, y), ((x - 1) % Ly, y))
            ])
            d = (d + 1) % 3

    Gr.add_edges_from([x for x in edge_groups[0] if (x[0] in whole.nodes) and (x[1] in whole.nodes)])
    Gb.add_edges_from([x for x in edge_groups[1] if (x[0] in whole.nodes) and (x[1] in whole.nodes)])
    Gy.add_edges_from([x for x in edge_groups[2] if (x[0] in whole.nodes) and (x[1] in whole.nodes)])
    return Gr, Gb, Gy, whole

def hexagonal_3_correctors(Lx: int, Ly: int) -> tuple[nx.DiGraph, ...]:
    if (Ly % 2):
        raise ValueError(f"Only even values of Ly are allowed, got: {Ly}")
    whole = nx.hexagonal_lattice_graph(Lx, Ly, periodic=True)
    pos = nx.get_node_attributes(whole, 'pos')

    Gbr, Gyr, Gyb = nx.DiGraph(), nx.DiGraph(), nx.DiGraph()
    for g in [Gbr, Gyr, Gyb]:
        g.add_nodes_from(whole.nodes)
        for node, p in pos.items():
            g.nodes[node]['pos'] = p

    edge_groups = [[], [], []]

    for i in range(Ly + 1):
        d = i % 2
        for j in range(Lx + 2):
            x = i
            y = 2 * j + (i + 1) % 2
            edge_groups[d].extend([
                ((x, (y + 2) % (2 * Lx)), (x, y)),
                (((x + 1) % Ly, (y - 1) % (2 * Lx)), (x, y)),
                (((x - 1) % Ly, (y - 1) % (2 * Lx)), (x, y))
            ])
            d = (d + 1) % 3

    Gbr.add_edges_from([x for x in edge_groups[0] if (x[0] in whole.nodes) and (x[1] in whole.nodes)])
    Gyr.add_edges_from([x[::-1] for x in edge_groups[2] if (x[0] in whole.nodes) and (x[1] in whole.nodes)])
    Gyb.add_edges_from([x for x in edge_groups[1] if (x[0] in whole.nodes) and (x[1] in whole.nodes)])
    return Gbr, Gyr, Gyb
def hexagonal_3_correctors_decomposed(Lx: int, Ly: int) -> tuple[list[nx.DiGraph, ...], ...]:
    if (Ly % 2):
        raise ValueError(f"Only even values of Ly are allowed, got: {Ly}")
    whole = nx.hexagonal_lattice_graph(Lx, Ly, periodic=True)
    pos = nx.get_node_attributes(whole, 'pos')

    Gbr, Gyr, Gyb = [nx.DiGraph(), nx.DiGraph(), nx.DiGraph()], [nx.DiGraph(), nx.DiGraph(), nx.DiGraph()],[nx.DiGraph(), nx.DiGraph(), nx.DiGraph()]
    for g in Gbr + Gyr + Gyb:
        g.add_nodes_from(whole.nodes)
        for node, p in pos.items():
            g.nodes[node]['pos'] = p

    edge_groups = [[[], [], []], [[],[],[]], [[],[],[]]]

    for i in range(Ly + 1):
        d = i % 2
        for j in range(Lx + 2):
            x = i
            y = 2 * j + (i + 1) % 2
            edge_groups[(i + d) % 3][d].extend([
                ((x, (y + 2) % (2 * Lx)), (x, y)),
                (((x + 1) % Ly, (y - 1) % (2 * Lx)), (x, y)),
                (((x - 1) % Ly, (y - 1) % (2 * Lx)), (x, y))
            ])
            d = (d + 1) % 3

    for i in range(3):
        Gbr[i].add_edges_from([x for x in edge_groups[i][0] if (x[0] in whole.nodes) and (x[1] in whole.nodes)])
        Gyr[i].add_edges_from([x[::-1] for x in edge_groups[i][2] if (x[0] in whole.nodes) and (x[1] in whole.nodes)])
        Gyb[i].add_edges_from([x for x in edge_groups[i][1] if (x[0] in whole.nodes) and (x[1] in whole.nodes)])
    return Gbr, Gyr, Gyb

def _compute_augmented_trotter_error_numerically(t: float, tau: float, Lx: int, Ly: int) -> float:
    _Gr, _Gb, _Gy, _whole = s3_tile_trotterization_graphs(Lx, Ly)
    _Gbr,_Gyr,_Gyb = hexagonal_3_correctors_decomposed(Lx, Ly)

    Gr, Gb, Gy, whole = cast_data_to_array(_Gr), cast_data_to_array(_Gb), cast_data_to_array(_Gy), cast_data_to_array(_whole)
    Gbr = [
        -(t*tau)**2 * cast_data_to_array(g) / 24 for g in _Gbr
    ]
    Gyr = [
        -(t*tau)**2 * cast_data_to_array(g) / 24 for g in _Gyr
    ]
    Gyb = [
        -(t*tau)**2 * cast_data_to_array(g) / 24 for g in _Gyb
    ]
    MGbr, MGyr, MGyb = [-g for g in Gbr[::-1]],[-g for g in Gyr[::-1]],[-g for g in Gyb[::-1]],

    sequence = [
        *Gbr,
        *Gyr,
        1j*(t*tau)*(1/2)*Gr,
        *Gyr,
        *Gbr,
        *Gyb,
        1j*(tau)*t*(1/2)*Gb,
        *Gyb,
        1j*(t*tau)*Gy,
        *MGyb,
        1j*(t*tau)*(1/2)*Gb,
        *MGyb,
        *MGbr,
        *MGyr,
        1j*(t*tau)*(1/2)*Gr,
        *MGyr,
        *MGbr
    ]

    target = 1j*t*tau*whole
    return error_between_exp_of_free_fermionic(sequence, target)