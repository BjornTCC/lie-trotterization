import math
import warnings

import numpy as np
import scipy as sp
import networkx as nx

from scipy.optimize import fsolve, minimize_scalar
from src.hubbard_models.free_fermionic_errors import error_between_exp_of_free_fermionic
from src.hubbard_models._free_fermionic_computations import cast_data_to_array, spectral_norm_of_free_fermionic_operator
from src.hubbard_models.split_operator_error_coefficients import fourth_order_augmented_split_operator_error_coefficients, fourth_order_suzuki_trotter_split_operator_error_coefficient

def compute_number_of_trotter_steps(t: float, eps: float, U: float, tau: float, Lx: int, Ly: int, type: str) -> int:
    if type == "tile":
        W = second_order_error_coefficient(U, tau, Lx, Ly)
        return math.ceil(
            np.sqrt(W / eps) * t**(3/2)
        )
    elif type == "tile 4th":
        W = fourth_order_suzuki_trotter_split_operator_error_coefficient(U, tau, 2*Lx*Ly,3)
        #W = _fourth_order_suzuki_trotter_error_coefficient(U, tau, L)
        return math.ceil(
            (W / eps)**(1/4) * t**(5/4)
        )
    elif type == "augmented tile":
        return _augmented_hexagonal_trotter_steps(t, eps, U, tau, Lx, Ly)
    else:
        raise ValueError(f"Invalid argument for type: {type}. Must be one of: \'plaquette\', \'plaquette suzuki-trotter\' or \'plaquette augmented\'.")

def compute_evolution_time_and_number_of_simulation_circuits_for_qpe(eps: float, U: float, tau: float, Lx: int, Ly: int, type: str, unitary_decomp: bool = True) -> tuple[float,int]:
    if type == "tile":
        W = second_order_error_coefficient(U, tau, Lx, Ly)
        return np.sqrt(eps / (3*W)), math.ceil(
            6.203 * W**(1/2) / (eps**(3/2))
        )
    elif type == "tile 4th":
        W = fourth_order_suzuki_trotter_split_operator_error_coefficient(U, tau, 2*Lx*Ly, 3)
       # W = _fourth_order_suzuki_trotter_error_coefficient(U, tau, L)
        return (eps / (5*W))**(1/4), math.ceil(
            4.463 * W ** (1 / 4) /(eps ** (5 / 4))
        )
    elif type == "augmented tile":
        return _augmented_plaquette_num_simulation_circuits(eps, U, tau, Lx, Ly, unitary_decomp = unitary_decomp)
    else:
        raise ValueError(f"Invalid argument for type: {type}. Must be one of: \'plaquette\', \'plaquette suzuki-trotter\' or \'plaquette augmented\'.")

def second_order_error_coefficient(U: float, tau: float, Lx: int, Ly: int) -> float:
    _, _, _, G = s3_tile_trotterization_graphs(Lx, Ly)
    R1 = 2*spectral_norm_of_free_fermionic_operator(G)

    return U*U*tau * R1 / 24 + 9.9*U*tau**2 * (2*Lx*Ly) / 12 + 0.8532 * tau**3 * (2*Lx*Ly)

def _augmented_hexagonal_trotter_steps(t: float, eps: float, U: float, tau: float, Lx: int, Ly: int) -> int:
    W5, W6, W7 = fourth_order_augmented_split_operator_error_coefficients(U, tau, 2*Lx*Ly, 3, unitary_decomp=True)

    f = lambda n: W5 * t**5 / n**4 + W6 * t**6 / n**6 + W7 * t**7 / n**6 + n * _compute_augmented_trotter_error_numerically(t/n, tau, Lx, Ly) - eps

    x0 = W5**(1/4) * t**(5/4) /(eps**(1/4))
    res = fsolve(f, x0, full_output = True)
    root = res[0][0]

    return math.ceil(root)

def _augmented_plaquette_num_simulation_circuits(eps: float, U: float, tau: float, Lx: int, Ly: int, unitary_decomp: bool = True) -> tuple[float, int]:
    W5, W6, W7 = fourth_order_augmented_split_operator_error_coefficients(U, tau, 2*Lx*Ly, 3, unitary_decomp=unitary_decomp)
    f = lambda t: W5 * t**5 + W6 * t**6 + W7 * t**7 + _compute_augmented_trotter_error_numerically(t, tau, Lx, Ly)

    opt_func = lambda t: 0.76*np.pi/(t*eps - f(t))

    t0 = (eps / (5*W5))**(1/4)
    tm = (eps / W5)**(1/4)

    bnds = [(0.00001,tm)]# These bounds ensure Npe > 0

    min_res = minimize_scalar(opt_func, bounds = bnds[0])

    if not min_res.success:
        warnings.warn(f"Npe optimizer failed with message: {min_res.message}")

    Npe = min_res.fun

    return  min_res.x, math.ceil(Npe)

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