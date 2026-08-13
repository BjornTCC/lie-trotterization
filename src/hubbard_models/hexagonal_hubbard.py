import math
import warnings

import numpy as np
import scipy as sp
import networkx as nx

from scipy.optimize import fsolve, minimize_scalar
from src.hubbard_models.free_fermionic_errors import error_between_exp_of_free_fermionic
from src.hubbard_models._free_fermionic_computations import cast_data_to_array, spectral_norm_of_free_fermionic_operator, ff_commutator
from src.hubbard_models.split_operator_error_coefficients import fourth_order_augmented_split_operator_error_coefficients, fourth_order_suzuki_trotter_split_operator_error_coefficient
from src.bch_formula.commutator_bound import fourth_order_suzuki_trotter_commutator_bound

def compute_number_of_trotter_steps(t: float, eps: float, U: float, tau: float, Lx: int, Ly: int, type: str) -> int:
    if type == "tile":
        W = second_order_error_coefficient(U, tau, Lx, Ly)
        return math.ceil(
            np.sqrt(W / eps) * t**(3/2)
        )
    elif type == "tile 4th":
        W = _fourth_order_suzuki_trotter_error_coefficient(U, tau, Lx, Ly)
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
        W = _fourth_order_suzuki_trotter_error_coefficient(U, tau, Lx, Ly)
        return (eps / (5*W))**(1/4), math.ceil(
            4.463 * W ** (1 / 4) /(eps ** (5 / 4))
        )
    elif type == "augmented tile":
        return _augmented_plaquette_num_simulation_circuits(eps, U, tau, Lx, Ly, unitary_decomp = unitary_decomp)
    else:
        raise ValueError(f"Invalid argument for type: {type}. Must be one of: \'plaquette\', \'plaquette suzuki-trotter\' or \'plaquette augmented\'.")

def second_order_error_coefficient(U: float, tau: float, Lx: int, Ly: int) -> float:
    _, _, _, G = tile_trotterization_hexagonal(Lx, Ly)
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

def tile_trotterization_hexagonal(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    if  (Lx % 2) or (Ly % 2):
        raise ValueError(f"Lx and Ly must be even, got: Lx = {Lx}, Ly = {Ly}")
    whole = nx.hexagonal_lattice_graph(m=Ly, n=2 * Lx, periodic=True, with_positions=True)

    lx, ly = 2 * Lx, 2 * Ly

    Gr, Gb, Gy = nx.Graph(), nx.Graph(), nx.Graph()
    for g in [Gr, Gb, Gy]:
        g.add_nodes_from(whole.nodes)

    for j in range(lx):
        for i in range(Ly // 2):
            x = j % lx
            d = (x + 1) % 4
            y = (4 * i + d) % ly
            Gr.add_edges_from([
                ((x, y), (x, (y + 1) % ly)),
                ((x, y), (x, (y - 1) % ly))
            ])

    for j in range(Lx):
        for i in range(Ly // 2):
            x = (2 * j) % lx
            d = (x + 3) % 4
            y = (4 * i + d) % ly
            Gb.add_edges_from([
                ((x, y), (x, (y + 1) % ly)),
                ((x, y), ((x - 1) % lx, y))
            ])
            x = (2 * j - 1) % lx
            d = (x + 2) % 4
            y = (4 * i + d) % ly
            Gb.add_edges_from([
                ((x, y), (x, (y + 1) % ly)),
                ((x, y), ((x + 1) % lx, y))
            ])

    for j in range(Lx):
        for i in range(Ly // 2):
            x = (2 * j) % lx
            d = (x + 2) % 4
            y = (4 * i + d) % ly
            Gy.add_edges_from([
                ((x, y), (x, (y + 1) % ly)),
                ((x, y), ((x + 1) % lx, y))
            ])
            x = (2 * j + 1) % lx
            d = (2 * j) % 4
            y = (4 * i + d) % ly
            Gy.add_edges_from([
                ((x, y), (x, (y + 1) % ly)),
                ((x, y), ((x - 1) % lx, y))
            ])

    return Gr, Gb, Gy, whole

def colored_hexagonal_trotterization(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    whole = nx.hexagonal_lattice_graph(m=Ly, n=2*Lx, periodic=True, with_positions=True)

    lx, ly = 2*Lx, 2*Ly

    Gr, Gb, Gy = nx.Graph(), nx.Graph(), nx.Graph()
    for g in [Gr, Gb, Gy]:
        g.add_nodes_from(whole.nodes)

    for i in range(ly):
        for j in range(Lx):
            x = 2 * j + (i % 2)
            y = i
            Gr.add_edge(
                    (x,y),((x+1) % lx, y % ly)
            )

            Gb.add_edge(
                    (x,y),(x % lx, (y+1) % ly)
            )

            Gy.add_edge(
                    (x,y),(x % lx, (y-1) % ly)
            )

    return Gr, Gb, Gy, whole

def _fourth_order_suzuki_trotter_error_coefficient(U: float, tau: float, Lx: int, Ly: int) -> float:

    # This bound is likely loose due to the non-uniformity of the hopping decomposition

    four_terms = fourth_order_suzuki_trotter_commutator_bound(4)
    H1, H2, H3, _ = tile_trotterization_hexagonal(Lx,Ly)
    hoppings = [H1, H2, H3]
    N = 2*Lx*Ly
    d = 2
    res = 0
    for term, coeff in four_terms.items():
        if (term[0], term[1],term[2], term[3]) == (4,4,4,4):
            res += coeff * d * tau * U ** 4 * N
        elif (term[0], term[1], term[3]) == (4,4,4):
            res += coeff * 384 * d * tau * U ** 4 * N
        elif (term[0], term[2],term[3]) == (4,4,4):
            res += coeff * 128 * d * tau * U ** 4 * N
        elif (term[1],term[2], term[3]) == (4,4,4):
            res += coeff * 6 * d * tau * U ** 4 * N
        elif (term[0],term[1], term[2]) == (4,4,4):
            res += coeff * 128 * d * tau * U ** 4 * N
        elif (term[0], term[3]) == (4,4):
            res += coeff * 192 * d * tau * U ** 4 * N
        elif (term[1], term[3]) == (4,4):
            res += coeff * 240 * d * tau * U ** 4 * N
        elif (term[2], term[3]) == (4,4):
            res += coeff * 80 * d * tau * U ** 4 * N
        elif (term[1], term[2]) == (4,4):
            res += coeff * 128 * d * tau * U ** 4 * N
        elif (term[0], term[2]) == (4,4):
            res += coeff * 192 * d * tau * U ** 4 * N
        elif (term[0], term[1]) == (4,4):
            res += coeff * 64 * d * tau * U ** 4 * N
        elif term[0] == 4:
            res += coeff * 32 * d * tau * U ** 4 * N
        elif term[1] == 4:
            res += coeff * 48 * d * tau * U ** 4 * N
        elif term[2] == 4:
            res += coeff * 64 * d * tau * U ** 4 * N
        elif term[3] == 4:
            res += coeff * 96 * d * tau * U ** 4 * N

        else:
            operator = ff_commutator(
                hoppings[term[0] - 1],
                ff_commutator(
                    hoppings[term[1] - 1],
                    ff_commutator(
                        hoppings[term[2] - 1],
                        ff_commutator(
                            hoppings[term[3] - 1],
                            hoppings[term[4] - 1]
                        )
                    )
                )
            )
            res += coeff * tau**5 * spectral_norm_of_free_fermionic_operator(operator)

    return res

"""
Augmented free fermionic:
"""

def rbr(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    whole = nx.hexagonal_lattice_graph(m=Ly, n=2 * Lx, periodic=True, with_positions=True)

    lx, ly = 2 * Lx, 2 * Ly

    res = nx.Graph()
    res.add_nodes_from(whole.nodes)

    for i in range(ly):
        for j in range(Lx):
            x = 2 * j + (i % 2)
            y = i
            res.add_edge(
                (x, y), ((x + 2) % lx, (y - 1) % ly), weight=2 / 24
            )

    return res

def ryr(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    whole = nx.hexagonal_lattice_graph(m=Ly, n=2 * Lx, periodic=True, with_positions=True)

    lx, ly = 2 * Lx, 2 * Ly

    res = nx.Graph()
    res.add_nodes_from(whole.nodes)

    for i in range(ly):
        for j in range(Lx):
            x = 2 * j + (i % 2)
            y = i
            res.add_edge(
                (x, y), ((x + 2) % lx, (y + 1) % ly), weight=2 / 24
            )

    return res


def byb(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    whole = nx.hexagonal_lattice_graph(m=Ly, n=2 * Lx, periodic=True, with_positions=True)

    lx, ly = 2 * Lx, 2 * Ly

    res = nx.Graph()
    res.add_nodes_from(whole.nodes)

    for i in range(ly):
        for j in range(Lx):
            x = 2 * j + (i % 2)
            y = i
            res.add_edge(
                (x, y), (x, (y + 3) % ly), weight=2 / 24
            )

    return res


def bbr(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    whole = nx.hexagonal_lattice_graph(m=Ly, n=2 * Lx, periodic=True, with_positions=True)

    lx, ly = 2 * Lx, 2 * Ly

    res = nx.Graph()
    res.add_nodes_from(whole.nodes)

    for i in range(ly):
        for j in range(Lx):
            x = 2 * j + (i % 2)
            y = i
            res.add_edge(
                (x, y), ((x - 1) % lx, (y + 2) % ly), weight=-2 / 12
            )

    return res

def yyr(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    whole = nx.hexagonal_lattice_graph(m=Ly, n=2 * Lx, periodic=True, with_positions=True)

    lx, ly = 2 * Lx, 2 * Ly

    res = nx.Graph()
    res.add_nodes_from(whole.nodes)

    for i in range(ly):
        for j in range(Lx):
            x = 2 * j + (i % 2)
            y = i
            res.add_edge(
                (x, y), ((x - 1) % lx, (y - 2) % ly), weight=-2 / 12
            )

    return res


def yyb(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    whole = nx.hexagonal_lattice_graph(m=Ly, n=2 * Lx, periodic=True, with_positions=True)

    lx, ly = 2 * Lx, 2 * Ly

    res = nx.Graph()
    res.add_nodes_from(whole.nodes)

    for i in range(ly):
        for j in range(Lx):
            x = 2 * j + (i % 2)
            y = i
            res.add_edge(
                (x, y), (x, (y - 3) % ly), weight=-2 / 12
            )

    return res


def triple_mixed_graphs(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    whole = nx.hexagonal_lattice_graph(m=Ly, n=2 * Lx, periodic=True, with_positions=True)

    lx, ly = 2 * Lx, 2 * Ly

    res1, res2, res3 = nx.Graph(), nx.Graph(), nx.Graph()
    for res in [res1, res2, res3]:
        res.add_nodes_from(whole.nodes)

    for i in range(ly):
        for j in range(Lx):
            x = 2 * j + (i % 2)
            y = i
            res1.add_edge(
                (x, y), ((x + 1) % lx, (y - 2) % ly), weight=2 / 12
            )
            res2.add_edge(
                (x, y), ((x + 1) % lx, (y + 2) % ly), weight=2 / 12
            )
            res3.add_edge(
                (x, y), ((x - 1) % lx, y), weight=-4 / 12
            )

    return res1, res2, res3

def color_correctors(Lx, Ly) -> tuple[list[float], list[nx.Graph]]:
    _corr = [c(Lx, Ly) for c in [rbr, ryr, byb, bbr, yyr, yyb]] + list(triple_mixed_graphs(Lx, Ly))
    return [1/3,1/12,-1/6], _corr

def _compute_augmented_trotter_error_numerically(t: float, tau: float, Lx: int, Ly: int) -> float:
    Gr, Gb, Gy, G = colored_hexagonal_trotterization(Lx, Ly)
    Gs = [Gr, Gb, Gy]
    coeffs, correctors = color_correctors(Lx, Ly)
    seq1 = [
        1j*(tau*t)*0.5*(1+c*(tau*t)**2)*cast_data_to_array(g) for c,g in zip(coeffs, Gs)
    ]
    seq2 = [
        1j * (tau*t) ** 3 * cast_data_to_array(g) for g in correctors
    ]

    seq = seq1 + seq2 + seq1[::-1]
    target = 1j*(tau*t)*cast_data_to_array(G)
    return error_between_exp_of_free_fermionic(seq, target)