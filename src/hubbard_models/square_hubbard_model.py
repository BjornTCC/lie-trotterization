"""
Functions for generating circuits to simulate hubbard models on the 2D square lattice with fully symmetric coulomb and hopping hamiltonians.

Second order circuits are copied directly from:
    Earl T Campbell 2022 Quantum Sci. Technol. 7 015007

In all calculations we assume a periodic lattice
"""

import math
import warnings

import numpy as np
import scipy as sp
import networkx as nx

from openfermion.ops import FermionOperator

from hamiltonians.primitives import hubbard_from_nx

from src.hubbard_models._free_fermionic_computations import spectral_norm_of_free_fermionic_operator, ff_commutator, cast_data_to_array
from src.hubbard_models.split_operator_error_coefficients import (
    plaquette_second_order_split_operator_error_coefficient,
    fourth_order_suzuki_trotter_split_operator_error_coefficient,
    fourth_order_augmented_split_operator_error_coefficients,
    square_augmented_so_coeffs
)
from src.hubbard_models.free_fermionic_errors import error_between_exp_of_free_fermionic
from src.bch_formula.commutator_bound import augmented_commutator_bound
from src.tools.single_variable_optimizer import one_d_optimizer
from src.tools.add_dict import add_dicts

from scipy.optimize import fsolve

def get_fermionic_operator(U: float, tau: float, L: int) -> FermionOperator:
    graph = nx.convert_node_labels_to_integers(nx.grid_2d_graph(L,L, periodic = True))
    return hubbard_from_nx(tau, U, graph)

def compute_number_of_trotter_steps_kwargs(t: float, eps: float, U: float, tau: float, L: int, type: str) -> dict:
    if type == "plaquette":
        W = second_order_error_coefficient_alt(U, tau, L)
        return {"num_trotter_steps": math.ceil(
            np.sqrt(W / eps) * t**(3/2)
        )}
    elif type == "plaquette suzuki-trotter":
        W = _fourth_order_suzuki_trotter_error_coefficient(U, tau, L)
        return {"num_trotter_steps": math.ceil(
            (W / eps)**(1/4) * t**(5/4)
        )}
    elif type == "augmented plaquette":
        return {"error_coefficients": _augmented_plaquette_error_coefficients(U, tau, L, unitary_decomp=True)}
    else:
        raise ValueError(f"Invalid argument for type: {type}. Must be one of: \'plaquette\', \'plaquette suzuki-trotter\' or \'plaquette augmented\'.")

def compute_evolution_time_and_number_of_simulation_circuits_for_qpe_kwargs(eps: float, U: float, tau: float, L: int, type: str, unitary_decomp: bool = True) -> dict:
    if type == "plaquette":
        W = second_order_error_coefficient_alt(U, tau, L)
        return {"num_simulation_steps": math.ceil(
            6.203 * W**(1/2) / (eps**(3/2))
        )}
    elif type == "plaquette suzuki-trotter":
        W = _fourth_order_suzuki_trotter_error_coefficient(U, tau, L)
        return {"num_simulation_steps": math.ceil(
            4.463 * W ** (1 / 4) /(eps ** (5 / 4))
        )}
    elif type == "augmented plaquette":
        return {"unitary_error_coefficients": _augmented_plaquette_error_coefficients(U, tau, L, unitary_decomp=unitary_decomp)}
    else:
        raise ValueError(f"Invalid argument for type: {type}. Must be one of: \'plaquette\', \'plaquette suzuki-trotter\' or \'plaquette augmented\'.")

def _augmented_plaquette_error_coefficients(U: float, tau: float, L: int, unitary_decomp: bool) -> tuple[float,...]:
    split_op_coeffs = square_augmented_so_coeffs(U, tau, L,  unitary_decomp)
    free_fermion_coeffs = _augmented_plaquette_free_fermionic_error_coefficients(tau/2, L)

    return add_dicts(split_op_coeffs, free_fermion_coeffs, scale_2nd=2)

def _fourth_order_suzuki_trotter_error_coefficient(U: float, tau: float, L: int) -> float:
    # From:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.108.195105
    return L*L * (2.1485 * tau**5 + 92.1642 * tau**4 * U + 14.3445 * tau**3 * U**2
                   + 1.0712 * tau**2 * U**3 + 0.07938 * tau*U**4)

def second_order_error_coefficient_alt(U: float, tau: float, L: int) -> float:
    # From: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.108.195105
    return L*L * (4.4142 * tau**3 + 8.0889 * tau**2 * U + 1.3062 * tau* U**2)

def second_order_error_coefficient(U: float, tau: float, L: int) -> float:
    G_red, G_blue, G = plaquette_decomposition_graphs(L)
    WSO = plaquette_second_order_split_operator_error_coefficient(U, tau, L**2, G)

    commutator_graph = ff_commutator(
        G_red, ff_commutator(G_blue, G_red)
    )

    return WSO + tau**2 * 3 / 24 * spectral_norm_of_free_fermionic_operator(commutator_graph)

def plaquette_decomposition_graphs(L: int) -> tuple[nx.Graph, nx.Graph, nx.Graph]:
    assert not (L % 2)
    red, blue = nx.Graph(), nx.Graph()

    whole = nx.grid_2d_graph(L, L, periodic=True)

    nodes = [(i, j) for i in range(L) for j in range(L)]

    for g in [red, blue]:
        g.add_nodes_from(nodes)
    for i in range(L):
        for j in range(L):
            if ((i + 1) % 2) and ((j + 1) % 2):
                red.add_edges_from([
                    ((i, j), ((i + 1) % L, j)),
                    ((i, j), (i, (j + 1) % L)),
                    (((i + 1) % L, j), ((i + 1) % L, (j + 1) % L)),
                    ((i, (j + 1) % L), ((i + 1) % L, (j + 1) % L))
                ])

            elif (i % 2) and (j % 2):
                blue.add_edges_from([
                    ((i, j), ((i + 1) % L, j)),
                    ((i, j), (i, (j + 1) % L)),
                    (((i + 1) % L, j), ((i + 1) % L, (j + 1) % L)),
                    ((i, (j + 1) % L), ((i + 1) % L, (j + 1) % L))
                ])

    return red, blue, whole

def plaquette_correctors(L: int) -> tuple[tuple[float], tuple[nx.Graph, ...]]:
    assert not (L % 2)
    m,n = L,L
    red, blue = nx.Graph(), nx.Graph()

    whole = nx.grid_2d_graph(m, n, periodic=True)

    nodes = [(i, j) for i in range(m) for j in range(n)]

    for g in [red, blue]:
        g.add_nodes_from(nodes)

    for i in range(m):
        for j in range(n):
            if ((i + 1) % 2) and ((j + 1) % 2):
                red.add_edges_from([
                    (((i - 1) % m, (j - 1) % n), ((i + 2) % m, (j - 1) % n)),
                    (((i - 1) % m, (j - 1) % n), ((i - 1) % m, (j + 2) % n)),
                    (((i + 2) % m, (j - 1) % n), ((i + 2) % m, (j + 2) % n)),
                    (((i - 1) % m, (j + 2) % n), ((i + 2) % m, (j + 2) % n)),
                ], weight=-2 / 12)

            elif (i % 2) and (j % 2):
                blue.add_edges_from([
                    (((i - 1) % m, (j - 1) % n), ((i + 2) % m, (j - 1) % n)),
                    (((i - 1) % m, (j - 1) % n), ((i - 1) % m, (j + 2) % n)),
                    (((i + 2) % m, (j - 1) % n), ((i + 2) % m, (j + 2) % n)),
                    (((i - 1) % m, (j + 2) % n), ((i + 2) % m, (j + 2) % n)),
                ], weight=2 / 24)

    return (1 / 6, -1 / 12), (red, blue)

def _augmented_plaquette_free_fermionic_error_coefficients(tau: float, L: int) -> tuple[float,...]:
    A, B, H = plaquette_decomposition_graphs(L)
    A = cast_data_to_array(A)
    B = cast_data_to_array(B)
    H = cast_data_to_array(H)
    Gs = [A, B]
    coeffs, correctors = plaquette_correctors(L)
    a, b = coeffs
    correctors = [cast_data_to_array(x) for x in correctors]

    F = sum(correctors)
    G = sum(correctors) + a * A + b * B

    W5, W6, W7, W9 = 0.0, 0.0, 0.0, 0.0

    augmented_coeffs = augmented_commutator_bound(2)
    for comm, c in augmented_coeffs.items():
        if "F" in comm:
            op = ff_commutator(Gs[comm[0] - 1], ff_commutator(Gs[comm[1] - 1], G))
        else:
            Nh = sum([1 for x in comm if x == 'H'])
            match Nh:
                case 0:
                    op = ff_commutator(Gs[comm[0] - 1],
                                       ff_commutator(Gs[comm[1] - 1],
                                                     ff_commutator(Gs[comm[2] - 1],
                                                                   ff_commutator(
                                                                       Gs[comm[3] - 1], Gs[comm[4] - 1],
                                                                   )
                                                                   )
                                                     )
                                       )

                case 1:
                    op = ff_commutator(H,
                                       ff_commutator(Gs[comm[1] - 1],
                                                     ff_commutator(Gs[comm[2] - 1],
                                                                   ff_commutator(
                                                                       Gs[comm[3] - 1], Gs[comm[4] - 1],
                                                                   )
                                                                   )
                                                     )
                                       )

                case 2:
                    op = ff_commutator(H,
                                       ff_commutator(H,
                                                     ff_commutator(Gs[comm[2] - 1],
                                                                   ff_commutator(
                                                                       Gs[comm[3] - 1], Gs[comm[4] - 1],
                                                                   )
                                                                   )
                                                     )
                                       )

                case 3:
                    op = ff_commutator(H,
                                       ff_commutator(H,
                                                     ff_commutator(H,
                                                                   ff_commutator(
                                                                       Gs[comm[3] - 1], Gs[comm[4] - 1],
                                                                   )
                                                                   )
                                                     )
                                       )

                case 4:
                    op = ff_commutator(H,
                                       ff_commutator(H,
                                                     ff_commutator(H,
                                                                   ff_commutator(
                                                                       H, Gs[comm[4] - 1],
                                                                   )
                                                                   )
                                                     )
                                       )
        W5 += tau ** 5 * spectral_norm_of_free_fermionic_operator(op) * c

    W5 += tau ** 5 * a * spectral_norm_of_free_fermionic_operator(ff_commutator(B, ff_commutator(B, A))) / 8

    for i in range(len(correctors)):
        for j in range(i + 1, len(correctors)):
            W6 += tau ** 6 * spectral_norm_of_free_fermionic_operator(ff_commutator(correctors[i], correctors[j])) / 2

    W7 += tau ** 7 * a ** 2 * spectral_norm_of_free_fermionic_operator(ff_commutator(A, ff_commutator(B, A))) / 8
    W7 += tau ** 7 * a * spectral_norm_of_free_fermionic_operator(ff_commutator(F + b*B, ff_commutator(B, A))) / 8

    W9 += tau ** 9 * a * spectral_norm_of_free_fermionic_operator(ff_commutator(F + b * B, ff_commutator(F + b * B, A))) / 12
    W9 += tau ** 9 * b * spectral_norm_of_free_fermionic_operator(ff_commutator(F, ff_commutator(F, A))) / 12
    W9 += tau ** 9 * a**2 * spectral_norm_of_free_fermionic_operator(ff_commutator(A, ff_commutator(F + b * B, B))) / 24
    W9 += tau ** 9 * b**2 * spectral_norm_of_free_fermionic_operator(ff_commutator(B, ff_commutator(F, B))) / 24
    return {5: 2*W5, 6: 2*W6, 7: 2*W7, 9: 2*W9}

def plaquette_decomposition_permutations(L: int) -> list[int]:
    assert not (L % 2)
    N = L**2
    red_order = {}
    blue_order = {}

    c1, c2 = 0,0
    for i in range(N // 4):
        red_order[4 * i + 0] = (
            c1,c2
        )
        red_order[4 * i + 1] = (
            c1 + 1,c2
        )
        red_order[4 * i + 2] = (
            c1,c2+1
        )
        red_order[4 * i + 3] = (
            c1+1,c2+1
        )

        blue_order[4 * i + 0] = (
            c1, c2
        )
        blue_order[4 * i + 1] = (
            (c1 - 1) % L,(c2) % L
        )
        blue_order[4 * i + 2] = (
            (c1 ) % L,(c2-1) % L
        )
        blue_order[4 * i + 3] = (
            (c1 - 1) % L,(c2-1) % L
        )
        if 2*(i+1) % L == 0:
            c1 = 0
            c2 += 2
        else:
            c1 += 2

    blue_order_rev = {y:x for x,y in blue_order.items()}
    permutation = []
    for i in range(N):
        permutation.append(
            blue_order_rev[red_order[i]]
        )

    return permutation

def augmented_permutations(L: int) -> dict[str: list]:
    res = {}
    assert not (L % 2)
    N = L**2
    red_order = {}
    blue_order = {}
    reflected_order = {}

    c1, c2 = 0,0
    for i in range(N // 4):
        red_order[4 * i + 0] = (
            c1,c2
        )
        red_order[4 * i + 1] = (
            c1 + 1,c2
        )
        red_order[4 * i + 2] = (
            c1,c2+1
        )
        red_order[4 * i + 3] = (
            c1+1,c2+1
        )

        reflected_order[4 * i + 0] = (
            c1, c2
        )
        reflected_order[4 * i + 1] = (
            c1, c2 + 1
        )
        reflected_order[4 * i + 2] = (
            c1 + 1, c2
        )
        reflected_order[4 * i + 3] = (
            c1 + 1, c2 + 1
        )

        blue_order[4 * i + 0] = (
            c1, c2
        )
        blue_order[4 * i + 1] = (
            (c1 + 3) % L,(c2) % L
        )
        blue_order[4 * i + 2] = (
            (c1 ) % L,(c2+3) % L
        )
        blue_order[4 * i + 3] = (
            (c1 +3) % L,(c2+3) % L
        )
        if 2*(i+1) % L == 0:
            c1 = 0
            c2 += 2
        else:
            c1 += 2

    blue_order_rev = {y:x for x,y in blue_order.items()}
    red_order_rev = {y:x for x,y in red_order.items()}
    reflected_order_rev = {y:x for x,y in reflected_order.items()}
    permutation_rta = []
    permutation_ref = []
    for i in range(N):
        permutation_rta.append(
            blue_order_rev[red_order[i]]
        )
        permutation_ref.append(
            reflected_order_rev[red_order[i]]
        )
    res["red_to_augred"] = permutation_rta
    res["augred_to_augblue"] = plaquette_decomposition_permutations(L)
    res["red_reflection"] = permutation_ref

    return res
