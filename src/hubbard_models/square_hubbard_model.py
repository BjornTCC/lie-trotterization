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
    second_order_split_operator_error_coefficient,
    fourth_order_suzuki_trotter_split_operator_error_coefficient,
    fourth_order_augmented_split_operator_error_coefficients
)

from scipy.optimize import fsolve, minimize_scalar

def get_fermionic_operator(U: float, tau: float, L: int) -> FermionOperator:
    graph = nx.convert_node_labels_to_integers(nx.grid_2d_graph(L,L, periodic = True))
    return hubbard_from_nx(tau, U, graph)

def compute_number_of_trotter_steps(t: float, eps: float, U: float, tau: float, L: int, type: str) -> int:
    if type == "plaquette":
        W = second_order_error_coefficient(U, tau, L)
        return math.ceil(
            np.sqrt(W / eps) * t**(3/2)
        )
    elif type == "plaquette suzuki-trotter":
        W = fourth_order_suzuki_trotter_split_operator_error_coefficient(U, tau, L**2,4)
        return math.ceil(
            (W / eps)**(1/4) * t**(5/4)
        )
        #return _plaquette_suzuki_trotter_steps(t,eps,U,tau, L)
    elif type == "augmented plaquette":
        return _augmented_plaquette_trotter_steps(t, eps, U, tau, L)
    else:
        raise ValueError(f"Invalid argument for type: {type}. Must be one of: \'plaquette\', \'plaquette suzuki-trotter\' or \'plaquette augmented\'.")

def compute_evolution_time_and_number_of_simulation_circuits_for_qpe(eps: float, U: float, tau: float, L: int, type: str, unitary_decomp: bool = True) -> tuple[float,int]:
    if type == "plaquette":
        W = second_order_error_coefficient(U, tau, L)
        return np.sqrt(eps / (3*W)), math.ceil(
            6.203 * W**(1/2) / (eps**(3/2))
        )
    elif type == "plaquette suzuki-trotter":
        W = fourth_order_suzuki_trotter_split_operator_error_coefficient(U, tau, L**2, 4)
        return (eps / (5*W))**(1/4), math.ceil(
            4.463 * W ** (1 / 4) /(eps ** (5 / 4))
        )
        #return _plaquette_suzuki_trotter_num_simulation_circuits(eps,U,tau, L)
    elif type == "augmented plaquette":
        return _augmented_plaquette_num_simulation_circuits(eps, U, tau, L, unitary_decomp=unitary_decomp)
    else:
        raise ValueError(f"Invalid argument for type: {type}. Must be one of: \'plaquette\', \'plaquette suzuki-trotter\' or \'plaquette augmented\'.")

def _augmented_plaquette_trotter_steps(t: float, eps: float, U: float, tau: float, L: int) -> int:
    W5, W6, W7 = fourth_order_augmented_split_operator_error_coefficients(U, tau, L**2, 4)

    f = lambda n: W5 * t**5 / n**4 + W6 * t**6 / n**6 + W7 * t**7 / n**6 + n * _augmented_free_fermionic_formula_error(t/n, U, tau, L) - eps

    x0 = W5**(1/4) * t**(5/4) /(eps**(1/4))
    res = fsolve(f, x0, full_output = True)
    root = res[0][0]

    return math.ceil(root)

def _augmented_plaquette_num_simulation_circuits(eps: float, U: float, tau: float, L: int, unitary_decomp: bool = True) -> tuple[float, int]:
    W5, W6, W7 = fourth_order_augmented_split_operator_error_coefficients(U, tau, L**2, 4, unitary_decomp=unitary_decomp)

    f = lambda t: W5 * t**5 + W6 * t**6 + W7 * t**7 + _augmented_free_fermionic_formula_error(t, U, tau, L)
    opt_func = lambda t: 0.76*np.pi/(t*eps - f(t))

    t0 = (eps / (5*W5))**(1/4)
    tm = (eps / W5)**(1/4)

    bnds = [(0.00001,tm)]# These bounds ensure Npe > 0

    min_res = minimize_scalar(opt_func, bounds = bnds[0])

    if not min_res.success:
        warnings.warn(f"Npe optimizer failed with message: {min_res.message}")

    Npe = min_res.fun

    return  min_res.x, math.ceil(Npe)


def _plaquette_suzuki_trotter_steps(t: float, eps: float, U: float, tau: float, L: int) -> int:
    ...

def _plaquette_suzuki_trotter_num_simulation_circuits(eps: float, U: float, tau: float, L: int) -> int:
    ...

def _augmented_free_fermionic_formula_error_coefficients(t: float, U: float, tau: float, L: int) -> float:
    Gr, Gb, G = plaquette_decomposition_graphs(L)

    ...

    return W5, W6, W7


def _augmented_free_fermionic_formula_error(t: float, U: float, tau: float, L: int) -> float:
    Gr, Gb, G = plaquette_decomposition_graphs(L)

    correction = ff_commutator(
        Gr, ff_commutator(Gb, Gr)
    )
    correction2 = ff_commutator(
        Gb, ff_commutator(Gb, Gr)
    )
    Gb_arr = cast_data_to_array(Gb)
    Gr_arr = cast_data_to_array(Gr)
    Mat_Gr = (t * tau/2 - t**3 * tau**3 / 12) * Gr_arr
    Mat_Gb = (t * tau/2 + t**3 * tau**3 / 24) * Gb_arr
    Mat_Cr = -t**3 * tau**3 * (cast_data_to_array(correction) + 2*cast_data_to_array(correction2) + 2*Gb_arr- 4*Gr_arr) / 24
    Mat_G = t * tau * cast_data_to_array(G)

    comp_mat = sp.linalg.expm(Mat_Gr) @ sp.linalg.expm(Mat_Gb) @ sp.linalg.expm(Mat_Cr) @ sp.linalg.expm(Mat_Gb) @ sp.linalg.expm(Mat_Gr)
    exact_mat = sp.linalg.expm(Mat_G)

    return 2*spectral_norm_of_free_fermionic_operator(
        comp_mat - exact_mat
    )

def second_order_error_coefficient(U: float, tau: float, L: int) -> float:
    G_red, G_blue, G = plaquette_decomposition_graphs(L)
    WSO = second_order_split_operator_error_coefficient(U, tau, L**2, G)

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
