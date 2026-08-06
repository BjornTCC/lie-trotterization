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
    fourth_order_augmented_split_operator_error_coefficients
)
from src.hubbard_models.free_fermionic_errors import error_between_exp_of_free_fermionic

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
        #W = fourth_order_suzuki_trotter_split_operator_error_coefficient(U, tau, L**2,4)
        W = _fourth_order_suzuki_trotter_error_coefficient(U, tau, L)
        return math.ceil(
            (W / eps)**(1/4) * t**(5/4)
        )
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
        #W = fourth_order_suzuki_trotter_split_operator_error_coefficient(U, tau, L**2, 4)
        W = _fourth_order_suzuki_trotter_error_coefficient(U, tau, L)
        return (eps / (5*W))**(1/4), math.ceil(
            4.463 * W ** (1 / 4) /(eps ** (5 / 4))
        )
    elif type == "augmented plaquette":
        return _augmented_plaquette_num_simulation_circuits(eps, U, tau, L, unitary_decomp=unitary_decomp)
    else:
        raise ValueError(f"Invalid argument for type: {type}. Must be one of: \'plaquette\', \'plaquette suzuki-trotter\' or \'plaquette augmented\'.")

def _augmented_plaquette_trotter_steps(t: float, eps: float, U: float, tau: float, L: int) -> int:
    #W5, W6, W7 = fourth_order_augmented_split_operator_error_coefficients(U, tau, L**2, 4, unitary_decomp=unitary_decomp)
    W5, W6, W7 = _augmented_plaquette_error_coefficients(U, tau, L,  unitary_decomp=True)

    #f = lambda n: W5 * t**5 / n**4 + W6 * t**6 / n**6 + W7 * t**7 / n**6 + n * _augmented_free_fermionic_formula_error(t/n, U, tau, L) - eps
    f = lambda n: W5 * t**5 / n**4 + W6 * t**6 / n**6 + W7 * t**7 / n**6 - eps

    x0 = W5**(1/4) * t**(5/4) /(eps**(1/4))
    res = fsolve(f, x0, full_output = True)
    root = res[0][0]

    return math.ceil(root)

def _augmented_plaquette_num_simulation_circuits(eps: float, U: float, tau: float, L: int, unitary_decomp: bool = True) -> tuple[float, int]:
    #W5, W6, W7 = fourth_order_augmented_split_operator_error_coefficients(U, tau, L**2, 4, unitary_decomp=unitary_decomp)
    W5, W6, W7 = _augmented_plaquette_error_coefficients(U, tau, L,  unitary_decomp=unitary_decomp)
    #f = lambda t: W5 * t**5 + W6 * t**6 + W7 * t**7 + _augmented_free_fermionic_formula_error(t, U, tau, L)
    f = lambda t: W5 * t**5 + W6 * t**6 + W7 * t**7
    opt_func = lambda t: 0.76*np.pi/(t*eps - f(t))

    t0 = (eps / (5*W5))**(1/4)
    tm = (eps / W5)**(1/4)

    bnds = [(0.00001,tm)]# These bounds ensure Npe > 0

    min_res = minimize_scalar(opt_func, bounds = bnds[0])

    if not min_res.success:
        warnings.warn(f"Npe optimizer failed with message: {min_res.message}")

    Npe = min_res.fun

    return  min_res.x, math.ceil(Npe)

def _augmented_plaquette_error_coefficients(U: float, tau: float, L: int, unitary_decomp: bool) -> tuple[float,...]:
    WSO5, WSO6, WSO7 = fourth_order_augmented_split_operator_error_coefficients(U, tau, L**2, 4, unitary_decomp)
    WFF5, WFF6, WFF7 = _augmented_plaquette_free_fermionic_error_coefficients(U, tau, L)

    return WFF5 + WSO5, WFF6 + WSO6, WFF7 + WSO7

def _augmented_plaquette_free_fermionic_error_coefficients(U: float, tau: float, L: int) -> tuple[float,...]:
    A, B, _ = plaquette_decomposition_graphs(L)

    a = -tau**2 / 6
    b = tau**2 / 12

    A = tau*cast_data_to_array(A)
    B = tau*cast_data_to_array(B)

    BA = ff_commutator(B,A)

    ABA = ff_commutator(
        A, BA
    )
    BBA = ff_commutator(
        B, BA
    )

    BBBA = ff_commutator(B, BBA)
    BABA = ff_commutator(B, ABA)

    AAABA = ff_commutator(A, ff_commutator(A, ABA))
    BBABA = ff_commutator(B, ff_commutator(B, ABA))
    BAABA = ff_commutator(B, ff_commutator(A, ABA))
    BBBBA = ff_commutator(B, ff_commutator(B, BBA))

    F1 = ABA/12 - b * B
    F2 = BBA/24 - a * A

    F12 = ff_commutator(F1, F2)

    F = F1 + F2

    FA = ff_commutator(F,A)
    FB = ff_commutator(F,B)
    FBA = ff_commutator(F, BA)
    FFA = ff_commutator(F, ff_commutator(F, A))
    FFB = ff_commutator(F, ff_commutator(F, B))
    FBBA = ff_commutator(F, BBA)
    FBBBA = ff_commutator(F, BBBA)

    W5 = (
        (3*b+a) * spectral_norm_of_free_fermionic_operator(ABA) / 20
        +(26*b + 13 *a) * spectral_norm_of_free_fermionic_operator(BBA) / 80
        +spectral_norm_of_free_fermionic_operator(AAABA) / 1920
        +spectral_norm_of_free_fermionic_operator(BBABA) / 320
        +spectral_norm_of_free_fermionic_operator(BAABA) / 480
        +spectral_norm_of_free_fermionic_operator(BBBBA) / 384
        +spectral_norm_of_free_fermionic_operator(FBA) / 20
    )
    W6 = (
        spectral_norm_of_free_fermionic_operator(F12) / 2
        +a*b*spectral_norm_of_free_fermionic_operator(BA)
        +b*spectral_norm_of_free_fermionic_operator(BBBA) / 48
        +b*spectral_norm_of_free_fermionic_operator(BABA) / 96
        +a*spectral_norm_of_free_fermionic_operator(FA) / 4
        +b*spectral_norm_of_free_fermionic_operator(FB) / 4
        +spectral_norm_of_free_fermionic_operator(FBBA)/96
    )
    W7 = (
        (41*a*b+27*b*b)*spectral_norm_of_free_fermionic_operator(BBA) / 224
        +(3*a*b+a*a)*spectral_norm_of_free_fermionic_operator(ABA) / 56
        +b*spectral_norm_of_free_fermionic_operator(BAABA) / 672
        +b*spectral_norm_of_free_fermionic_operator(BBABA)/224
        +b*spectral_norm_of_free_fermionic_operator(BBBBA)/224
        +spectral_norm_of_free_fermionic_operator(FFA)/28
        +spectral_norm_of_free_fermionic_operator(FFB)/28
        +(3*a+b)*spectral_norm_of_free_fermionic_operator(FBA) / 28
        +spectral_norm_of_free_fermionic_operator(FBBBA) / 672
    )

    return 2*W5, 2*W6, 2*W7 # Factor of 2 is due to spin sectors.

def _fourth_order_suzuki_trotter_error_coefficient(U: float, tau: float, L: int) -> float:
    H1, H2, _ = plaquette_decomposition_graphs(L)

    hoppings = [H1, H2]

    N = L**2
    d = 2
    res = 0
    for term, coeff in childs_commutator_error_coeffs.items():
        if (term[0], term[1],term[2], term[3]) == (3,3,3,3):
            res += coeff * d * tau * U ** 4 * N
        elif (term[0], term[1], term[3]) == (3,3,3):
            res += coeff * 384 * d * tau * U ** 4 * N
        elif (term[0], term[2],term[3]) == (3,3,3):
            res += coeff * 128 * d * tau * U ** 4 * N
        elif (term[1],term[2], term[3]) == (3,3,3):
            res += coeff * 6 * d * tau * U ** 4 * N
        elif (term[0],term[1], term[2]) == (3,3,3):
            res += coeff * 128 * d * tau * U ** 4 * N
        elif (term[0], term[3]) == (3,3):
            res += coeff * 192 * d * tau * U ** 4 * N
        elif (term[1], term[3]) == (3,3):
            res += coeff * 240 * d * tau * U ** 4 * N
        elif (term[2], term[3]) == (3,3):
            res += coeff * 80 * d * tau * U ** 4 * N
        elif (term[1], term[2]) == (3,3):
            res += coeff * 128 * d * tau * U ** 4 * N
        elif (term[0], term[2]) == (3,3):
            res += coeff * 192 * d * tau * U ** 4 * N
        elif (term[0], term[1]) == (3,3):
            res += coeff * 64 * d * tau * U ** 4 * N
        elif term[0] == 3:
            res += coeff * 32 * d * tau * U ** 4 * N
        elif term[1] == 3:
            res += coeff * 48 * d * tau * U ** 4 * N
        elif term[2] == 3:
            res += coeff * 64 * d * tau * U ** 4 * N
        elif term[3] == 3:
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
    Mat_Gr = 1j*(t * tau/2 - t**3 * tau**3 / 12) * Gr_arr
    Mat_Gb = 1j*(t * tau/2 + t**3 * tau**3 / 24) * Gb_arr
    Mat_Cr = 1j*t**3 * tau**3 * (cast_data_to_array(correction) - 2*Gb_arr) / 24
    Mat_Cb = 1j*t**3 * tau**3 * (cast_data_to_array(correction2) + 2*Gr_arr) / 12
    Mat_G = 1j*t * tau * cast_data_to_array(G)

    seq = [Mat_Gr, Mat_Gb, Mat_Cr, Mat_Cb, Mat_Gb, Mat_Gr]

    return 2*error_between_exp_of_free_fermionic(seq, Mat_G)

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

childs_commutator_error_coeffs = {
    (1,1,1,2,1): 0.0047,
    (1,1,2,2,1): 0.0057,
    (1,1,3,2,1): 0.0057,
    (1,2,1,2,1): 0.0046,
    (1,2,2,2,1): 0.0074,
    (1,2,3,2,1): 0.0082,
    (1,3,1,2,1): 0.0046,
    (1,3,2,2,1): 0.0070,
    (1,3,3,2,1): 0.0082,
    (2,1,1,2,1): 0.0150,
    (2,1,2,2,1): 0.0161,
    (2,1,3,2,1): 0.0161,
    (2,2,1,2,1): 0.0239,
    (2,2,2,2,1): 0.0315,
    (2,2,3,2,1): 0.0303,
    (2,3,1,2,1): 0.0179,
    (2,3,2,2,1): 0.0232,
    (2,3,3,2,1): 0.0259,
    (3,1,1,2,1): 0.0204,
    (3,1,2,2,1): 0.0225,
    (3,1,3,2,1): 0.0225,
    (3,2,1,2,1): 0.0423,
    (3,2,2,2,1): 0.0585,
    (3,2,3,2,1): 0.0502,
    (3,3,1,2,1): 0.0423,
    (3,3,2,2,1): 0.0681,
    (3,3,3,2,1): 0.0648,
    (1,1,1,3,1): 0.0047,
    (1,1,2,3,1): 0.0057,
    (1,1,3,3,1): 0.0057,
    (1,2,1,3,1): 0.0046,
    (1,2,2,3,1): 0.0070,
    (1,2,3,3,1): 0.0082,
    (1,3,1,3,1): 0.0046,
    (1,3,2,3,1): 0.0058,
    (1,3,3,3,1): 0.0074,
    (2,1,1,3,1): 0.0150,
    (2,1,2,3,1): 0.0161,
    (2,1,3,3,1): 0.0161,
    (2,2,1,3,1): 0.0239,
    (2,2,2,3,1): 0.0306,
    (2,2,3,3,1): 0.0303,
    (2,3,1,3,1): 0.0179,
    (2,3,2,3,1): 0.0206,
    (2,3,3,3,1): 0.0241,
    (3,1,1,3,1): 0.0204,
    (3,1,2,3,1): 0.0225,
    (3,1,3,3,1): 0.0225,
    (3,2,1,3,1): 0.0423,
    (3,2,2,3,1): 0.0571,
    (3,2,3,3,1): 0.0502,
    (3,3,1,3,1): 0.0423,
    (3,3,2,3,1): 0.0641,
    (3,3,3,3,1): 0.0621,
    (1,1,1,3,2): 0.0043,
    (1,1,2,3,2): 0.0057,
    (1,1,3,3,2): 0.0057,
    (1,2,1,3,2): 0.0035,
    (1,2,2,3,2): 0.0062,
    (1,2,3,3,2): 0.0082,
    (1,3,1,3,2): 0.0035,
    (1,3,2,3,2): 0.0046,
    (1,3,3,3,2): 0.0074,
    (2,1,1,3,2): 0.0141,
    (2,1,2,3,2): 0.0161,
    (2,1,3,3,2): 0.0161,
    (2,2,1,3,2): 0.0212,
    (2,2,2,3,2): 0.0290,
    (2,2,3,3,2): 0.0303,
    (2,3,1,3,2): 0.0153,
    (2,3,2,3,2): 0.0179,
    (2,3,3,3,2): 0.0241,
    (3,1,1,3,2): 0.0186,
    (3,1,2,3,2): 0.0217,
    (3,1,3,3,2): 0.0225,
    (3,2,1,3,2): 0.0377,
    (3,2,2,3,2): 0.0537,
    (3,2,3,3,2): 0.0502,
    (3,3,1,3,2): 0.0377,
    (3,3,2,3,2): 0.0601,
    (3,3,3,3,2): 0.0628,
}