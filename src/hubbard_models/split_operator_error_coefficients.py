import networkx as nx
import numpy as np

from src.hubbard_models._free_fermionic_computations import spectral_norm_of_free_fermionic_operator, ff_commutator

def tight_binding_square_hopping_norm(tau: float, L: int) -> float:
    ks = [2*np.pi *m / L for m in range(L)]
    energies = [tau*abs(np.cos(kx) + np.cos(ky) + np.cos(kz)) for kx in ks for ky in ks for kz in ks]
    return 2*sum(energies)

def plaquette_second_order_split_operator_error_coefficient(U: float, tau: float, N: int, hopping_graph: nx.Graph) -> float:
    Hh = tau * 2*spectral_norm_of_free_fermionic_operator(hopping_graph)
    return min(
        (U**2 / 12) * Hh + (U * tau**2 / 12) * N * (np.sqrt(5) + 8),
        (U * tau**2 / 6) * N * (np.sqrt(5) + 8) + (U**2 / 24) * Hh
    )

def cubic_second_order_split_operator_error_coefficient(U: float, tau: float, L: int) -> float:
    Hh = tight_binding_square_hopping_norm(tau, L)
    BIHI = U ** 2 * Hh

    N = L**3

    _HHI_term = N*(6**2 * 4 + 2*24)

    BHHI = (U / 2) * tau ** 2 * _HHI_term

    return min(BHHI / 12 + BIHI / 24, BHHI / 24 + BIHI / 12)

def second_order_split_operator_error_coefficient(U: float, tau: float, N: int, hopping_graph: nx.Graph) -> float:
    Hh = tau * 2*spectral_norm_of_free_fermionic_operator(hopping_graph)
    BIHI = U**2 * Hh

    _HHI_term = 0.0
    for node in hopping_graph.nodes:
        TI = nx.Graph()
        TI.add_nodes_from(hopping_graph.nodes(data = True))
        for edge in hopping_graph.edges:
            if node in edge:
                if nx.is_weighted(hopping_graph):
                    TI.add_edge(*edge, weight = hopping_graph[edge[0]][edge[1]]["weight"])
                else:
                    TI.add_edge(*edge, weight = 1.0)

        _HHI_term += 2*spectral_norm_of_free_fermionic_operator(ff_commutator(TI, hopping_graph)) + 2*spectral_norm_of_free_fermionic_operator(TI)**2

    BHHI = (U / 2) * tau**2 * _HHI_term

    return min(BHHI / 12 + BIHI / 24, BHHI / 24 + BIHI / 12)

def fourth_order_suzuki_trotter_split_operator_error_coefficient(U: float, tau: float, N: int, d: int) -> float:
    return d*tau*U*N*min(
        0.0284*U**3 + 3.5004*d*tau*U**2 + 2.6864*d**2*tau**2*U + 0.4512 * d**3*tau**3,
        0.0047*U**3 + 1.3766*d*tau*U**2 + 3.5808*d**2*tau**2*U + 2.7264 * d**3*tau**3,
    )

def fourth_order_suzuki_trotter_split_operator_error_coefficient_square(U: float, tau: float, L: int) -> float:
    return tau * U * L**2 * min(
        0.0696 * U ** 3 + 0.3935 * tau * U ** 2 + 7.623 * tau ** 2 * U + 30.6 * tau ** 3,
        0.0116 * U ** 3 + 0.3073 * tau * U ** 2 + 10.2 * tau ** 2 * U + 184.9 * tau ** 3,
    )

def fourth_order_suzuki_trotter_split_operator_error_coefficient_cubic(U: float, tau: float, L: int) -> float:
    return tau * U * L**3 * min(
        0.0804*U**3 + 0.783*tau*U**2 + 31.92*tau**2*U + 176.76 *tau**3,
        0.0133*U**3 + 0.5766*tau*U**2 + 44.86*tau**2*U + 1068.1 *tau**3,
    )

def fourth_order_augmented_split_operator_error_coefficients(U: float, tau: float, N: int, d: int, unitary_decomp: bool = True) -> tuple[float, float, float]:
    Wso = (1 / 120 * U**3 + 9091 / 2880 * d*tau*U**2 + 257 / 90 * d**2 * tau**2 * U + 3 / 5 * d**3 * tau**3) * d*tau*U*N
    WSO6 = 11 * d ** 3 * tau ** 3 * U ** 3 * N / 9
    WH2H1 = (5 / 72 * (d**3 + d**2) * tau + 1 / 12 * d * tau + (5*d + 1) / 144 * U)*d*tau**2*U**2*N
    WH2H2H1 = 1 / 288 * d**2 * tau**2 * U**4 * N
    WSO7 = (67 / 2688 * U + 55 / 21 * d * tau) * d**2*tau**2 * U**4 * N
    if unitary_decomp:
        return {5: Wso + WH2H1, 6: WH2H2H1 + WSO6, 7: WSO7}
    return {5: Wso, 6: WH2H2H1 + WSO6, 7: WSO7}

def square_augmented_so_coeffs(U: float, tau: float, L: int, unitary_decomp: bool = True) -> tuple[float]:
    d = 4
    N = L**2
    WSO5 = (
                   0.009356*U**3 +
                   0.1613*tau*U**2 +
                   3.751 * tau**2 *U +
                   14.70*tau**3
           )*tau*U*N
    WH2H1 = (1.793*tau + 0.128 * U)*tau**2*U**2*N
    WH2H1_7 = (1.703*d**3 - 0.0006 * d**2 + 0.0004 * d)*tau**3*U**4*N
    WH2H2H1 = d*(d-1)*tau**2*U**4*N / 1152
    WSO7 = (1.64*U**2 + 6.06 * d*tau * U + 4.77 * d**2*tau**2) * d**2*tau**2 * U**3 * N
    WSO9 = 190*d**4*tau**4*U**5
    if unitary_decomp:
        return {5: WSO5 + WH2H1, 6: WH2H2H1, 7: WSO7 + WH2H1_7, 9: WSO9}
    return {5: WSO5, 6: WH2H2H1, 7: WSO7, 9: WSO9}

def hexagonal_augmented_so_coeffs(U: float, tau: float, Lx: int, Ly: int, unitary_decomp: bool = True) -> tuple[float, float, float]:
    N = 2*Lx* Ly
    d = 3
    WSO5 = (
                  0.005991 * U**3
                  + 0.05626 * tau*U**2
                  + 0.2836 * tau**2 * U
                  + 0.8377 *  tau**3
          ) * tau*U*N
    WH2H1 = (0.499 * tau + 0.064 * U)*tau**2*U**2*N
    WH2H1_7 = (1.703*d**3 - 0.0006 * d**2 + 0.0004 * d)*tau**3*U**4*N
    WH2H2H1 = d*(d-1)*tau**2*U**4*N / 1152
    WSO7 = (1.64*U**2 + 6.06 * d*tau * U + 4.77 * d**2*tau**2) * d**2*tau**2 * U**3 * N
    WSO9 = 190*d**4*tau**4*U**5
    if unitary_decomp:
        return {5: WSO5 + WH2H1, 6: WH2H2H1, 7: WSO7 + WH2H1_7, 9: WSO9}
    return {5: WSO5, 6: WH2H2H1, 7: WSO7, 9: WSO9}

def cubic_augmented_so_coeffs(U: float, tau: float, L: int, unitary_decomp: bool = True) -> tuple[float]:
    d = 6
    N = L**3
    WSO5 = (
                   0.01081*U**3
                   + 0.3216*tau*U**2
                   + 15.09 * tau**2 *U
                   + 84.88*tau**3
           )*tau*U*N
    WH2H1 = (6.736*tau + 0.285*U)*tau**2*U**2*N
    WH2H1_7 = (1.703*d**3 - 0.0006 * d**2 + 0.0004 * d)*tau**3*U**4*N
    WH2H2H1 = d*(d-1)*tau**2*U**4*N / 1152
    WSO7 = (1.64*U**2 + 6.06 * d*tau * U + 4.77 * d**2*tau**2) * d**2*tau**2 * U**3 * N
    WSO9 = 190*d**4*tau**4*U**5
    if unitary_decomp:
        return {5: WSO5 + WH2H1, 6: WH2H2H1, 7: WSO7 + WH2H1_7, 9: WSO9}
    return {5: WSO5, 6: WH2H2H1, 7: WSO7, 9: WSO9}
