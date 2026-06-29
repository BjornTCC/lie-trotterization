import networkx as nx
import numpy as np

from src.hubbard_models._free_fermionic_computations import spectral_norm_of_free_fermionic_operator

def second_order_split_operator_error_coefficient(U: float, tau: float, N: int, hopping_graph: nx.Graph) -> float:
    Hh = tau * spectral_norm_of_free_fermionic_operator(hopping_graph)
    return min(
        (U**2 / 12) * Hh + (U * tau**2 / 12) * N * (np.sqrt(5) + 8),
        (U * tau**2 / 6) * N * (np.sqrt(5) + 8) + (U**2 / 24) * Hh
    )

def fourth_order_suzuki_trotter_split_operator_error_coefficient(U: float, tau: float, N: int, d: int) -> float:
    return d*tau*U*N*min(
        0.9088*U**3 + 15.1296*d*tau*U**2 + 5.6256*d**2*tau**2*U + 1.2032 * d**3*tau**3,
        0.1504*U**3 + 9.0624*d*tau*U**2 + 5.5168*d**2*tau**2*U + 7.274 * d**3*tau**3
    )

def fourth_order_augmented_split_operator_error_coefficients(U: float, tau: float, N: int, d: int, unitary_decomp: bool = True) -> tuple[float, float, float]:
    Wso = (4 / 15 * U**3 + 25877 / 2880 * d*tau*U**2 + 36 / 5 * d**2 * tau**2 * U + 8 / 5 * d**3 * tau**3) * d*tau*U*N
    WH2H1 = (5 / 9 * (d**2 + d) * tau + 2 / 3 * d**3 * tau + (5*d + 1) / 36 * U)*d*tau**2*U**2*N
    WH2H2H1 = 1 / 72 * d**2 * tau**2 * U**4 * N
    WSO7 = (1 / 7 * U + 832 / 63 * d * tau) * d**2*tau**2 * U**4 * N
    if unitary_decomp:
        return Wso + WH2H1, WH2H2H1, WSO7
    return Wso, WH2H2H1, WSO7