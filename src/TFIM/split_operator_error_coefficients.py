import numpy as np

def _compute_moment(interaction: list[float], p: int) -> float:
    """
    :param interaction: list of interaction coefficients for a single site
    :param p: moment
    :return:
        The value of the p'th moment of the interaction, to the power of p
    """
    return sum([abs(x)**p for x in interaction])

def second_order_coefficient(sites: np.ndarray, interaction: callable, U: float) -> float:
    N, d = sites.shape
    interaction_vals = [interaction(sites[0],y) for y in sites[1:]]

    I = _compute_moment(interaction_vals, 1)
    return (2*N*I*U/3) * min(2*U + I, U + 2*I)


def fourth_order_coefficient(sites: np.ndarray, interaction: callable, U: float) -> float:
    N, d = sites.shape
    interaction_vals = [interaction(sites[0],y) for y in sites[1:]]

    I = _compute_moment(interaction_vals, 1)
    return N * min(
        7.2704 * I * U**4 + 8.8064 * I**2 * U**3 + (7.6032 + 2.6368)*I**3 * U**2 + 1.2032 * I**4 * U,
        1.2032 * I * U**4 + 5.12 * I**2 * U**3 + (10.71 + 6.912)*I**3 * U**2 + 7.2704 * I**4 * U
    )

def augmented_coefficients(sites: np.ndarray, interaction: callable, U: float, unitary_decomp: bool = True) -> float:
    N, d = sites.shape
    interaction_vals = [interaction(sites[0],y) for y in sites[1:]]

    I = _compute_moment(interaction_vals, 1)

    W5 = N * (1.65 * I * U**4 + 3.65 * I**2 * U**3 + (3.6 + 0.38)*I**3 * U**2 + 0.58 * I**4 * U)
    W7 = N * (14.3 * I**2 * U**5 + (1.78 + 1.78) * I**3 * U**4 + (0.508 + 1.53) * I**4 * U**3)
    W9 = N * ((57.5 + 58.7) * I**3 * U**6 + (0.24 + 0.47 + 0.24 + 0.35 + 1.16) * I**4 * U**5)

    if unitary_decomp:
        I2 = _compute_moment(interaction_vals, 2)
        I3 = _compute_moment(interaction_vals, 3)
        W5 += N*(2*I**2 * U**3 + 2*I**3 * U**2 + 2 * I2 * U**3 + I3* U**2 + 3 * I * I2 * U**2)*2/9
        W7 += U**4 * N * ((0.593 * 0.593) * I**3 + 0.149 * I3 + 0.445 * I * I2)
    return {5: W5, 7: W7, 9: W9}