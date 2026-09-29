
def second_order_coefficient(
    interaction_strength: float,
    hopping_strength: float,
    connectivity: int,
    num_sites: int,
) -> float:
    N = num_sites
    d = connectivity
    v = hopping_strength / 2
    V = interaction_strength / 2
    return min(
        (4*(8*N + 1)*d*v * V**2 * N**3) / 12 + ((8*N + 1)*d**2*v**2*V*N**2)/24,
        (4*(8*N + 1)*d*v * V**2 * N**3) / 24 + ((8*N + 1)*d**2*v**2*V*N**2)/12,
    )

def augmented_coefficients(
    interaction_strength: float,
    hopping_strength: float,
    connectivity: int,
    num_sites: int,
) -> dict[int: float]:
    N = num_sites
    d = connectivity
    v = hopping_strength / 2
    V = interaction_strength / 2
    W5 = (
        d*v*V**4*N**5 * 256 /180
        + d**2 * v**2 * V**3 * N**4 * (1024/160 + 1024 * 17/2880 + 48*11/2880)
        + d**3 * v**3 * V**2 * N**3 * (768 * 37/5760 + 320 /5760 + 320/960)
        + d**4 * v**4 * V * N**2 * 192 * 13 / 5760
    )
    W7 = (
        d**2*v**2 * V**5 * N**6 * 134/2016
        + d**3 * v**3 * V**4 * N**5 * 37120 / 2304
        + d**4 * v**4 * V**3 * N**4 * 36352/4032
    )
    W9 = 12160 * d**4 * v**4 * V**5 * N**6

    W5 += (
            d**3 * v**3 * V**2 * N**3 * (768 + 320) / 1152
            + d**2* v**2 * V**3 *N**4 * (48 + 1024) / 576
            + (2*N**2 *d + N*d)*(2*v*V)*12*d*v**2*V * 2/1152
            + (2*N**2 *d + N*d)*(2*v*V)*16*d*v*V**2*N * 2/576
            + (2*N**2 *d + N*d)**2*(2*v*V)*12*d*v**2*V /1152
            + (2*N**2 *d + N*d)**2*(2*v*V)*16*d*v*V**2*N /576
    )

    W6 = (
        (2*N**2 + N*d)*(2*v*V)**2 * (8*N+1)*d*v*N*4/13824
        + (2*N**2 + N*d)**2*(2*v*V)**2 * (8*N+1)*d*v*N*2/13824
    )

    return {
        5: W5,
        6: W6,
        7: W7,
        9: W9
    }
