import math
import warnings

import numpy as np
import scipy as sp
import networkx as nx

from scipy.optimize import fsolve, minimize
from src.hubbard_models.free_fermionic_errors import error_between_exp_of_free_fermionic
from src.hubbard_models._free_fermionic_computations import cast_data_to_array, spectral_norm_of_free_fermionic_operator, ff_commutator
from src.hubbard_models.split_operator_error_coefficients import (
    fourth_order_augmented_split_operator_error_coefficients,
    fourth_order_suzuki_trotter_split_operator_error_coefficient,
    hexagonal_augmented_so_coeffs
)
from src.bch_formula.commutator_bound import fourth_order_suzuki_trotter_commutator_bound, augmented_commutator_bound
from src.tools.single_variable_optimizer import one_d_optimizer
from src.tools.add_dict import add_dicts

class NegativeTrotterStepError(ValueError):
    """Raised when an optimizer converges, but the solution violates physical/domain constraints."""
    pass

def compute_number_of_trotter_steps(t: float, eps: float, U: float, tau: float, Lx: int, Ly: int, type: str) -> dict:
    if type == "tile":
        W = second_order_error_coefficient(U, tau, Lx, Ly)
        return {"num_trotter_steps": math.ceil(
            np.sqrt(W / eps) * t ** (3 / 2)
        )}
    elif type == "tile 4th":
        W = _fourth_order_suzuki_trotter_error_coefficient(U, tau, Lx, Ly)
        return {"num_trotter_steps": math.ceil(
            (W / eps) ** (1 / 4) * t ** (5 / 4)
        )}
    elif type == "augmented tile":
        return {"error_coefficients": _augmented_hex_error_coefficients(U, tau, Lx, Ly, unitary_decomp=True)}
    else:
        raise ValueError(f"Invalid argument for type: {type}. Must be one of: \'plaquette\', \'plaquette suzuki-trotter\' or \'plaquette augmented\'.")

def compute_evolution_time_and_number_of_simulation_circuits_for_qpe(eps: float, U: float, tau: float, Lx: int, Ly: int, type: str, unitary_decomp: bool = True) -> dict:
    if type == "tile":
        W = second_order_error_coefficient(U, tau, Lx, Ly)
        return {"num_simulation_steps": math.ceil(
            6.203 * W ** (1 / 2) / (eps ** (3 / 2))
        )}
    elif type == "tile 4th":
        W = _fourth_order_suzuki_trotter_error_coefficient(U, tau, Lx, Ly)
        return {"num_simulation_steps": math.ceil(
            4.463 * W ** (1 / 4) / (eps ** (5 / 4))
        )}
    elif type == "augmented tile":
        return {"unitary_error_coefficients": _augmented_hex_error_coefficients(U, tau, Lx, Ly, unitary_decomp=unitary_decomp)}
    else:
        raise ValueError(f"Invalid argument for type: {type}. Must be one of: \'plaquette\', \'plaquette suzuki-trotter\' or \'plaquette augmented\'.")

def second_order_error_coefficient(U: float, tau: float, Lx: int, Ly: int) -> float:
    _, _, _, G = tile_trotterization_hexagonal(Lx, Ly)
    R1 = 2*spectral_norm_of_free_fermionic_operator(G)

    return U*U*tau * R1 / 24 + 9.9*U*tau**2 * (2*Lx*Ly) / 12 + 0.8532 * tau**3 * (2*Lx*Ly)

def _augmented_hexagonal_trotter_steps(t: float, eps: float, U: float, tau: float, Lx: int, Ly: int) -> int:
    W5, W6, W7, W9 = _augmented_hex_error_coefficients(U, tau, Lx,Ly , unitary_decomp=unitary_decomp)

    f = lambda n: W5 * t**5 / n**4 + W6 * t**6 / n**6 + W7 * t**7 / n**6 + W9* t**9 / n**8 - eps

    x0 = W5**(1/4) * t**(5/4) /(eps**(1/4))
    res = fsolve(f, x0, full_output = True)
    root = res[0][0]

    return math.ceil(root)

def _augmented_hexagonal_num_simulation_circuits(eps: float, U: float, tau: float, Lx: int, Ly: int, unitary_decomp: bool = True) -> tuple[float, int]:
    W5, W6, W7, W9 = _augmented_hex_error_coefficients(U, tau, Lx,Ly , unitary_decomp=unitary_decomp)

    f = lambda t: W5 * t**5  + W6 * t**6  + W7 * t**7 + W9* t**9

    opt_func = lambda t: 0.76*np.pi/(t*eps - f(t))

    tm = (eps / W5)**(1/4)

    t, Npe = one_d_optimizer(opt_func, 0.00001, tm, strict_positive=True)

    return  t, max(math.ceil(Npe),1)

def _augmented_hex_error_coefficients(U: float, tau: float, Lx: int, Ly: int, unitary_decomp: bool) -> tuple[float,...]:
    split_op_coeffs = hexagonal_augmented_so_coeffs(U, tau, Lx, Ly,  unitary_decomp)
    free_fermionic_coeffs = _compute_augmented_trotter_error_coeficients_free_fermionic(tau/2, Lx, Ly)

    return add_dicts(split_op_coeffs, free_fermionic_coeffs, scale_2nd=2)

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
    return 2*Lx*Ly *tau*min((12.49922069753601*tau**4
                         +61.53483968248771 * tau**3 * U
                         +5.299926184200563 * tau**2 * U**2
                         +0.3702416057697711 * tau * U**3
                         +0.011193645493917945 * U**4
                         ),(25.41646912677477*tau**4
                         +19.79553438848779 * tau**3 * U
                         +11.027334180525738 * tau**2 * U**2
                         +3.015549518249892 * tau * U**3
                         +0.691517376259029 * U**4
                         ))

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

def _compute_augmented_trotter_error_coeficients_free_fermionic(tau: float, Lx: int, Ly: int) -> tuple[float]:

    A, B, C, H = colored_hexagonal_trotterization(Lx, Ly)
    A = cast_data_to_array(A)
    B = cast_data_to_array(B)
    C = cast_data_to_array(C)
    H = cast_data_to_array(H)
    Gs = [A, B, C]
    coeffs, correctors = color_correctors(Lx, Ly)
    a,b,c = coeffs
    correctors = [cast_data_to_array(x) for x in correctors]

    F = sum(correctors)
    G = sum(correctors) + a*A + b*B + c*C

    W5, W6, W7, W9 = 0.0, 0.0, 0.0, 0.0

    augmented_coeffs = augmented_commutator_bound(3)
    for comm, c in augmented_coeffs.items():
        if "F" in comm:
            op = ff_commutator(Gs[comm[0] - 1], ff_commutator(Gs[comm[1] - 1], G))
        else:
            Nh = sum([1 for x in comm if x == 'H'])
            match Nh:
                case 0:
                    op = ff_commutator( Gs[comm[0] - 1],
                        ff_commutator( Gs[comm[1] - 1],
                            ff_commutator( Gs[comm[2] - 1],
                                ff_commutator(
                                    Gs[comm[3] - 1], Gs[comm[4] - 1],
                                )
                            )
                        )
                    )

                case 1:
                    op = ff_commutator(H,
                        ff_commutator( Gs[comm[1] - 1],
                            ff_commutator( Gs[comm[2] - 1],
                                ff_commutator(
                                    Gs[comm[3] - 1], Gs[comm[4] - 1],
                                )
                            )
                        )
                    )

                case 2:
                    op = ff_commutator( H,
                        ff_commutator( H,
                            ff_commutator( Gs[comm[2] - 1],
                                ff_commutator(
                                    Gs[comm[3] - 1], Gs[comm[4] - 1],
                                )
                            )
                        )
                    )

                case 3:
                    op = ff_commutator( H,
                        ff_commutator( H,
                            ff_commutator( H,
                                ff_commutator(
                                    Gs[comm[3] - 1], Gs[comm[4] - 1],
                                )
                            )
                        )
                    )

                case 4:
                    op = ff_commutator( H,
                        ff_commutator( H,
                            ff_commutator( H,
                                ff_commutator(
                                    H, Gs[comm[4] - 1],
                                )
                            )
                        )
                    )
        W5 += tau**5 * spectral_norm_of_free_fermionic_operator(op) * c

    W5 += tau**5 * a * spectral_norm_of_free_fermionic_operator(ff_commutator(C, ff_commutator(C, A))) / 8
    W5 += tau**5 * b * spectral_norm_of_free_fermionic_operator(ff_commutator(C, ff_commutator(C, B))) / 8
    W5 += tau**5 * a * spectral_norm_of_free_fermionic_operator(ff_commutator(B, ff_commutator(B, A))) / 8

    for i in range(len(correctors)):
        for j in range(i+1,len(correctors)):
            W6 += tau**6 * spectral_norm_of_free_fermionic_operator(ff_commutator(correctors[i], correctors[j])) / 2

    W7 += tau**7 * b**2 * spectral_norm_of_free_fermionic_operator(ff_commutator(B, ff_commutator(C, B))) / 8
    W7 += tau**7 * a**2 * spectral_norm_of_free_fermionic_operator(ff_commutator(A, ff_commutator(C, A))) / 8
    W7 += tau**7 * a**2 * spectral_norm_of_free_fermionic_operator(ff_commutator(A, ff_commutator(B, A))) / 8
    W7 += tau**7 * a * spectral_norm_of_free_fermionic_operator(ff_commutator(F + b*B + c*C, ff_commutator(C, A))) / 4
    W7 += tau**7 * a * spectral_norm_of_free_fermionic_operator(ff_commutator(F + b*B + c*C, ff_commutator(B, A))) / 4
    W7 += tau**7 * b * spectral_norm_of_free_fermionic_operator(ff_commutator(F +  c*C, ff_commutator(C, B))) / 4

    W9 += tau**9 * a * spectral_norm_of_free_fermionic_operator(ff_commutator(F +a*A + b* B, ff_commutator(F + b*B + c*C, A))) / 12
    W9 += tau**9 * a**2 * spectral_norm_of_free_fermionic_operator(ff_commutator(A, ff_commutator(F + b*B + c*C, A))) / 24
    W9 += tau**9 * b * spectral_norm_of_free_fermionic_operator(ff_commutator(F +c*C, ff_commutator(F + c*C, B))) / 12
    W9 += tau**9 * c * spectral_norm_of_free_fermionic_operator(ff_commutator(F, ff_commutator(F, C))) / 12
    W9 += tau**9 * b**2 * spectral_norm_of_free_fermionic_operator(ff_commutator(B, ff_commutator(F + c*C, B))) / 24
    W9 += tau**9 * c**2 * spectral_norm_of_free_fermionic_operator(ff_commutator(C, ff_commutator(F, C))) / 24

    return {5: 2*W5, 6: 2*W6, 7: 2*W7, 9: 2*W9}