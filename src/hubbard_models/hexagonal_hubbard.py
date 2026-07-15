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
from src.hubbard_models.free_fermionic_errors import error_between_exp_of_free_fermionic

from scipy.optimize import fsolve, minimize_scalar

def s1_tile_trotterization_graphs(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
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

def s2_tile_trotterization_graphs(Lx: int, Ly: int) -> tuple[nx.Graph, ...]:
    if (Ly % 2):
        raise ValueError(f"For this tiling Ly({Ly}) must be even")
    whole = nx.hexagonal_lattice_graph(m=Ly, n=2*Lx, periodic=True, with_positions=True)

    lx, ly = 2*Lx, 2*Ly

    Gr, Gb, Gy = nx.Graph(), nx.Graph(), nx.Graph()
    for g in [Gr, Gb, Gy]:
        g.add_nodes_from(whole.nodes)

    # Add red edges
    for i in range(Lx):
        for j in range(Ly // 2):
            x = (2*i) % lx
            y = (4 * j + 1) % ly
            Gr.add_edge(
                    (x,y),(x, (y + 1) % ly)
            )
            Gr.add_edge(
                    (x,y),(x, (y - 1) % ly)
            )

            x = (2*i + 1) % lx
            y = (4 * j + 2) % ly
            Gr.add_edge(
                    (x,y),(x % lx, (y + 1) % ly)
            )
            Gr.add_edge(
                    (x,y),(x % lx, (y - 1) % ly)
            )
    # Add blue edges
    for i in range(Lx):
        for j in range(Ly // 2):
            s = 2*(i % 2)
            x = (2*i) % lx
            y = (4 * j + 1 + s) % ly
            Gb.add_edge(
                    (x,y),(x, (y + 1) % ly)
            )
            Gb.add_edge(
                    (x,y),((x - 1) % lx, y)
            )

            x = (2*i + 1) % lx
            y = (4 * j + 1 + s) % ly
            Gb.add_edge(
                    (x,y),(x % lx, (y + 1) % ly)
            )
            Gb.add_edge(
                (x, y), ((x + 1) % lx, y)
            )
    # Add yellow edges
    for i in range(Lx):
        for j in range(Ly // 2):
            x = (2 * i) % lx
            y = (4 * j + 2) % ly
            Gy.add_edge(
                (x, y), (x, (y + 1) % ly)
            )
            Gy.add_edge(
                (x, y), ((x +1) % lx, y)
            )

            x = (2 * i + 1) % lx
            y = (4 * j) % ly
            Gy.add_edge(
                (x, y), (x % lx, (y + 1) % ly)
            )
            Gy.add_edge(
                (x, y), ((x - 1) % lx, y)
            )

    return Gr, Gb, Gy, whole