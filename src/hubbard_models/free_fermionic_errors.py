import numpy as np
import scipy as sp

import networkx as nx

from src.hubbard_models._free_fermionic_computations import (
    ff_commutator,
    spectral_norm_of_free_fermionic_operator,
    cast_data_to_array,
    T as FF
)

def second_order_error_coefficient_of_free_fermionic_operator(hopping_graphs: list[FF]) -> float:
    hopping_mats = [cast_data_to_array(x) for x in hopping_graphs]
    res = 0
    L = len(hopping_graphs)
    for i in range(L):
        mat = np.zeros(hopping_mats[0].shape, dtype = np.complex128)
        for j in range(i+1, L):
            mat += ff_commutator(hopping_mats[i], ff_commutator(hopping_mats[i], hopping_mats[j]))
        res += spectral_norm_of_free_fermionic_operator(mat) / 24

    for i in range(L):
        mat = np.zeros(hopping_mats[0].shape, dtype = np.complex128)
        for j in range(i + 1, L):
            for k in range(i+1, L):
                mat += ff_commutator(hopping_mats[k], ff_commutator(hopping_mats[j], hopping_mats[i]))
        res += spectral_norm_of_free_fermionic_operator(mat) / 12
    return res

def error_between_exp_of_free_fermionic(hopping_graphs: list[FF], target: FF, scale: float = 1.0, hop_scales: list[float] = None) -> float:
    if hop_scales is None:
        hop_scales = [1.0 for _ in hopping_graphs]
    hopping_mats = [scale * y*cast_data_to_array(x) for x, y in zip(hopping_graphs, hop_scales, strict = True)]
    target_mat = scale*cast_data_to_array(target)

    unitary = np.eye(hopping_mats[0].shape[0], dtype = np.complex128)
    for m in hopping_mats:
        unitary = unitary @ sp.linalg.expm(m)

    return spectral_norm_of_free_fermionic_operator(sp.linalg.logm(unitary) - target_mat)