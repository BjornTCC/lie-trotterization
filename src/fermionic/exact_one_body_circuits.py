from copy import deepcopy

import numpy as np

from openfermion import FermionOperator

from qiskit import QuantumCircuit

from plane_rotations import PlaneRotation

def synthesize_one_body_circuit(
        H: FermionOperator,
        t: float,
        real: bool = False,
        imag: bool = False,
) -> QuantumCircuit:
    if (real and imag) or (not real and not imag):
        raise ValueError("Not implemented yet")

    if real:
        return _synthesize_kappa_type_rotation(H, t)

    if imag:
        return _synthesize_happa_type_rotation(H, t)

