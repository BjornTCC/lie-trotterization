from copy import deepcopy

import numpy as np

from src.fermionic_circuits.circuits import fswap
from qiskit import QuantumCircuit, QuantumRegister
from openfermion.ops import QubitOperator, FermionOperator

def mixed_control_kappa_baseline(tau: float) -> QuantumCircuit:
    register = QuantumRegister(4, name="sim")
    qc = QuantumCircuit(register)
    qc.cx(3,2)
    qc.ccx(0,2,3)
    qc.ry(-tau, 3)
    qc.ccx(0,2,3)
    qc.ccx(1,2,3)
    qc.ry(tau, 3)
    qc.ccx(1,2,3)
    qc.cx(3,2)
    return qc.reverse_bits()

def mixed_control_happa_baseline(tau: float) -> QuantumCircuit:
    register = QuantumRegister(4, name="sim")
    qc = QuantumCircuit(register)
    qc.cx(3,2)
    qc.h(3)
    qc.ccx(1,2,3)
    qc.ccx(0,2,3)
    qc.rz(tau, 3)
    qc.ccx(0,2,3)
    qc.ccx(1,2,3)
    qc.rz(-tau, 3)
    qc.h(3)
    qc.cx(3,2)
    return qc.reverse_bits()


def _boolean_change_of_basis() -> QuantumCircuit:
    register = QuantumRegister(4, name="sim")
    bool_qc = QuantumCircuit(register)
    bool_qc.x([1, 3])
    bool_qc.cx(2, 0)
    bool_qc.cx(1, 2)
    bool_qc.cx(3, 1)
    bool_qc.cx(2, 0)
    bool_qc.cx(2, 1)
    return bool_qc

def _kappa_spin_change_of_basis() -> QuantumCircuit:
    register = QuantumRegister(4, name="sim")
    qc = QuantumCircuit(register)
    qc.cx([1,3], [0,2])
    qc.cx([0,2], [1,3])
    qc.sdg([1,3])
    qc.h([1, 3])
    qc.tdg([1, 3])
    qc.h([1, 3])
    qc.s([1, 3])

    qc.cx([0,2], [1,3])
    qc.sdg([1,3])
    qc.h([1, 3])
    qc.t([1, 3])
    qc.h([1, 3])
    qc.s([1, 3])
    qc.cx([1,3], [0,2])

    return qc

def spin_conjoined_mixed_control_happa_baseline(tau: float) -> QuantumCircuit:
    register = QuantumRegister(4, name="sim")
    qc = QuantumCircuit(register, global_phase = 0*tau)

    cob = _boolean_change_of_basis()
    qc.compose(cob, qubits = register, inplace = True)

    qc.h([2,3])
    qc.mcx([0,1,3],2, ctrl_state = 0)
    qc.rz(2*tau,2)
    qc.mcx([0,1,3],2, ctrl_state = 0)
    qc.rz(-2*tau,2)
    qc.h([2,3])

    qc.compose(cob.inverse(), qubits = register, inplace = True)
    return qc.reverse_bits()

def spin_conjoined_mixed_control_kappa_baseline(tau: float) -> QuantumCircuit:
    register = QuantumRegister(4, name="sim")
    qc = QuantumCircuit(register, global_phase = 0*tau)

    cob = _boolean_change_of_basis()
    qc.compose(cob, qubits = register, inplace = True)

    qc.h([2,3])
    qc.mcx([0,1,3],2, ctrl_state = 0)
    qc.ry(2*tau,2)
    qc.mcx([0,1,3],2, ctrl_state = 0)
    qc.ry(-2*tau,2)
    qc.h([2,3])

    qc.compose(cob.inverse(), qubits = register, inplace = True)
    return qc.reverse_bits()