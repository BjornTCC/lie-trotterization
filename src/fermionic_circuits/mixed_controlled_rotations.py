from copy import deepcopy

import numpy as np

from circuits import fswap
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

def _boolean_change_of_basis() -> QuantumCircuit:
    register = QuantumRegister(4, name="sim")
    bool_qc = QuantumCircuit(register)
    bool_qc.x([1, 3])
    bool_qc.cx(3, 0)
    bool_qc.cx(3, 1)
    bool_qc.cx(3, 2)
    bool_qc.cx(1, 2)
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
    ancilla_register = QuantumRegister(1, name="c-a")
    qc = QuantumCircuit(register, ancilla_register, global_phase = 0*tau)

    cob = _boolean_change_of_basis()
    qc.compose(cob, qubits = register, inplace = True)

    aq = ancilla_register[0]

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
    ancilla_register = QuantumRegister(1, name="c-a")
    qc = QuantumCircuit(register, ancilla_register, global_phase = 0*tau)

    cob = _boolean_change_of_basis()
    qc.compose(cob, qubits = register, inplace = True)

    aq = ancilla_register[0]

    qc.h([2,3])
    qc.mcx([0,1,3],2, ctrl_state = 0)
    qc.ry(2*tau,2)
    qc.mcx([0,1,3],2, ctrl_state = 0)
    qc.ry(-2*tau,2)
    qc.h([2,3])

    qc.compose(cob.inverse(), qubits = register, inplace = True)
    return qc.reverse_bits()



if __name__ == "__main__":
    import scipy as sp
    from openfermion.ops import QubitOperator, FermionOperator
    from openfermion.linalg import get_sparse_operator
    from openfermion.transforms import jordan_wigner, normal_ordered
    from openfermion.utils import hermitian_conjugated

    from qiskit.quantum_info import Operator

    def print_non_zero_matrix_elements(matrix: np.ndarray, nbits: int, tol:float = 1e-10) -> None:


        #print(mat)
        d = 2**nbits
        for i in range(d):
            for j in range(d):
                element = matrix[i,j]
                if abs(element) > tol:
                    print(f"{np.binary_repr(i, nbits)}, {np.binary_repr(j, nbits)}: {element}")

    def fermion_to_matrix(
            op: FermionOperator,
            theta: float,
            n_qubits: int
    ) -> np.ndarray:
        qubit_op = jordan_wigner(op)
        # print(qubit_op)
        sparse_matrix = get_sparse_operator(qubit_op, n_qubits=n_qubits)
        return sp.linalg.expm(theta * sparse_matrix.todense())


    tau = np.pi/8

    #kappa_op = 1j*(FermionOperator("3^ 2", 1.0) + FermionOperator("2^ 3", 1.0)) * (FermionOperator("0^ 0", 1.0) - FermionOperator("1^ 1", 1.0))**2
    #kappa_op = FermionOperator("1^ 0") + FermionOperator("3^ 2") - FermionOperator("0^ 1") - FermionOperator("2^ 3")
    #kappa_op = (FermionOperator("3^ 2", 1.0) - FermionOperator("2^ 3", 1.0)) * (FermionOperator("0^ 0", 1.0) - FermionOperator("1^ 1", 1.0)) +  (FermionOperator("1^ 0", 1.0) - FermionOperator("0^ 1", 1.0)) * (FermionOperator("2^ 2", 1.0) - FermionOperator("3^ 3", 1.0))
    kappa_op = 1j*(FermionOperator("3^ 2", 1.0) + FermionOperator("2^ 3", 1.0)) * (
                FermionOperator("0^ 0", 1.0) - FermionOperator("1^ 1", 1.0))**2 + 1j*(
                           FermionOperator("1^ 0", 1.0) + FermionOperator("0^ 1", 1.0)) * (
                           FermionOperator("2^ 2", 1.0) - FermionOperator("3^ 3", 1.0))**2

    kappa_mat = fermion_to_matrix(kappa_op, tau, 4)

    qc = spin_conjoined_mixed_control_happa_baseline(tau)
    qc_mat = Operator(qc).to_matrix()
    print(qc)
    print_non_zero_matrix_elements(qc_mat[::2,::2], 4)
    print()
    print_non_zero_matrix_elements(kappa_mat,4)
    print(np.allclose(qc_mat[::2,::2], kappa_mat))


