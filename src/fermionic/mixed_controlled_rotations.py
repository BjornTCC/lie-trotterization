from copy import deepcopy

import numpy as np

from qiskit import QuantumCircuit, QuantumRegister

def mixed_control_kappa_baseline(tau: float) -> QuantumCircuit:
    qc = QuantumCircuit(4)
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
    qc = QuantumCircuit(4)
    qc.cx(3,2)
    qc.h(3)
    qc.ccx(0,2,3)
    qc.ccx(1,2,3)
    qc.rz(tau, 3)
    qc.ccx(1,2,3)
    qc.ccx(0,2,3)
    qc.rz(-tau, 3)
    qc.h(3)
    qc.cx(3,2)
    return qc.reverse_bits()

if __name__ == "__main__":
    import scipy as sp
    from openfermion.ops import QubitOperator, FermionOperator
    from openfermion.linalg import get_sparse_operator
    from openfermion.transforms import jordan_wigner, normal_ordered
    from openfermion.utils import hermitian_conjugated

    from qiskit.quantum_info import Operator

    def print_non_zero_matrix_elements(matrix: np.ndarray, nbits: int, tol:float = 1e-10) -> None:
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


    tau = 0.1

    kappa_op = 1j*(FermionOperator("3^ 2", 1.0) + FermionOperator("2^ 3", 1.0)) * (FermionOperator("0^ 0", 1.0) - FermionOperator("1^ 1", 1.0))**2
    kappa_mat = fermion_to_matrix(kappa_op, tau, 4)

    qc = mixed_control_happa_baseline(tau)
    print(qc)
    print_non_zero_matrix_elements(Operator(qc).to_matrix(), 4)
    print()
    print_non_zero_matrix_elements(kappa_mat,4)
