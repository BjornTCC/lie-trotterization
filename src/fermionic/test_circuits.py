import numpy as np
import scipy as sp

from openfermion.ops import QubitOperator, FermionOperator
from openfermion.linalg import get_sparse_operator
from openfermion.transforms import jordan_wigner, normal_ordered
from openfermion.utils import hermitian_conjugated

from qiskit.quantum_info import Operator

from circuits import single_no_fermion_operator_circuit


def fermion_to_matrix(
        op: FermionOperator,
        theta: float,
        n_qubits: int
) -> np.ndarray:
    qubit_op = jordan_wigner(op)
    sparse_matrix = get_sparse_operator(qubit_op, n_qubits=n_qubits)
    return sp.linalg.expm(theta*sparse_matrix.todense())

def fermion_circuit_matrix(
        op: FermionOperator,
        theta: float,
        n_qubits: int,
        real: bool,
) -> np.ndarray:
    # must fix endianess
    qc = single_no_fermion_operator_circuit(op, theta, n_qubits, real)
    print(qc)
    mat = Operator(qc).to_matrix()
    if "c-a" in [reg.name for reg in qc.qregs]:
        # pick out the section in which the ancilla qubit = 0
        return mat[:2**(n_qubits+1):2, :2**(n_qubits+1):2]
    return mat

def compare_circuit_to_matrix(
        operator: FermionOperator,
        theta: float,
        n_qubits: int,
        real: bool,
        verbose: int = 1
) -> float:

    circuit_matrix = fermion_circuit_matrix(operator, theta, n_qubits, real)
    if operator == -normal_ordered(hermitian_conjugated(operator)):
        op = operator
    else:
        op = operator- hermitian_conjugated(operator)

    check_matrix = fermion_to_matrix(op, theta, n_qubits)
    if verbose:
        print("check matrix:")
        print(check_matrix)
        print("\ncircuit matrix:")
        print(circuit_matrix)
    return np.linalg.norm(check_matrix - circuit_matrix)

nq = 2
theta = 0.1

op_part = FermionOperator("1^ 0", 1j)

print(compare_circuit_to_matrix(op_part, theta, nq, False))