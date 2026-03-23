import numpy as np
import scipy as sp

from itertools import combinations

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
    mat = Operator(qc).to_matrix()
    if "c-a" in [reg.name for reg in qc.qregs]:
        # pick out the section in which the ancilla qubit = 0
        return mat[::2, ::2]
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
        print("\ncircuit:")
        print(single_no_fermion_operator_circuit(operator, theta, n_qubits, real))
        nq = n_qubits
        print("\nNon-zero matrix elements of check matrix:")
        for i in range(2 ** nq):
            for j in range(2 ** nq):
                if abs(check_matrix[i, j]) > 1e-10:
                    print(format(i, f'0{nq}b') + "," + format(j, f'0{nq}b') + f": {check_matrix[i, j]}")
        print("\nNon-zero matrix elements of circuit mat:")
        for i in range(2 ** nq):
            for j in range(2 ** nq):
                if abs(circuit_matrix[i, j]) > 1e-10:
                    print(format(i, f'0{nq}b') + "," + format(j, f'0{nq}b') + f": {circuit_matrix[i, j]}")
    return np.linalg.norm(check_matrix - circuit_matrix)

def get_sorted_subsets_of_size(iterable, size):
  """
  Generates all sorted subsets (combinations) of a given size from an iterable.

  Args:
    iterable: The input iterable (list, tuple, string, etc.).
    size: The desired size (length) of the subsets.

  Yields:
    Each subset as a sorted tuple.
  """
  # The combinations function handles the generation and sorting internally
  for subset in combinations(iterable, size):
    yield subset
"""
nq = 2
theta = 1.0

op_part = FermionOperator("1^ 0^ 1 0", 1j)

print(compare_circuit_to_matrix(op_part, theta, nq, False))
"""

if __name__ == "__main__":
    max_order = 4
    max_qubits = 4
    theta = np.pi/8

    operator = FermionOperator("3^ 2^ 1 0")
    print(operator)
    print(compare_circuit_to_matrix(operator, theta, max_qubits, True))
    operator = FermionOperator("3^ 1^ 2 0")
    print(operator)
    print(compare_circuit_to_matrix(operator, theta, max_qubits, True))

    #exit()
    qubits_ind = list(range(max_qubits))

    for order in range(1, max_order + 1):
        if order % 2 != 0:
            continue
        for r,coeff in zip([True, False], [1.0, 1j]):
            #for dec_order in [order // 2]:
            for dec_order in range(order + 1):
                inc_order = order - dec_order
                inc_inds = list(get_sorted_subsets_of_size(qubits_ind, inc_order))
                dec_inds = list(get_sorted_subsets_of_size(qubits_ind, dec_order))
                for inc_ind in inc_inds:
                    for dec_ind in dec_inds:
                        string = "".join([f"{i}^ " for i in inc_ind[::-1]]) + " " + "".join([f"{i} " for i in dec_ind[::-1]])
                        op_part = FermionOperator(string, coeff)
                        try:
                            val = compare_circuit_to_matrix(op_part, theta, max_qubits, r, verbose=False)

                            if val > 1e-12:
                                print(f"Operator {op_part} has an error:")
                                print(compare_circuit_to_matrix(op_part, theta, max_qubits, r, verbose=True))
                                exit()

                        except Exception as e:
                            print(f"Operator {op_part} failed with: {e}")
                            exit()
    print("success")