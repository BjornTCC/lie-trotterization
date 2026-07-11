import numpy as np
import scipy as sp

from itertools import combinations

from openfermion.ops import QubitOperator, FermionOperator
from openfermion.linalg import get_sparse_operator
from openfermion.transforms import jordan_wigner, normal_ordered
from openfermion.utils import hermitian_conjugated

from qiskit.quantum_info import Operator

from circuits import single_normal_ordered_fermion_operator_circuit

def fermion_circuit_matrix(
        op: FermionOperator,
        theta: float,
        n_qubits: int,
) -> np.ndarray:
    # must fix endianess
    qc = single_normal_ordered_fermion_operator_circuit(op, theta, n_qubits)
    mat = Operator(qc).to_matrix()
    if "c-a" in [reg.name for reg in qc.qregs]:
        # pick out the section in which the ancilla qubit = 0
        return mat[::2, ::2]
    return mat

def compare_circuit_to_matrix(
        operator: FermionOperator,
        theta: float,
        n_qubits: int,
        verbose: int = 1
) -> float:

    circuit_matrix = fermion_circuit_matrix(operator, theta, n_qubits)
    if operator == -normal_ordered(hermitian_conjugated(operator)):
        op = operator
    else:
        op = operator- normal_ordered(hermitian_conjugated(operator))

        if len(op.terms) == 1:
            op = op / 2
    check_matrix = fermion_circuit_matrix(op, theta, n_qubits)
    if verbose:
        print("\nop:")
        print(op)
        print("circuit:")
        print(single_normal_ordered_fermion_operator_circuit(operator, theta, n_qubits))
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
    max_order = 6
    max_qubits = 6
    theta = np.pi/8

    for order in range(1, max_order + 1):
        if order % 2 != 0:
            continue
        for coeff in [1.0, 1j, 2 + 3*1j]:
            for dec_order in [order // 2]:
            #for dec_order in range(order + 1):
                qubits_ind = list(range(max_qubits))
                inc_order = order - dec_order
                inc_inds = list(get_sorted_subsets_of_size(qubits_ind, inc_order))
                dec_inds = list(get_sorted_subsets_of_size(qubits_ind, dec_order))
                for inc_ind in inc_inds:
                    for dec_ind in dec_inds:
                        string = "".join([f"{i}^ " for i in inc_ind[::-1]]) + " " + "".join([f"{i} " for i in dec_ind[::-1]])
                        op_part = FermionOperator(string, coeff) - normal_ordered(hermitian_conjugated(FermionOperator(string, coeff)))
                        if len(op_part.terms) == 1:
                            op_part /= 2
                        print(op_part)
                        try:
                            val = compare_circuit_to_matrix(op_part, theta, max_qubits, verbose=False)

                            if val > 1e-12:
                                print(f"Operator {op_part} has an error:")
                                print(compare_circuit_to_matrix(op_part, theta, max_qubits, verbose=True))
                                exit()

                        except Exception as e:
                            print(f"Operator {op_part} failed with: {e}")
                            exit()
    print("success")