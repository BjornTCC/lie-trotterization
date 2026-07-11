from copy import deepcopy

import numpy as np

from qiskit import QuantumCircuit, QuantumRegister

from openfermion.ops import FermionOperator

def single_normal_ordered_fermion_operator_circuit(
        operator: FermionOperator,
        angle: float,
        num_qubits: int,
        tol: float = 1e-12,
        type: str = ""
) -> QuantumCircuit:
    if len(operator.terms) == 0:
        register = QuantumRegister(num_qubits, name="sim")
        return QuantumCircuit(register)
    assert len(operator.terms) <= 2, "Operator must contain at most 2 terms"
    assert operator.is_normal_ordered(), "Operator must be normal ordered"

    ads, coeff = list(operator.terms.items())[0]

    int_list = [x[0] for x in ads]

    ctrl_up = set([x[0] for x in ads if x[1]])
    ctrl_down = set([x[0] for x in ads if not x[1]])

    if ctrl_up == ctrl_down:
        return n_type_circuit(list(ctrl_up), num_qubits, np.imag(coeff), angle).reverse_bits()

    if abs(np.real(coeff)) < tol or type == "imag":
        _type = "imag"
    elif abs(np.imag(coeff)) < tol or type == "real":
        _type = "real"
    else:
        _type = "mixed"

    n_bits = ctrl_up.intersection(ctrl_down)

    return xy_type_circuit(
        list(ctrl_up - n_bits),
        list(ctrl_down - n_bits),
        list(n_bits),
        num_qubits,
        angle * coeff,
        type = _type,
    ).reverse_bits()

def fswap(qc: QuantumCircuit, i: int, j: int, inplace: bool = True) -> QuantumCircuit | None:
    if inplace:
        qc.swap(i,j)
        qc.cz(i,j)
    else:
        qc_new = qc.copy()
        fswap(qc_new, i,j)
        return qc_new

def get_fswap_network(
        register: QuantumRegister,
        permutation: list[int],
) -> QuantumCircuit:

    # Compute a circuit implementing the permutation with fermionic algebra
    res = QuantumCircuit(register)
    _permutation = deepcopy(permutation)
    N = len(_permutation)

    while True:
        swapped = False
        for i in range(1,N):
            if _permutation[i-1] > _permutation[i]:
                _permutation[i], _permutation[i - 1] = _permutation[i - 1], _permutation[i]
                fswap(res, i, i-1, inplace=True)
                swapped = True
        if not swapped:
            break
    return res

def compute_signs_due_to_n_bits(
        dag_indices: list[int],
        ndag_indices: list[int],
        n_bits: list[int],
) -> int:
    res = 1
    for i in n_bits:
        res *= (-1)**(len([x for x in dag_indices if x < i]) + len([x for x in ndag_indices if x > i]) + len([x for x in n_bits if x < i]))
    return res

def compute_fermionic_permutation(
    register: QuantumRegister,
    dag_indices: list[int],
    ndag_indices: list[int],
    n_bits: list[int],
) -> tuple[QuantumCircuit, tuple, tuple, list]:

    all_indices = dag_indices + ndag_indices + n_bits
    max_pos, min_pos = np.argmax(all_indices), np.argmin(all_indices)
    max_index, min_index = all_indices[max_pos], all_indices[min_pos]

    n_list = list(range(min_index, min_index + len(n_bits)))
    ndag_tup = (min_index + len(n_bits), min_index + len(n_bits) + len(ndag_indices) - 1)
    dag_tup = (ndag_tup[0] + len(ndag_indices), ndag_tup[1] + len(ndag_indices))

    # Compute the desired permutation
    inverse_permutation = {}
    for i, j in enumerate(n_bits):
        inverse_permutation[min_index + i] = j
    for i, j in enumerate(sorted(ndag_indices)):
        inverse_permutation[min_index + len(n_bits) + i] = j
    for i, j in enumerate(sorted(dag_indices)):
        inverse_permutation[min_index + len(n_bits) + len(ndag_indices) + i] = j

    remaining_indices = list(range(min_index)) + list(range(min_index + len(n_bits) + len(ndag_indices) + len(dag_indices), max_index + 1))
    permutation_dict = {x: y for y, x in inverse_permutation.items()}
    permutation = []

    rem_index = 0
    for i in range(max_index + 1):
        if i in permutation_dict.keys():
            permutation.append(permutation_dict[i])
        else:
            permutation.append(remaining_indices[rem_index])
            rem_index += 1

    fswap_network = get_fswap_network(register, permutation)

    return fswap_network, dag_tup, ndag_tup, n_list

def xy_type_circuit(
        dag_indices: list[int],
        ndag_indices: list[int],
        n_bits: list[int],
        num_qubits: int,
        angle: float | complex,
        type: str = "real"
) -> QuantumCircuit:
    register = QuantumRegister(num_qubits, name="sim")
    fswap_network, dag_tup, ndag_tup, n_list = compute_fermionic_permutation(register, dag_indices, ndag_indices, n_bits)

    # compute some signs due to permutation
    n_sign = compute_signs_due_to_n_bits(dag_indices, ndag_indices, n_bits)
    x_sign = (-1)**(len(ndag_indices) * (len(ndag_indices)-1) // 2)
    baseline_circ = xy_type_circuit_baseline(dag_tup, ndag_tup, n_list, num_qubits, x_sign*n_sign*angle, type)
    res = fswap_network.compose(baseline_circ)
    return res.compose(fswap_network.inverse())

def xy_type_circuit_baseline(
        dag_indices: tuple[int,int],
        ndag_indices: tuple[int,int],
        n_type_bits: list[int],
        num_qubits: int,
        angle: float | complex,
        type: str = "imag",
) -> QuantumCircuit:
    match type:
        case "real":
            return xy_type_circuit_baseline_real(dag_indices, ndag_indices, n_type_bits, num_qubits, np.real(angle))
        case "imag":
            return xy_type_circuit_baseline_imag(dag_indices, ndag_indices, n_type_bits, num_qubits, np.imag(angle))
        case "mixed":
            return xy_type_circuit_baseline_mixed(dag_indices, ndag_indices, n_type_bits, num_qubits, angle)
        case _:
            raise ValueError("Not implemented")

def get_change_of_basis_circuit(
        dag_indices: tuple[int,int],
        ndag_indices: tuple[int,int],
        n_type_bits: list[int],
        register: QuantumRegister
) -> QuantumCircuit:
    circuit = QuantumCircuit(register)

    target_qubit = dag_indices[1]

    change_of_basis_circuit = QuantumCircuit(register)
    for i in range(*dag_indices):
        change_of_basis_circuit.cx(target_qubit, i)
    for i in range(ndag_indices[0],ndag_indices[1]+1):
        change_of_basis_circuit.x(i)
        change_of_basis_circuit.cx(target_qubit, i)

    for i in n_type_bits:
        change_of_basis_circuit.x(i)
    return change_of_basis_circuit

def xy_type_circuit_baseline_real(
        dag_indices: tuple[int,int],
        ndag_indices: tuple[int,int],
        n_type_bits: list[int],
        num_qubits: int,
        angle: float,
) -> QuantumCircuit:
    register = QuantumRegister(num_qubits, name="sim")
    change_of_basis_circuit = get_change_of_basis_circuit(dag_indices, ndag_indices, n_type_bits, register)

    cr_circuit = QuantumCircuit(register)
    ctrl_qubits = list(range(ndag_indices[0], dag_indices[1])) + n_type_bits
    cr_circuit.mcx(ctrl_qubits, dag_indices[1], ctrl_state = 0)
    cr_circuit.ry(-angle, dag_indices[1])
    cr_circuit.mcx(ctrl_qubits, dag_indices[1], ctrl_state = 0)
    cr_circuit.ry(angle, dag_indices[1])

    res = change_of_basis_circuit.compose(cr_circuit)
    return res.compose(change_of_basis_circuit.inverse())

def xy_type_circuit_baseline_imag(
        dag_indices: tuple[int, int],
        ndag_indices: tuple[int, int],
        n_type_bits: list[int],
        num_qubits: int,
        angle: float,
) -> QuantumCircuit:
    register = QuantumRegister(num_qubits, name="sim")
    change_of_basis_circuit = get_change_of_basis_circuit(dag_indices, ndag_indices, n_type_bits, register)

    cr_circuit = QuantumCircuit(register)
    ctrl_qubits = list(range(ndag_indices[0], dag_indices[1])) + n_type_bits
    cr_circuit.h(dag_indices[1])
    cr_circuit.mcx(ctrl_qubits, dag_indices[1], ctrl_state=0)
    cr_circuit.rz(angle, dag_indices[1])
    cr_circuit.mcx(ctrl_qubits, dag_indices[1], ctrl_state=0)
    cr_circuit.rz(-angle, dag_indices[1])
    cr_circuit.h(dag_indices[1])

    res = change_of_basis_circuit.compose(cr_circuit)
    return res.compose(change_of_basis_circuit.inverse())

def xy_type_circuit_baseline_mixed(
        dag_indices: tuple[int, int],
        ndag_indices: tuple[int, int],
        n_type_bits: list[int],
        num_qubits: int,
        angle: complex,
) -> QuantumCircuit:
    register = QuantumRegister(num_qubits, name="sim")
    change_of_basis_circuit = get_change_of_basis_circuit(dag_indices, ndag_indices, n_type_bits, register)

    alpha = np.real(angle)
    beta = np.imag(angle)

    t = np.sqrt(alpha**2 + beta**2)
    theta = 0.5*np.atan2(alpha,beta)

    cr_circuit = QuantumCircuit(register)
    ctrl_qubits = list(range(ndag_indices[0], dag_indices[1])) + n_type_bits
    cr_circuit.h(dag_indices[1])
    cr_circuit.rx(theta*2, dag_indices[1])
    cr_circuit.mcx(ctrl_qubits, dag_indices[1], ctrl_state=0)
    cr_circuit.rz(t, dag_indices[1])
    cr_circuit.mcx(ctrl_qubits, dag_indices[1], ctrl_state=0)
    cr_circuit.rz(-t, dag_indices[1])
    cr_circuit.rx(-theta*2, dag_indices[1])
    cr_circuit.h(dag_indices[1])

    res = change_of_basis_circuit.compose(cr_circuit)
    return res.compose(change_of_basis_circuit.inverse())

def n_type_circuit(
        n_qubits: list[int],
        num_qubits: int,
        coeff: float,
        angle: float,
) -> QuantumCircuit:
    register = QuantumRegister(num_qubits, name="sim")

    if len(n_qubits) == 0:
        cr_circuit = QuantumCircuit(register, global_phase=coeff * angle)
    elif len(n_qubits) == 1:
        cr_circuit = QuantumCircuit(register, global_phase=coeff * angle / 2)
        cr_circuit.rz(coeff*angle, n_qubits)
    else:
        sign = (-1) ** ((len(n_qubits) * (len(n_qubits) - 1)) // 2)
        clean_ancilla = QuantumRegister(1, name = "c-a")
        cr_circuit = QuantumCircuit(register, clean_ancilla, global_phase=sign*coeff * angle / 2)
        cr_circuit.mcx([register[i] for i in n_qubits], clean_ancilla[0])
        cr_circuit.rz(sign*coeff*angle, clean_ancilla[0])
        cr_circuit.mcx([register[i] for i in n_qubits], clean_ancilla[0])
    return cr_circuit

def compose_fermionic_circuits(
        circuit1: QuantumCircuit,
        circuit2: QuantumCircuit,
        inplace: bool = False,
) -> QuantumCircuit | None:
    qreg_names_1 = set([reg.name for reg in circuit1.qregs])
    qreg_names_2 = set([reg.name for reg in circuit2.qregs])
    if qreg_names_2 <= qreg_names_1:
        return circuit1.compose(circuit2, inplace = inplace)
    else:
        if inplace:
            res = circuit1
        else:
            res = deepcopy(circuit1)

        for qr in circuit2.qregs:
            if qr not in res.qregs:
                res.add_register(qr)
        res.compose(circuit2, inplace = True)
        if not inplace:
            return res

if __name__ == "__main__":
    op = FermionOperator("0^ 1", 1j) + FermionOperator("1^ 0", 1j)
    print(single_normal_ordered_fermion_operator_circuit(op, np.pi/8, 2))