from copy import deepcopy

import numpy as np

from qiskit import QuantumCircuit, QuantumRegister

from openfermion.ops import FermionOperator

def fswap(qc: QuantumCircuit, i: int, j: int, inplace: bool = True) -> QuantumCircuit | None:
    if inplace:
        qc.cx(i,j)
        qc.cz(i,j)
    else:
        qc_new = qc.copy()
        fswap(qc_new, i,j)
        return qc_new

def single_no_fermion_operator_circuit(
        operator: FermionOperator,
        angle: float,
        num_qubits: int,
        real: bool
) -> QuantumCircuit:
    assert len(operator.terms) == 1, "Operator must contain a single term"
    assert operator.is_normal_ordered(), "Operator must be normal ordered"

    ads, coeff = list(operator.terms.items())[0]

    ctrl_up = set([x[0] for x in ads if x[1]])
    ctrl_down = set([x[0] for x in ads if not x[1]])

    if ctrl_up == ctrl_down:
        return n_type_circuit(list(ctrl_up), num_qubits, np.imag(coeff), angle)

    n_bits = ctrl_up.intersection(ctrl_down)

    return xy_type_circuit(
        num_qubits,
        ctrl_up - n_bits,
        ctrl_down - n_bits,
        n_bits,
        coeff,
        angle,
        real
    ).reverse_bits()

def xy_type_circuit(
        num_qubits: int,
        ctrl_up: set[int],
        ctrl_down: set[int],
        n_bits: set[int],
        coeff: float,
        angle: float,
        real: bool,
) -> QuantumCircuit:
    # Needs to be adjusted for phase on passive qubits
    if len(ctrl_up) == 0:
        return xy_type_circuit(
            num_qubits,
            ctrl_down,
            ctrl_up,
            n_bits,
            -coeff if real else coeff,
            angle,
            real
        )
    else:
        target_qubit = list(ctrl_up)[0]
        ctrl_up_mod = ctrl_up - set([target_qubit])

        # Compute the qubits that introduce a phase

        z_intervals = sorted(list(ctrl_up) + list(ctrl_down))

        phase_qubits = []
        for k in range(len(z_intervals) - 1):
            if (k % 2) == 0:
                phase_qubits.extend([i for i in range(z_intervals[k] + 1, z_intervals[k+1])])

        register = QuantumRegister(num_qubits, name="sim")
        cr_circuit = QuantumCircuit(register)
        ctrl_list = list(ctrl_up_mod) + list(ctrl_down) + list(n_bits)
        if real:
            if len(ctrl_list) > 0:
                cr_circuit.mcx(ctrl_list, target_qubit)
                cr_circuit.ry(-np.real(coeff) * angle, target_qubit)
                cr_circuit.mcx(ctrl_list, target_qubit)
                cr_circuit.ry(np.real(coeff) * angle, target_qubit)
            else:
                cr_circuit.ry(2*np.real(coeff) * angle, target_qubit)
        else:
            if len(ctrl_list) > 0:
                cr_circuit.mcx(ctrl_list, target_qubit)
                cr_circuit.rz(np.imag(coeff) * angle, target_qubit)
                cr_circuit.mcx(ctrl_list, target_qubit)
                cr_circuit.rz(-np.imag(coeff) * angle, target_qubit)
            else:
                cr_circuit.rz(-2*np.imag(coeff) * angle, target_qubit)

        change_of_basis_circuit = QuantumCircuit(register)
        for i in ctrl_up_mod:
            change_of_basis_circuit.cx(target_qubit, i)
            change_of_basis_circuit.x(i)

        for i in ctrl_down:
            change_of_basis_circuit.cx(target_qubit, i)
        if not real:
            change_of_basis_circuit.h(target_qubit)

        for i in phase_qubits:
            change_of_basis_circuit.cx(i, target_qubit)


        res = change_of_basis_circuit.compose(cr_circuit)
        res.compose(change_of_basis_circuit.inverse(), inplace = True)
        return res

def n_type_circuit(
        n_qubits: list[int],
        num_qubits: int,
        coeff: float,
        angle: float,
) -> QuantumCircuit:
    register = QuantumRegister(num_qubits, name="sim")

    if len(n_qubits) == 1:
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

    from qiskit.quantum_info import Operator

    operator = FermionOperator("1^ 1", 1)
    operator2 = FermionOperator("2^ 1", 1)
    nq = 3
    theta = 0.1

    qc1 = single_no_fermion_operator_circuit(operator, theta, nq, True)
    qc2 = single_no_fermion_operator_circuit(operator2, theta, nq, True)
    print(qc1)
    print(qc2)

    print(dir(qc1.qubits[0]._register))
    print(qc1.qubits[0]._register.name)
    print([q._index for q in qc1.qubits if q._register.name == "sim"])
    qc = compose_fermionic_circuits(qc2,qc1, inplace=True)
    qc = qc2
    print(qc)

    print(f"(cos, sin) = ({np.cos(theta), np.sin(theta)})")
    print(f"e^i*theta/2 = {np.exp(1j*theta/2)}")

    mat = Operator(qc).to_matrix()
    print("Non-zero matrix elements:")
    for i in range(2**nq):
        for j in range(2**nq):
            if abs(mat[i,j]) > 1e-10:
                print(format(i, f'0{nq}b') +"," + format(j, f'0{nq}b') + f": {mat[i,j]}")


    Operator(qc).to_matrix()