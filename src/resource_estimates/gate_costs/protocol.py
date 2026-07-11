from __future__ import annotations

from abc import ABC

from qiskit import QuantumCircuit, QuantumRegister

gate_labels = [
    "rz", "ry", "rx",
    "t", "tdg", "toffoli",
    "s", "sdg", "h",
    "cx", "fswap",
    "x", "y", "z"
]

depth_labels = [
    "t_depth",
    "toffoli_depth",
    "rotation_depth"
]

class ResourceGate(ABC):

    def __init__(self, *args, **kwargs) -> None:
        self._ancilla_qubits = 0
        for gate in gate_labels:
            setattr(self, "_" + gate, 0)

        for gate in depth_labels:
            setattr(self, "_" + gate, 0)
        self.__post_init__(*args, **kwargs)

    @property
    def ancilla_qubits(self) -> int:
        return self._ancilla_qubits

    @property
    def rz(self) -> int:
        return self._rz

    @property
    def ry(self) -> int:
        return self._ry

    @property
    def rx(self) -> int:
        return self._rx

    @property
    def t(self) -> int:
        return self._t

    @property
    def tdg(self) -> int:
        return self._tdg

    @property
    def toffoli(self) -> int:
        return self._toffoli

    @property
    def cx(self) -> int:
        return self._cx

    @property
    def fswap(self) -> int:
        return self._fswap

    @property
    def h(self) -> int:
        return self._h

    @property
    def s(self) -> int:
        return self._s

    @property
    def sdg(self) -> int:
        return self._sdg

    @property
    def x(self) -> int:
        return self._x

    @property
    def y(self) -> int:
        return self._y

    @property
    def z(self) -> int:
        return self._z

    @property
    def cliffords(self) -> int:
        return self.cx + self.s + self.sdg + self.h + self.x + self.y + self.z

    @property
    def arbitrary_rotations(self) -> int:
        return self.rx + self.ry + self.rz

    @property
    def t_depth(self) -> int:
        return self._t_depth

    @property
    def toffoli_depth(self) -> int:
        return self._toffoli_depth

    @property
    def rotation_depth(self) -> int:
        return self._rotation_depth

    def symmetric_controlled(self) -> ResourceGate:
        raise NotImplementedError(f"Resource gate {self}, does not have symmetric_controlled implemented.")

    def to_qiskit_circuit(self, register: QuantumRegister, qubits: list, parameters: list) -> QuantumCircuit:
        raise NotImplementedError(f"Qiskit circuit for ResourceGate {self} not implemented")

    def gate_costs_as_dict(self) -> dict[str: int]:
        return {x: getattr(self, x) for x in gate_labels if getattr(self, x) > 0}

    def non_clifford_depths(self) -> dict[str: int]:
        return {x: getattr(self, x) for x in depth_labels if getattr(self, x) > 0}