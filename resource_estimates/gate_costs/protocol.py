
from abc import ABC

from qiskit import QuantumCircuit, QuantumRegister

class ResourceGate(ABC):

    @property
    def ancilla_qubits(self) -> int | None:
        return 0

    @property
    def rotations(self) -> int | None:
        return 0

    @property
    def t(self) -> int:
        return 0

    @property
    def tdg(self) -> int:
        return 0

    @property
    def toffoli(self) -> int:
        return 0

    @property
    def cx(self) -> int | None:
        return 0

    @property
    def h(self) -> int | None:
        return 0

    @property
    def s(self) -> int | None:
        return 0

    @property
    def sdag(self) -> int | None:
        return 0

    @property
    def x(self) -> int | None:
        return 0

    @property
    def y(self) -> int | None:
        return 0

    @property
    def z(self) -> int:
        return 0

    @property
    def cliffords(self) -> int:
        return self.cx + self.s + self.sdg + self.h + self.x + self.y + self.z

    def to_qiskit_circuit(self, register: QuantumRegister, qubits: list, parameters: list) -> QuantumCircuit:
        raise NotImplementedError(f"Qiskit circuit for ResourceGate {self} not implemented")