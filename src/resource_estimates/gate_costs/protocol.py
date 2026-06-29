
from abc import ABC

from qiskit import QuantumCircuit, QuantumRegister

gate_labels = [
    "rz", "ry", "rx",
    "t", "tdg", "toffoli",
    "s", "sdg", "h",
    "cx", "fswap",
    "x", "y", "z"
]

class ResourceGate(ABC):

    @property
    def ancilla_qubits(self) -> int:
        return 0

    @property
    def rz(self) -> int:
        return 0

    @property
    def ry(self) -> int:
        return 0

    @property
    def rx(self) -> int:
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
    def cx(self) -> int:
        return 0

    @property
    def fswap(self) -> int:
        return 0

    @property
    def h(self) -> int:
        return 0

    @property
    def s(self) -> int:
        return 0

    @property
    def sdg(self) -> int:
        return 0

    @property
    def x(self) -> int:
        return 0

    @property
    def y(self) -> int:
        return 0

    @property
    def z(self) -> int:
        return 0

    @property
    def cliffords(self) -> int:
        return self.cx + self.s + self.sdg + self.h + self.x + self.y + self.z

    @property
    def arbitrary_rotations(self) -> int:
        return self.rx + self.ry + self.rz

    def to_qiskit_circuit(self, register: QuantumRegister, qubits: list, parameters: list) -> QuantumCircuit:
        raise NotImplementedError(f"Qiskit circuit for ResourceGate {self} not implemented")

    def gate_costs_as_dict(self) -> dict[str: int]:
        return {x: getattr(self, x) for x in gate_labels if getattr(self, x) > 0}