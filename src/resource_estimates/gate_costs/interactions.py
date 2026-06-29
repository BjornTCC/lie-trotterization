from src.resource_estimates.gate_costs.protocol import ResourceGate

class AncillaOccupationPair(ResourceGate):

    @property
    def ancilla_qubits(self) -> int:
        return 1

    @property
    def rz(self) -> int:
        return 1

    @property
    def toffoli(self) -> int:
        return 2

class AncillaFreeOccupationPair(ResourceGate):

    @property
    def rz(self) -> int:
        return 3

    @property
    def cx(self) -> int:
        return 2


class ShiftedOccupationPair(ResourceGate):

    @property
    def rz(self) -> int:
        return 1

    @property
    def cx(self) -> int:
        return 2