from src.resource_estimates.gate_costs.protocol import ResourceGate

class AncillaOccupationPair(ResourceGate):

    def __post_init__(self) -> None:
        self._ancilla_qubits = 1
        self._rz = 1
        self._toffoli = 2

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledAncillaOccupationPair()

class ControlledAncillaOccupationPair(ResourceGate):

    def __post_init__(self) -> None:
        self._ancilla_qubits = 1
        self._rz = 1
        self._cx = 2
        self._toffoli = 2
class AncillaFreeOccupationPair(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 3
        self._cx = 2

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledAncillaFreeOccupationPair()
class ControlledAncillaFreeOccupationPair(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 3
        self._cx = 4

class ShiftedOccupationPair(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 1
        self._cx = 2

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledShiftedOccupationPair()
class ControlledShiftedOccupationPair(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 1
        self._cx = 4