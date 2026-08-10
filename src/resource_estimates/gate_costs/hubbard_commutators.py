from src.resource_estimates.gate_costs.protocol import ResourceGate

class MixedControlledKappa(ResourceGate):

    def __post_init__(self) -> None:
        self._ry = 2
        self._toffoli = 4
        self._cx = 2

        self._rotation_depth = 2
        self._toffoli_depth = 4

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledMixedControlledKappa()

class ControlledMixedControlledKappa(ResourceGate):

    def __post_init__(self) -> None:
        self._ry = 2
        self._toffoli = 4
        self._cx = 4

        self._rotation_depth = 2
        self._toffoli_depth = 4


class MixedControlledHappa(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._h = 2
        self._toffoli = 4
        self._cx = 2

        self._rotation_depth = 2
        self._toffoli_depth = 4

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledMixedControlledHappa()

class ControlledMixedControlledHappa(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._h = 2
        self._toffoli = 4
        self._cx = 4

        self._rotation_depth = 2
        self._toffoli_depth = 4

class SpinSymmetricMixedControlledKappa(ResourceGate):

    def __post_init__(self) -> None:
        self._ancilla_qubits = 2
        self._x = 4
        self._h = 6
        self._ry = 2
        self._tdg = 0
        self._t = 8
        self._toffoli = 0
        self._cx = 22

        self._rotation_depth = 1
        self._t_depth = 3


    def symmetric_controlled(self) -> ResourceGate:
        return ControlledSpinSymmetricMixedControlledKappa()

class ControlledSpinSymmetricMixedControlledKappa(ResourceGate):

    def __post_init__(self) -> None:
        self._ancilla_qubits = 2
        self._x = 4
        self._h = 6
        self._ry = 2
        self._tdg = 0
        self._t = 8
        self._toffoli = 0
        self._cx = 24

        self._rotation_depth = 1
        self._t_depth = 3

class SpinSymmetricMixedControlledHappa(ResourceGate):

    def __post_init__(self) -> None:
        self._ancilla_qubits = 2
        self._x = 4
        self._h = 4
        self._rz = 2
        self._tdg = 0
        self._t = 8
        self._cx = 22

        self._rotation_depth = 1
        self._t_depth = 3

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledMixedControlledHappa()

class ControlledSpinSymmetricMixedControlledHappa(ResourceGate):

    def __post_init__(self) -> None:
        self._ancilla_qubits = 2
        self._x = 4
        self._h = 4
        self._rz = 2
        self._tdg = 0
        self._t = 8
        self._cx = 24

        self._rotation_depth = 1
        self._t_depth = 3