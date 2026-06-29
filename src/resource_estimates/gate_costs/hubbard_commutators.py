from src.resource_estimates.gate_costs.protocol import ResourceGate

class MixedControlledKappa(ResourceGate):

    @property
    def ry(self) -> int:
        return 2

    @property
    def toffoli(self) -> int:
        return 4

    @property
    def cx(self) -> int:
        return 2

class MixedControlledHappa(ResourceGate):

    @property
    def rz(self) -> int:
        return 2

    @property
    def h(self) -> int:
        return 2

    @property
    def toffoli(self) -> int:
        return 4

    @property
    def cx(self) -> int:
        return 2


class SpinSymmetricMixedControlledKappa(ResourceGate):

    @property
    def ancilla_qubits(self) -> int:
        return 1

    @property
    def x(self) -> int:
        return 2

    @property
    def h(self) -> int:
        return 4 * 2

    @property
    def s(self) -> int:
        return 4 * 1

    @property
    def sdg(self) -> int:
        return 4 * 1

    @property
    def ry(self) -> int:
        return 2

    @property
    def tdg(self) -> int:
        return 4*1

    @property
    def t(self) -> int:
        return 4*1

    @property
    def toffoli(self) -> int:
        return 4

    @property
    def cx(self) -> int:
        return 6 + 4*4

    @property
    def fswap(self) -> int:
        return 2

class SpinSymmetricMixedControlledHappa(ResourceGate):

    @property
    def ancilla_qubits(self) -> int:
        return 1

    @property
    def x(self) -> int:
        return 2

    @property
    def h(self) -> int:
        return 4 * 2

    @property
    def s(self) -> int:
        return 4 * 1

    @property
    def sdg(self) -> int:
        return 4 * 1

    @property
    def ry(self) -> int:
        return 2

    @property
    def tdg(self) -> int:
        return 4 * 1

    @property
    def t(self) -> int:
        return 4 * 1

    @property
    def toffoli(self) -> int:
        return 4

    @property
    def cx(self) -> int:
        return 6 + 4 * 4
