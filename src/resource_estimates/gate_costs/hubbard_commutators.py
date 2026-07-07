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
    def x(self) -> int:
        return 4

    @property
    def h(self) -> int:
        return 6

    @property
    def ry(self) -> int:
        return 2

    @property
    def tdg(self) -> int:
        return 8

    @property
    def t(self) -> int:
        return 8

    @property
    def toffoli(self) -> int:
        return 0

    @property
    def cx(self) -> int:
        return 22


class SpinSymmetricMixedControlledHappa(ResourceGate):

    @property
    def x(self) -> int:
        return 4

    @property
    def h(self) -> int:
        return 4

    @property
    def rz(self) -> int:
        return 2

    @property
    def tdg(self) -> int:
        return 8

    @property
    def t(self) -> int:
        return 8

    @property
    def toffoli(self) -> int:
        return 0

    @property
    def cx(self) -> int:
        return 22
