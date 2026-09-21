from src.resource_estimates.gate_costs.protocol import ResourceGate
from src.resource_estimates.gate_costs.pauli_string_rotation import PauliStringGate

class XField(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_xs = 1)

class YField(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_ys=1)

class ZField(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_zs = 1)


class ZInteraction(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_zs=2)

class YInteraction(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_ys=2)

class ZYInteraction(ResourceGate):

    def __post_init__(self) -> None:
        self._cx = 2
        self._rz = 1
        self._ry = 1
        self._h = 2
        self._s = 1
        self._sdg = 1