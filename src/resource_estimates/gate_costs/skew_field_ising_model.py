from src.resource_estimates.gate_costs.protocol import ResourceGate
from src.resource_estimates.gate_costs.pauli_string_rotation import PauliStringGate

class XField(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_xs = 1)

    def symmetric_controlled(self) -> ResourceGate:
        return SymmetricControlledXField()

class SymmetricControlledXField(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_xs = 1, num_zs = 1)

class YField(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_ys=1)

    def symmetric_controlled(self) -> ResourceGate:
        return SymmetricControlledYField()

class SymmetricControlledYField(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_ys = 1, num_zs = 1)

class ZField(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_zs = 1)


    def symmetric_controlled(self) -> ResourceGate:
        return SymmetricControlledZField()

class SymmetricControlledZField(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_zs = 2)

class ZInteraction(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_zs=2)

    def symmetric_controlled(self) -> ResourceGate:
        return SymmetricControlledZInteraction()

class SymmetricControlledZInteraction(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_zs = 3)

class YInteraction(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_ys=2)

    def symmetric_controlled(self) -> ResourceGate:
        return SymmetricControlledYInteraction()

class SymmetricControlledYInteraction(PauliStringGate):

    def __init__(self) -> None:
        super().__init__(num_ys = 2, num_zs = 1)

class ZYInteraction(ResourceGate):

    def __post_init__(self) -> None:
        self._cx = 2
        self._rz = 1
        self._ry = 1
        self._h = 2
        self._s = 1
        self._sdg = 1

    def symmetric_controlled(self) -> ResourceGate:
        return SymmetricControlledZYInteraction()

class SymmetricControlledZYInteraction(ResourceGate):

    def __post_init__(self) -> None:
        self._cx = 4
        self._rz = 1
        self._ry = 1
        self._h = 2
        self._s = 1
        self._sdg = 1

