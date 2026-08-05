from src.resource_estimates.gate_costs.protocol import ResourceGate

class PauliStringGate:

    def __post_init__(
            self,
            num_xs: int = 0,
            num_ys: int = 0,
            num_zs: int = 0,
            width: int = 1,
    ) -> None:
        if (num_xs + num_ys + num_zs) > 0:
            width = max(num_xs + num_ys + num_zs, width)
        self._cx = 2 * (width - 1)

        self._rz = 1
        self._h += 2*num_xs
        self._s += num_ys
        self._sdg += num_ys
