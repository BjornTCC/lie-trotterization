from copy import deepcopy

from src.resource_estimates.gate_costs.protocol import ResourceGate

class Fswap(ResourceGate):

    def __post_init__(self) -> None:

        self._cx = 2
        self._h = 8
        self._s = 2
        self._sdg = 6


def fswap_network_count(
        permutation: list[int]
) -> int:
    res = 0
    _permutation = deepcopy(permutation)
    N = len(_permutation)
    while True:
        swapped = False
        for i in range(1,N):
            if _permutation[i-1] > _permutation[i]:
                _permutation[i], _permutation[i - 1] = _permutation[i - 1], _permutation[i]
                res += 1
                swapped = True
        if not swapped:
            break
    return res