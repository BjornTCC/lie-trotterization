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
) -> dict[Fswap: int]:
    permutation = {i: j for i,j in enumerate(permutation)}

    res = 0
    _permutation = deepcopy(permutation)
    N = len(_permutation)
    for i in range(N - 1):
        j = i
        if j < N - 1:
            while _permutation[j + 1] < _permutation[j]:
                res += 1
                _permutation[j], _permutation[j + 1] = _permutation[j + 1], _permutation[j]
                if j == N - 2:
                    break
                else:
                    j += 1
        j = i
        if j > 0:
            while _permutation[j - 1] > _permutation[j]:
                res += 1
                _permutation[j], _permutation[j - 1] = _permutation[j - 1], _permutation[j]
                if j == 1:
                    break
                else:
                    j -= 1
    return res
