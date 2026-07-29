from src.resource_estimates.gate_costs.protocol import ResourceGate, gate_labels, depth_labels

def _add_dicts_inplace(dct1, dct2, scale2: int = 1):
    for x,y in dct2.items():
        if x in dct1.keys():
            dct1[x] += y*scale2
        else:
            dct1[x] = y*scale2

def compute_FKN_counts(logn: int) -> dict[tuple: int]:
    if logn == 1:
        return {(0,2): 1}

    n = 2**logn
    res = {
        (i, n): 1 for i in range(n // 2)
    }
    rec_res = compute_FKN_counts(logn - 1)
    _add_dicts_inplace(res, rec_res, 2)
    return res

class FFFTOneDimension(ResourceGate):

    def __post_init__(self, logn: int) -> None:
        n = 2**logn
        self._fswap = (n // 2)*(2*n - logn - 3)
        self._FKN_counts = compute_FKN_counts(logn)

        for (k,N), c in self._FKN_counts.items():
            self._h += 6*c
            self._s += 3*c
            self._t += 1*c
            self._tdg += 1*c
            self._cx += 3*c

            if ((4*k) % N) == 0:
                self._s += c * (4*k) // N
            elif ((8*k) % N) == 0:
                self._t += c
            elif ((16*k) % N) == 0:
                self._t += c*5 / 2
            else:
                self._rz += 1
        self._t = int(self._t)

        self._t_depth = 2*logn
        if logn >= 3:
            self._t_depth += 1
            self._rotation_depth = logn - 3

class FFFT(ResourceGate):

    def __post_init__(self, logn: int, d: int, spin: bool) -> None:
        base_application = FFFTOneDimension(logn)
        num_applications = d*(int(spin) + 1) * (2**((d-1)*logn))

        for gate in gate_labels + depth_labels:
            val = num_applications * getattr(base_application, "_" + gate)
            setattr(self, "_" + gate, val)