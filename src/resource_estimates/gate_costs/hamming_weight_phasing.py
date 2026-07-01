import math
import numpy as np

from src.resource_estimates.gate_costs.protocol import ResourceGate, gate_labels

"""
Gate method for handling Hamming Weight Phasing.

Warning: Does not check if the circuit can actually be hamming weight-phased using the specified parameters!
"""

def _compute_num_rotations(dct: dict[ResourceGate: int], m: int) -> int:
    NR = sum([n*G.arbitrary_rotations for G,n in dct.items()])
    m_rounds = NR // m
    remaining_round_size = NR % m
    if m_rounds > 0 and remaining_round_size == 0:
        return m_rounds * math.floor(np.log2(m) + 1)
    elif m_rounds > 0:
        return m_rounds * math.floor(np.log2(m) + 1) + math.floor(np.log2(remaining_round_size) + 1)
    else:
        return math.floor(np.log2(remaining_round_size) + 1)
class HWPGate(ResourceGate):

    def __init__(
            self,
            composite_gates: dict[ResourceGate: int],
            m: int
    ) -> None:
        self.composite_gates = composite_gates
        self.m = m
        self._ancilla_qubits = self.m - 1
        self._extra_toffolis = self.m - 1
        self._toffoli = self._extra_toffolis + sum([n*G.toffoli for G,n in composite_gates.items()])
        self._rz = _compute_num_rotations(self.composite_gates, m)
        self._ry = 0
        self._rx = 0

        for gate in gate_labels:
            if gate in ["rz", "rx", "ry", "toffoli"]:
                continue
            gate_count = sum([n*getattr(G, gate) for G,n in composite_gates.items()])
            setattr(self, "_" + gate, gate_count)

        self._h += sum([2*n*(G.rx + G.ry) for G,n in composite_gates.items()])
        self._s += sum([n*G.ry for G,n in composite_gates.items()])
        self._sdg += sum([n*G.ry for G,n in composite_gates.items()])