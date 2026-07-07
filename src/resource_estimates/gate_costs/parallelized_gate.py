from src.resource_estimates.gate_costs.protocol import ResourceGate, gate_labels, depth_labels

class ParallelizedResourceGate(ResourceGate):

    def __init__(
            self,
            gates: dict[ResourceGate: int],
    ) -> None:
        self.gates = gates
        for gate in gate_labels:
            attr = "_" + gate
            val = sum({v * getattr(g, gate) for g,v in gates.items()})
            setattr(self, attr, val)

        for gate_depth in depth_labels:
            val = max([getattr(g, gate_depth) for g in gates.keys()])
            setattr(self, "_" + gate_depth, val)