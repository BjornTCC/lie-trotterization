from src.resource_estimates.gate_costs.protocol import ResourceGate

class Kappa(ResourceGate):

    @property
    def ry(self) -> int:
        return 2

    @property
    def cx(self) -> int:
        return 4

    @property
    def rotation_depth(self) -> int:
        return 2

class Happa(ResourceGate):

    @property
    def rz(self) -> int:
        return 2

    @property
    def h(self) -> int:
        return 2

    @property
    def cx(self) -> int:
        return 4

    @property
    def rotation_depth(self) -> int:
        return 2
class Occupation(ResourceGate):

    @property
    def rz(self) -> int:
        return 1

    @property
    def rotation_depth(self) -> int:
        return 1

class Fij(ResourceGate):

    @property
    def t(self) -> int:
        return 1

    @property
    def tdg(self) -> int:
        return 1

    @property
    def h(self) -> int:
        return 6

    @property
    def s(self) -> int:
        return 3

    @property
    def cx(self) -> int:
        return 3

    @property
    def t_depth(self) -> int:
        return 2

fgate = Fij()

class FreeFermionicS1Tile(ResourceGate):

    @property
    def rz(self) -> int:
        return 2

    @property
    def h(self) -> int:
        return 8

    @property
    def s(self) -> int:
        return 6

    @property
    def cx(self) -> int:
        return 2

    @property
    def rotation_depth(self) -> int:
        return 1

class FreeFermionicS2Tile(ResourceGate):

    @property
    def rz(self) -> int:
        return 2 + 2*fgate.rz

    @property
    def t(self) -> int:
        return 2*fgate.t

    @property
    def tdg(self) -> int:
        return 2*fgate.tdg

    @property
    def h(self) -> int:
        return 8 + 2*fgate.h

    @property
    def s(self) -> int:
        return 6 + 2*fgate.s

    @property
    def cx(self) -> int:
        return 2 + 2*fgate.cx

    @property
    def t_depth(self) -> int:
        return 2*fgate.t_depth

    @property
    def rotation_depth(self) -> int:
        return 1

class FreeFermionicS3Tile(ResourceGate):

    @property
    def rz(self) -> int:
        return 2 + 2 * fgate.t + 2*fgate.tdg + 2*fgate.rz

    @property
    def t(self) -> int:
        return 4 * fgate.t

    @property
    def tdg(self) -> int:
        return 4 * fgate.tdg

    @property
    def h(self) -> int:
        return 8 + 4 * fgate.h

    @property
    def s(self) -> int:
        return 6 + 4 * fgate.s

    @property
    def cx(self) -> int:
        return 2 + 4 * fgate.cx


    @property
    def t_depth(self) -> int:
        return 2 + 2*fgate.t_depth
    @property
    def rotation_depth(self) -> int:
        return 1 + 2*fgate.t_depth


class FreeFermionicS4Tile(ResourceGate):

    @property
    def rz(self) -> int:
        return 2 + 6 * fgate.rz

    @property
    def t(self) -> int:
        return 6 * fgate.t

    @property
    def tdg(self) -> int:
        return 6 * fgate.tdg

    @property
    def h(self) -> int:
        return 8 + 6 * fgate.h

    @property
    def s(self) -> int:
        return 6 + 6 * fgate.s

    @property
    def cx(self) -> int:
        return 2 + 6 * fgate.cx

    @property
    def fswap(self) -> int:
        return 2

    @property
    def t_depth(self) -> int:
        return 4*fgate.t_depth

    @property
    def rotation_depth(self) -> int:
        return 1

class FreeFermionicC4Tile(ResourceGate):

    @property
    def rz(self) -> int:
        return 2 + 4 * fgate.rz

    @property
    def t(self) -> int:
        return 4 * fgate.t

    @property
    def tdg(self) -> int:
        return 4 * fgate.tdg

    @property
    def h(self) -> int:
        return 8 + 4 * fgate.h

    @property
    def s(self) -> int:
        return 6 + 4 * fgate.s

    @property
    def cx(self) -> int:
        return 2 + 4 * fgate.cx

    @property
    def t_depth(self) -> int:
        return 2*fgate.t_depth

    @property
    def rotation_depth(self) -> int:
        return 1