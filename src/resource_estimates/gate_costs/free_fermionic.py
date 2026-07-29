from src.resource_estimates.gate_costs.protocol import ResourceGate

class Kappa(ResourceGate):

    def __post_init__(self) -> None:
        self._ry = 2
        self._cx = 4
        self._rotation_depth = 2


    def symmetric_controlled(self) -> ResourceGate:
        return ControlledKappa()

class ControlledKappa(ResourceGate):

    def __post_init__(self) -> None:
        self._ry = 2
        self._cx = 6
        self._rotation_depth = 2

class Happa(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._h = 2
        self._cx = 4
        self._rotation_depth = 2

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledHappa()

class ControlledHappa(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._h = 2
        self._cx = 6
        self._rotation_depth = 2

class Occupation(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 1
        self._rotation_depth = 1

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledOccupation()

class ControlledOccupation(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 1
        self._cx = 2
        self._rotation_depth = 1

class Fij(ResourceGate):

    def __post_init__(self) -> None:
        self._t = 1
        self._tdg = 1
        self._h = 1
        self._s = 3
        self._cx = 3
        self._t_depth = 2

fgate = Fij()

class FreeFermionicS1Tile(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._h = 8
        self._s = 6
        self._cx = 2
        self._rotation_depth = 1

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledFreeFermionicS1Tile()

class ControlledFreeFermionicS1Tile(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._h = 8
        self._s = 6
        self._cx = 6
        self._rotation_depth = 1

class FreeFermionicS2Tile(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._t = 2 * fgate.t
        self._tdg = 2 * fgate.tdg
        self._h = 8 + 2*fgate.h
        self._s = 6 + 2*fgate.s
        self._cx = 2 + 2*fgate.cx

        self._t_depth = 2*fgate.t_depth
        self._rotation_depth = 1

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledFreeFermionicS2Tile()

class ControlledFreeFermionicS2Tile(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._t = 2 * fgate.t
        self._tdg = 2 * fgate.tdg
        self._h = 8 + 2 * fgate.h
        self._s = 6 + 2 * fgate.s
        self._cx = 6 + 2 * fgate.cx

        self._t_depth = 2 * fgate.t_depth
        self._rotation_depth = 1

class FreeFermionicS3Tile(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2 + 2 * fgate.t + 2*fgate.tdg
        self._t = 4*fgate.t
        self._tdg = 4* fgate.tdg
        self._h = 8 + 4*fgate.h
        self._s = 6 + 4*fgate.s
        self._cx = 2 + 4*fgate.cx

        self._t_depth = 2 + 2*fgate.t_depth
        self._rotation_depth = 1 + 2*fgate.t_depth

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledFreeFermionicS3Tile()

class ControlledFreeFermionicS3Tile(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2 + 2 * fgate.t + 2 * fgate.tdg
        self._t = 4 * fgate.t
        self._tdg = 4 * fgate.tdg
        self._h = 8 + 4 * fgate.h
        self._s = 6 + 4 * fgate.s
        self._cx = 6 + 4 * fgate.cx

        self._t_depth = 2 + 2 * fgate.t_depth
        self._rotation_depth = 1 + 2 * fgate.t_depth

class FreeFermionicS4Tile(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._t = 6*fgate.t
        self._tdg = 6*fgate.tdg
        self._h = 8+ 6*fgate.h
        self._s = 6 + 6 * fgate.s
        self._cx = 2 + 6*fgate.cx
        self._fswap = 2

        self._t_depth = 4*fgate.t_depth
        self._rotation_depth = 1

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledFreeFermionicS4Tile()

class ControlledFreeFermionicS4Tile(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._t = 6 * fgate.t
        self._tdg = 6 * fgate.tdg
        self._h = 8 + 6 * fgate.h
        self._s = 6 + 6 * fgate.s
        self._cx = 6 + 6 * fgate.cx
        self._fswap = 2

        self._t_depth = 4 * fgate.t_depth
        self._rotation_depth = 1

class FreeFermionicC4Tile(ResourceGate):
    def __post_init__(self) -> None:
        self._rz = 2
        self._t = 4*fgate.t
        self._tdg = 4*fgate.tdg
        self._h = 8 + 4*fgate.h
        self._s = 6 + 4*fgate.s
        self._cx = 2 + 4* fgate.cx

        self._t_depth = 2*fgate.t_depth
        self._rotation_depth = 1

    def symmetric_controlled(self) -> ResourceGate:
        return ControlledFreeFermionicC4Tile()

class ControlledFreeFermionicC4Tile(ResourceGate):

    def __post_init__(self) -> None:
        self._rz = 2
        self._t = 4 * fgate.t
        self._tdg = 4 * fgate.tdg
        self._h = 8 + 4 * fgate.h
        self._s = 6 + 4 * fgate.s
        self._cx = 2 + 4 * fgate.cx

        self._t_depth = 2 * fgate.t_depth
        self._rotation_depth = 1