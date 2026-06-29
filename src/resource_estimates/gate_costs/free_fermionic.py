from src.resource_estimates.gate_costs.protocol import ResourceGate

class SingleExcitation(ResourceGate):

    @property
    def ry(self) -> int:
        return 2

    @property
    def cx(self) -> int:
        return 4

class SingleHermitianExcitation(ResourceGate):

    @property
    def rz(self) -> int:
        return 2

    @property
    def h(self) -> int:
        return 2

    @property
    def cx(self) -> int:
        return 4

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

class FreeFermionicS2Tile(ResourceGate):

    @property
    def rz(self) -> int:
        return 2

    @property
    def t(self) -> int:
        return 2


    @property
    def tdg(self) -> int:
        return 2

    @property
    def h(self) -> int:
        return 8 + 2*6

    @property
    def s(self) -> int:
        return 6 + 2*3

    @property
    def cx(self) -> int:
        return 2 + 2*3

    class FreeFermionicS3Tile(ResourceGate):

        @property
        def rz(self) -> int:
            return 2 + 2*2

        @property
        def t(self) -> int:
            return 2

        @property
        def tdg(self) -> int:
            return 2

        @property
        def h(self) -> int:
            return 8 + 4 * 6

        @property
        def s(self) -> int:
            return 6 + 4 * 3

        @property
        def cx(self) -> int:
            return 2 + 4 * 3

    class FreeFermionicS4Tile(ResourceGate):

        @property
        def rz(self) -> int:
            return 2

        @property
        def t(self) -> int:
            return 1*6

        @property
        def tdg(self) -> int:
            return 1*6

        @property
        def h(self) -> int:
            return 8 + 6 * 6

        @property
        def s(self) -> int:
            return 6 + 6 * 3

        @property
        def cx(self) -> int:
            return 2 + 6 * 3

        @property
        def fswap(self) -> int:
            return 2

    class FreeFermionicC4Tile(ResourceGate):

        @property
        def rz(self) -> int:
            return 2

        @property
        def t(self) -> int:
            return 4*1

        @property
        def tdg(self) -> int:
            return 4*1

        @property
        def h(self) -> int:
            return 8 + 4 * 6

        @property
        def s(self) -> int:
            return 6 + 4 * 3

        @property
        def cx(self) -> int:
            return 2 + 4 * 3