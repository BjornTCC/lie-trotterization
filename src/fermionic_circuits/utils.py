from openfermion import FermionOperator, normal_ordered

def check_operators_proportional(A: FermionOperator, B: FermionOperator) -> bool:
    if A == FermionOperator() or B == FermionOperator:
        return A.isclose(B)

    A_norm = normal_ordered(A)
    A_norm /= A_norm.induced_norm()
    B_norm = normal_ordered(B)
    B_norm /= B_norm.induced_norm()

    return (A_norm.isclose(B_norm)) or (A_norm.isclose(-B_norm))
