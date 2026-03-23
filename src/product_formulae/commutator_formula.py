from src.bch_formula import compute_bch_terms

from src.product_formulae.suzuki_trotter import suzuki_trotter

def commutator_formula(
        H: any,
        t: float,
        commutator_func: any,
        decomposer: any,
        order: int = 3,
) -> list[any]:
    terms = decomposer(H)
    if order == 1:
        return suzuki_trotter(terms, t)

    bch_terms = decomposer(compute_bch_terms(terms, commutator_func)[0])
    linear_terms = [t*x/2 for x in terms]
    cubed_terms = [-t**3 * x for x in bch_terms]
    return linear_terms + cubed_terms + linear_terms[::-1]