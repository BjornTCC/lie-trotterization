from copy import deepcopy

from src.bch_formula import compute_bch_terms, unitary_modified_bch_3_v1

from src.product_formulae.suzuki_trotter import suzuki_trotter

def combine_cubic_into_linear(
        linear_terms: list[any],
        cubic_terms: list[any],
        t: float,
        check_combinable: callable,
) -> tuple[list[any], list[any]]:
    new_linear_terms = deepcopy(linear_terms)
    new_cubic_terms = []

    for cubic_term in cubic_terms:
        combined = False
        for i,linear_term in enumerate(linear_terms):
            if check_combinable(cubic_term, linear_term):
                new_linear_terms[i] = linear_term + cubic_term * t**3/2
                combined = True
                break

        if not combined:
            new_cubic_terms.append(cubic_term*t**3)

    return new_linear_terms, new_cubic_terms

def commutator_formula(
        H: any,
        t: float,
        commutator_func: any,
        decomposer: any,
        order: int = 3,
        check_combinable: callable = None,
) -> list[any]:
    terms = decomposer(H)
    if order == 1:
        return suzuki_trotter(terms, t)

    bch_terms = compute_bch_terms(terms, commutator_func, max_term=order)
    linear_terms = [t*x/2 for x in terms]
    if order == 3:
        if check_combinable is not None:
            linear_terms, cubic_terms = combine_cubic_into_linear(linear_terms, decomposer(bch_terms[0]), t, check_combinable)
            return linear_terms + cubic_terms + linear_terms[::-1]
        cubed_terms = [t**3 * x for x in decomposer(bch_terms[0])]
        return linear_terms + cubed_terms + linear_terms[::-1]

    terms5 = [t**5 * x for x in decomposer(bch_terms[1])]
    cubed_terms = [t**3 * x/2 for x in decomposer(bch_terms[0])]
    return linear_terms + cubed_terms + terms5 + cubed_terms[::-1] + linear_terms[::-1]

def suzuki_trotter_with_unitary_correction(
        terms: list[any],
        t: float,
        commutator_func: any,
        decomposer: any,
) -> tuple[list[any], list[any], list[any]]:
    # Returns on the format (Middle terms, U, U^dag)

    linear_terms = [t * x / 2 for x in terms]
    _, unitary_generator = unitary_modified_bch_3_v1(terms, commutator_func)
    middle_terms = linear_terms[:-1] + [linear_terms[-1]] + linear_terms[-1::-1]

    unitary_decomposition = decomposer(unitary_generator)
    L_unitary = [t ** 2 * x for x in unitary_decomposition]
    R_unitary = [-t ** 2 * x for x in unitary_decomposition[::-1]]
    return middle_terms, L_unitary, R_unitary


def commutator_formula_with_unitary_correction(
        H: any,
        t: float,
        commutator_func: any,
        decomposer: any,
        check_combinable: callable = None,
) -> tuple[list[any], list[any], list[any]]:
    # Returns on the format (Middle terms, U, U^dag)

    terms = decomposer(H)
    linear_terms = [t*x/2 for x in terms]
    cubic_terms, unitary_generator = unitary_modified_bch_3_v1(terms, commutator_func)

    unitary_decomposition = decomposer(unitary_generator)
    L_unitary = [t**2 * x for x in unitary_decomposition]
    R_unitary = [-t ** 2 * x for x in unitary_decomposition[::-1]]
    if check_combinable is not None:
        linear_terms, cubic_terms = combine_cubic_into_linear(linear_terms, decomposer(-cubic_terms), t, check_combinable)
        middle_terms = linear_terms + cubic_terms + linear_terms[::-1]
    else:
        middle_terms = linear_terms + [-t**3 * x for x in cubic_terms] + linear_terms[::-1]
    return middle_terms, L_unitary, R_unitary
