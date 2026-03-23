"""
Collection of standard suzuki-trotter formulae
"""

def suzuki_trotter(
        terms: list[any],
        t: float,
        recursion_order: int = 1,
        bounded_parameter: bool = True
) -> list[any]:
    if recursion_order == 1:
        return [t/2 * x for x in terms[:-1]] + [terms[-1] * t] + [t/2 * x for x in terms[:-1][::-1]]

    elif bounded_parameter:
        pm = (4 - 4**(1 / (2*recursion_order - 1)))**(-1)
        rec1 = suzuki_trotter(terms, pm*t, recursion_order-1, bounded_parameter=bounded_parameter)
        rec2 = suzuki_trotter(terms, (1-4*pm)*t, recursion_order-1, bounded_parameter=bounded_parameter)
        return rec1 * 2 + rec2 + rec1*2
    else:
        pm = (2 - 2**(1 / (2*recursion_order - 1)))**(-1)
        rec1 = suzuki_trotter(terms, pm*t, recursion_order-1, bounded_parameter=bounded_parameter)
        rec2 = suzuki_trotter(terms, (1-2*pm)*t, recursion_order-1, bounded_parameter=bounded_parameter)
        return rec1  + rec2 + rec1
