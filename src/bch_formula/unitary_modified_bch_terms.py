def unitary_modified_bch_3_v1(
        terms: list[any],
        commutator_func: any
) -> tuple[any, any]:
    # Returns in format (3'rd order terms, unitary generator)
    M = len(terms)

    unitary_generator = 0
    for i in range(M):
        for j in range(i+1, M):
            unitary_generator += commutator_func(terms[j], terms[i]) * 1/24

    mid_term = 0
    for k in range(M):
        for j in range(k+1, M):
            mid_term += commutator_func(terms[j], commutator_func(terms[j], terms[k])) * 1/24

    for k in range(M):
        for j in range(k+1, M):
            for i in range(j+1, M):
                mid_term += commutator_func(terms[i], commutator_func(terms[j], terms[k])) * 1/12
    return mid_term, unitary_generator