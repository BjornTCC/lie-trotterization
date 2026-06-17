

def approximate_disentangling_formula(
        terms: list[any],
        t: float,
        commutator_func: any,
        decomposer: any,
) -> list[any]:
    res = []
    res_rev = []

    M = len(terms)

    outer_unitary = 0
    for j in range(1, M):
        outer_unitary += t**2*commutator_func(terms[j], terms[0])/24

    outer_left = decomposer(outer_unitary)
    outer_right = [-x for x in outer_left[::-1]]

    left_res = []
    right_res = []

    for j in range(2*M-2):
        if j % 2 == 0:
          left_res.append(t*terms[j // 2]/2)
          right_res.append(t*terms[j // 2]/2)
        else:
            c = 0
            for i in range((j-1)//2, M):
                c += t**2*commutator_func(terms[i], terms[(j-1) // 2])/24
            for i in range((j + 1) // 2, M):
                c += t ** 2 * commutator_func(terms[i], terms[(j + 1) // 2]) / 24

            cterms = decomposer(c)
            left_res.extend(cterms)
            right_res.extend([-x for x in cterms])

    return outer_left + left_res+ [t*terms[-1]] + right_res[::-1] + outer_right