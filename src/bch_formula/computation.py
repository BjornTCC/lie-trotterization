import json

from copy import deepcopy
"""
Computing the symmetric BCH formula using results from the repo:
https://github.com/HaraldHofstaetter/BCH?tab=readme-ov-file
"""

def table_to_json(table_file: str) -> None:

    res = []

    with open(table_file) as file:
        next(file)
        for line in file:
            data = line.rstrip().split("\t")

            index = int(data[0])
            order = int(data[1])
            i1 = int(data[2])
            i2 = int(data[3])

            c = data[4].split("/")
            coeff = int(c[0]), int(c[1])

            res.append({
                "order": order,
                "i1": i1,
                "i2": i2,
                "coeff": tuple(coeff)
            })

    with open(table_file.split(".")[0] + ".json", "w") as f:
        json.dump(res,f, indent = 4)

def compute_bch_terms(
        sum: list[any],
        commutator_func: any,
        max_term: int = 3,
        bch_terms: dict = None,
) -> list[any]:

    if bch_terms is None:
        with open("_symmetric_bch.json", "r") as f:
            bch_terms = json.load(f)
    if max_term == 3:
        return [compute_bch_terms_3(sum, commutator_func)]
    else:
        raise ValueError("not implemented yet")

def compute_bch_terms_3(
        sum: list[any],
        commutator_func: any,
) -> any:
    if len(sum) == 1:
        return 0.0

    A = sum[0]
    B = deepcopy(sum[1])
    for t in sum[2:]:
        B += t

    C = commutator_func(A,B)

    T1 = -commutator_func(A,C) * 1/24
    T2 = -commutator_func(B,C) * 1/12

    return T1 + T2 + compute_bch_terms_3(sum[1:], commutator_func)

def compute_bch_terms_rec(
        sum: list[any],
        terms: list[any],
        commutator_func: any,
        max_term: int,
        bch_terms: dict,
) -> list[list[any]]:
    ...
