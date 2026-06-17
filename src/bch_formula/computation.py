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
        with open("src/bch_formula/_symmetric_bch.json", "r") as f:
            bch_terms = json.load(f)
    if max_term == 3:
        return [compute_bch_terms_3(sum, commutator_func)]
    elif max_term == 5:
        return [compute_bch_terms_3(sum, commutator_func), compute_bch_terms_5(sum, commutator_func)]
    else:
        raise ValueError("not implemented yet")

def compute_bch_terms_3(
        terms: list[any],
        commutator_func: any,
) -> any:
    c = commutator_func

    res = 0.0

    for i in range(len(terms)):
        A = -terms[i]
        B = sum(terms[i:])
        res -= c(A,c(A,B)) / 24 + c(B,c(A,B)) / 12

    #print("\ncub:",res)
    return res

def compute_bch_terms_5(
        terms: list[any],
        commutator_func: any,
) -> any:
    c = commutator_func

    res = 0.0

    for i in range(len(terms)):
        A = -terms[i]
        B = sum(terms[i:])
        #print("\nB:")
        #print(B)

        cubic_term = compute_bch_terms_3(terms[i+1:], c)

        #print(f"\ncubic_term:")
        #print(cubic_term)

        T0 = c(A,c(A,cubic_term)) * (-1/24)

        T1 = c(A,c(A,c(A,c(A,B)))) * (7.0/5760)
        T2 = c(A,c(A,c(c(A,B),B))) * (-7.0/1440)
        T3 = c(c(A,c(A,B)),c(A,B)) * (1.0/360)
        T4 = c(A,c(c(c(A,B),B),B)) * (1.0/180)
        T5 = c(c(A,B),c(c(A,B),B)) * (1.0/120)
        T6 = c(c(c(c(A,B),B),B),B) * (-1.0/720)

        #for i,t in enumerate([T0,T1,T2,T3,T4,T5,T6]):
        #    print(f"\nT{i}:")
        #    print(t)

        res += T0 + T1 + T2 + T3 + T4 + T5 + T6

    #print("\n5th:")
    #print(res)
    return res

def compute_bch_terms_rec(
        sum: list[any],
        terms: list[any],
        commutator_func: any,
        max_term: int,
        bch_terms: dict,
) -> list[list[any]]:
    ...
