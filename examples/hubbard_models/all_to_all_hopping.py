import numpy as np

from ._primitives import site_hopping

from openfermion import FermionOperator

def all_to_all_random(size: int) -> FermionOperator:

    terms = []

    for i in range(size):
        for j in range(i+1, size):
            terms.append(
                np.random.uniform(-1,1) * site_hopping(i,j)
            )
    return sum(terms)