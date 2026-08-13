import math
import numpy as np
from scipy.special import binom

def fourth_order_suzuki_trotter_commutator_bound(num_terms: int) -> dict[tuple[int, ...]: float]:
    # A wrapper for the code below special to my use case
    rule = SplittingMethod.suzuki(num_terms, 2)

    s = (rule.num_layers + 1) // 2

    res = commutator_bound(rule, s)
    return {tuple([y + 1 for y in x.commidx[::-1]]): x.weight.item() for x in res}

"""
This code is entirely copied from:
https://github.com/qc-tum/fermi_hubbard_commutators

Article citation:
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.108.195105
"""

def multinomial(params):
    """
    Multinomial coefficient.
    """
    if len(params) == 1:
        return 1
    return binom(sum(params), params[-1]) * multinomial(params[:-1])


def integer_sum_tuples(s: int, nbins: int):
    """
    Generate lexicographically sorted non-negative integer lists
    such that the sum of the integers equals 's'.
    """
    if nbins <= 0:
        raise ValueError(f"'nbins' must be at least 1, received {nbins}")
    if s < 0:
        raise ValueError(f"'s' cannot be negative, received {s}")
    c = nbins * [0]
    c[-1] = s
    yield tuple(c)
    done = False
    while not done:
        done = True
        for i in reversed(range(1, nbins)):
            if c[i] > 0:
                c[i] -= 1
                # swap c[-1] <-> c[i]
                c[-1], c[i] = c[i], c[-1]
                c[i - 1] += 1
                yield tuple(c)
                done = False
                break

def _construct_suzuki_indices_coeffs(nterms: int, k: int):
    """
    Recursively construct the Suzuki product rule indices and coefficients for order `2 k`.
    """
    if k <= 0:
        raise ValueError(f"`k` must be a positive integer, received {k}")
    if k == 1:
        indices = list(range(nterms)) + list(reversed(range(nterms)))
        coeffs = (2*nterms) * (0.5,)
    else:
        uk = 1./(4 - 4**(1./(2*k-1)))
        ik1, ck1 = _construct_suzuki_indices_coeffs(nterms, k - 1)
        ck1_uk = [uk*c for c in ck1]
        ck1_14uk = [(1 - 4*uk)*c for c in ck1]
        indices = ik1 + ik1 + ik1 + ik1 + ik1
        coeffs = ck1_uk + ck1_uk + ck1_14uk + ck1_uk + ck1_uk
    return merge_layers(indices, coeffs)


def merge_layers(indices, coeffs):
    """
    Merge neighboring layers with the same index.
    """
    assert len(coeffs) == len(indices)
    mindices = [indices[0]]
    mcoeffs  = [coeffs[0]]
    for i, c in zip(indices[1:], coeffs[1:]):
        if mindices[-1] == i:
            mcoeffs[-1] += c
        else:
            mindices.append(i)
            mcoeffs.append(c)
    return mindices, mcoeffs

class SplittingMethod:
    """
    Splitting method described by the number of (Hamiltonian) terms
    (typically two, as for even-odd splitting), indices into these terms,
    and corresponding coefficients (time sub-step coefficients).
    """
    def __init__(self, nterms: int, indices, coeffs, order: int):
        # consistency check
        for i in indices:
            if i < 0 or i >= nterms:
                raise ValueError(f"index {i} out of range, must be between 0 and {nterms}")
        if len(coeffs) != len(indices):
            raise ValueError("length of coefficient list and index list must agree")
        weights = np.zeros(nterms)
        for i, c in zip(indices, coeffs):
            weights[i] += c
        if not np.allclose(weights, np.ones(nterms)):
            raise ValueError("weights for each term do not sum to 1")
        self.nterms  = nterms
        self.indices = list(indices)
        self.coeffs  = list(coeffs)
        self.order   = order

    @classmethod
    def trotter(cls, nterms: int):
        """
        Construct the first-order Lie-Trotter splitting formula.
        """
        return cls(nterms, range(nterms), nterms * [1], 1)

    @classmethod
    def suzuki(cls, nterms: int, k: int):
        """
        Construct the Suzuki product rule for order `2 k`.
        """
        indices, coeffs = _construct_suzuki_indices_coeffs(nterms, k)
        return cls(nterms, indices, coeffs, 2*k)

    @property
    def num_terms(self):
        """
        Number of (Hamiltonian) terms.
        """
        return self.nterms

    @property
    def num_layers(self):
        """
        Number of layers (substeps).
        """
        return len(self.coeffs)

    def __str__(self):
        """
        String representation of the product rule.
        """
        return f"Splitting method of order {self.order} for {self.nterms} terms using {self.num_layers} layers,\n  indices: {self.indices}\n  coeffs:  {self.coeffs}"

class WeightedNestedCommutator:
    """
    Symbolic weighted nested commutator.

    A commutation index `i` is interpreted as
    [A_{i[-1]}, ...[A_{i[2]}, [A_{i[1]}, A_{i[0]}]]...]
    """
    def __init__(self, commidx, weight):
        self.commidx = commidx
        self.weight = weight

    def __str__(self):
        """
        String representation of the weighted nested commutator.
        """
        s = f"H_{self.commidx[0]}"
        for i in self.commidx[1:]:
            s = f"[H_{i}, {s}]"
        s = str(self.weight) + " * " + s
        return s


def commutator_bound(splitting: SplittingMethod, s: int):
    """
    Coefficients for commutator bounds on a splitting method (product rule).
    """
    weights = np.zeros((splitting.order + 1) * (splitting.num_terms,))
    for j in range(1, splitting.num_layers):
        bcoeff = np.zeros(splitting.num_terms)
        for i, c in zip(splitting.indices[:j], splitting.coeffs[:j]):
            bcoeff[i] += c
        for q in integer_sum_tuples(splitting.order, s-j if j < s else j-s+1):
            if q[0] == 0:
                continue
            mq = multinomial(q)
            for k in range(splitting.num_terms):
                if bcoeff[k] == 0:
                    continue
                commidx = (k,)
                w = bcoeff[k]
                for i in range(len(q)):
                    l = j + i if j < s else j - i
                    commidx += q[i] * (splitting.indices[l],)
                    w *= abs(splitting.coeffs[l])**q[i]
                assert len(commidx) == splitting.order + 1
                if commidx[0] == commidx[1]:
                    # [A, A] = 0
                    continue
                if commidx[0] > commidx[1]:
                    # [A, B] = -[B, A], and absolute values agree
                    commidx = (commidx[1], commidx[0]) + commidx[2:]
                weights[commidx] += mq * w
    weights /= math.factorial(splitting.order + 1)
    # assemble return value
    res = []
    for idx, w in np.ndenumerate(weights):
        if w != 0:
            res.append(WeightedNestedCommutator(idx, w))
    return res
