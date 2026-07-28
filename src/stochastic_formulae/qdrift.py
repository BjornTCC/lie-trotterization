import math

import numpy as np
from scipy.special import lambertw

import warnings

from qiskit.quantum_info import Kraus, SparsePauliOp

def solve_N_from_eps(eps: float, lambd: float, time: float, exact: bool) -> int:
    if exact:
        fl_res = np.real_if_close(2*lambd*time/lambertw(eps/(lambd*time)))
    else:
        fl_res = 2*(lambd*time)**2/eps
    return math.ceil(fl_res)

def qdrift(
        terms: list[any],
        weights: list[any],
        time: float,
        *args,
        eps: float = None,
        N: float = None,
        exact_error_bound: bool = False,
        seed: int = None
) -> list[any]:
    #print("Qdrift terms:", terms)
    if (eps is None and N is None) or (eps is not None and N is not None):
        raise ValueError("Exactly one of kwargs \"eps\" or \"N\" must be specified")

    if min(weights) < 0:
        raise ValueError("weights must be non-negative")

    lambd = sum(weights)
    probabilities = [h/lambd for h in weights]

    if eps is not None:
        N = solve_N_from_eps(eps, lambd, time, exact_error_bound)
    elif exact_error_bound:
        warnings.warn("Warning: kwarg \"exact_error_bound\" causes no change when \"N\" is specified")

    M = len(terms)
    tau = lambd*time/N
    rng = np.random.default_rng(seed)
    samples = rng.choice(M, size = N, replace = True, p = probabilities)
    return [tau*terms[i] for i in samples]

def qdrift_kraus(
        operator: SparsePauliOp,
        time: float,
        *args,
        eps: float = None,
        N: float = None,
        exact_error_bound: bool = False,
        seed: int = None
) -> Kraus:
    #print("Qdrift terms:", terms)
    if (eps is None and N is None) or (eps is not None and N is not None):
        raise ValueError("Exactly one of kwargs \"eps\" or \"N\" must be specified")

    weights = [abs(x) for x in operator.coeffs]
    signs = [np.sign(x) for x in operator.coeffs]

    lambd = sum(weights)
    sqrt_probabilities = [np.sqrt(h/lambd) for h in weights]

    if eps is not None:
        N = solve_N_from_eps(eps, lambd, time, exact_error_bound)
    elif exact_error_bound:
        warnings.warn("Warning: kwarg \"exact_error_bound\" causes no change when \"N\" is specified")

    ID = SparsePauliOp("I"*operator.num_qubits)

    tau = lambd*time/N
    base_kraus_ops = [
        sp * (ID * np.cos(tau) - 1j*sign*np.sin(tau) * SparsePauliOp(op))
        for sp, sign, op in zip(sqrt_probabilities, signs, operator.paulis)
    ]
    return Kraus(base_kraus_ops)