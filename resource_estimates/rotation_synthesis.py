"""
Functions to calculate the T gate count for rotation synthesis
"""

import math
import mpmath
import numpy as np

from pygridsynth.gridsynth import gridsynth_gates

def gridsynth_rotation_cost(theta: float | list[float], eps: float,  max_samples: int = 100) -> int:
    if isinstance(theta, float):
        gates = gridsynth_gates(theta=mpmath.mpf(theta), epsilon=mpmath.mpf(eps))
        return gates.count("T")
    elif max_samples == -1 or len(theta) < max_samples:
        res = 0
        for t in theta:
            res += gridsynth_rotation_cost(t, eps / len(theta))
        return res
    else:
        res = 0
        angles = np.random.choice(theta, size = max_samples, replace = False)
        for t in angles:
            res += gridsynth_rotation_cost(t, eps / len(theta))
        return math.ceil(res * len(theta) / max_samples)

def rus_rotation_cost(theta: float | list[float], eps: float) -> int:
    if isinstance(theta, float):
        nr = 1
    else:
        nr = len(theta)
    return math.ceil(nr * (1.149 * np.log2(nr / eps) + 9.2))