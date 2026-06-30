"""
Functions to calculate the T gate count for rotation synthesis
"""

from copy import deepcopy
import math
import mpmath
import numpy as np
import warnings

from pygridsynth.gridsynth import gridsynth_gates

def gridsynth_rotation_cost(theta: float | list[float], eps: float,  max_samples: int = 100) -> tuple[int,int]:
    if isinstance(theta, float):
        gates = gridsynth_gates(theta=mpmath.mpf(theta), epsilon=mpmath.mpf(eps))
        return gates.count("T")
    elif max_samples == -1 or len(theta) < max_samples:
        res = 0
        for t in theta:
            res += gridsynth_rotation_cost(t, eps / len(theta))
        return res, 0
    else:
        res = 0
        angles = np.random.choice(theta, size = max_samples, replace = False)
        for t in angles:
            res += gridsynth_rotation_cost(t, eps / len(theta))
        return math.ceil(res * len(theta) / max_samples), 0

def rus_rotation_cost(theta: float | list[float], eps: float) -> tuple[int,int]:
    if isinstance(theta, float):
        nr = 1
    else:
        nr = len(theta)
    return math.ceil(nr * (1.149 * np.log2(nr / eps) + 9.2)), math.ceil(2 * nr*np.log(nr/eps) / np.log(5))

def synthesize_resource_dict_rotations(
        resources: dict[str: int],
        target_error: float,
        synthesis_strategy: float
) -> dict[str: int]:
    res = deepcopy(resources)
    nr = sum([res["rz"], res["rx"], res["ry"]])
    warnings.warn(
        "When synthesizing rotation into clifford plus T, this function does not compute the single-qubit clifford count due to arbitrary rotations.")
    match synthesis_strategy:
        case "gridsynth":
            t_cost, cx_cost = gridsynth_rotation_cost(np.random.uniform(size=nr), target_error)
        case "RUS":
            t_cost, cx_cost = rus_rotation_cost([0.0], target_error/nr)
            t_cost *= nr
            cx_cost *= nr
        case _:
            raise ValueError(f"Rotation synthesis method: {synthesis_strategy} not recognized/implemented.")

    del res["rx"]
    del res["ry"]
    del res["rz"]

    if "t" in res.keys():
        res["t"] += t_cost
    else:
        res["t"] = t_cost
    if "cx" in res.keys():
        res["cx"] += cx_cost
    else:
        res["cx"] = cx_cost
    return {x: y for x, y in res.items() if y != 0}