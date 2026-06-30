import math
import warnings

from scipy.optimize import fsolve

from src.resource_estimates.gate_costs.protocol import ResourceGate, gate_labels
from src.resource_estimates.rotation_synthesis import synthesize_resource_dict_rotations

def hamiltonian_simulation_cost(
        time: float,
        target_error: float,
        gates_per_trotter_step: dict[ResourceGate: int],
        outer_gates: dict[ResourceGate: int] = {},
        num_trotter_steps: int | None = None,
        trotter_steps_from_time_and_error: callable = None,
        error_coefficients: dict[int: float] = None,
        x: float | None = None,
        synthesize_rotation_with: str = "RUS"
) -> dict[str: int]:

    """
    Function that computes the gate counts for hamiltonian simulation.
    :param time: Simulation time
    :param target_error: Desired error
    :param gates_per_trotter_step: A dictionary of ResourceGate that make up the trotter
    :param outer_gates: Any gates that are executed once at the beginning and or end. Specify the total number of gates that are executed
    :param num_trotter_steps: The number of trotter steps to execute
    :param trotter_steps_from_time_and_error: A function that takes simulation time and target error and returns an integer representing the number of trotter steps
    :param error_coefficients: coefficents in a polynomial that estimates the trotter error
    :param x: Proportion of the error to allocate to rotation synthesis, must be between 0 and 1, or None
    :param synthesize_rotation_with: method to synthesize rotations. Options: RUS, gridsynth
    :return: dict[str: int] gate counts in {clifford + T + Tdg} or {clifford + T + Tdg + Rx + Ry + Rz} if x is None
    """

    trotter_error = (1 - x) * target_error if x is not None else target_error

    _, n = _compute_time_step_and_num_trotter_steps(time, trotter_error, num_trotter_steps, trotter_steps_from_time_and_error, error_coefficients)
    res = {}

    for nlayers, gate_bundle in zip([1, n], [outer_gates, gates_per_trotter_step]):
        for gate, count in gate_bundle.items():
            for g in gate_labels:
                if g in res.keys():
                    res[g] += nlayers*count*getattr(gate, g)
                else:
                    res[g] = nlayers*count*getattr(gate, g)

    if x is None:
        return {x: y for x, y in res.items() if y!= 0}

    rotation_error = x * target_error
    return synthesize_resource_dict_rotations(res, rotation_error, synthesize_rotation_with)

def _compute_time_step_and_num_trotter_steps(
        time: float,
        target_error: float,
        num_trotter_steps: int | None = None,
        trotter_steps_from_time_and_error: callable = None,
        error_coefficients: dict[int: float] = None,
) -> tuple[float, int]:
    if num_trotter_steps is not None:
        n = num_trotter_steps
    elif trotter_steps_from_time_and_error is not None:
        n = trotter_steps_from_time_and_error(time, target_error)
    elif error_coefficients is not None:
        def f(n) -> float:
            res = 0
            for power, coeff in error_coefficients.items():
                res += coeff * time ** power / (n ** (power - 1))
            return res - target_error

        min_power = min(error_coefficients)
        n0 = (error_coefficients[min_power] / target_error) ** (1 / (min_power - 1)) * time ** (
                    1 + 1 / (min_power - 1))

        root, _, ier, mesg = fsolve(f, n0, full_output=True)
        if not ier:
            warnings.warn(
                f"fsolve failed to converge with msg:\n{mesg}"
            )

        n = math.ceil(root[0])
    else:
        raise ValueError(
            f"One of \'num_trotter_steps\', \'trotter_steps_from_time_and_error\' or \'error_coefficients\' must be specified.")
    return time / n, n
