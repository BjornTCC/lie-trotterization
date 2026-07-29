import math
import warnings
import numpy as np

from scipy.optimize import minimize

from src.resource_estimates.gate_costs.protocol import ResourceGate, gate_labels, depth_labels
from src.resource_estimates.rotation_synthesis import synthesize_resource_dict_rotations

def adaptive_phase_estimation_resources(
        target_error: float,
        gates_per_trotter_step: dict[ResourceGate: int],
        unitary_gates: dict[ResourceGate: int] = {},
        num_simulation_steps: int | None = None,
        simulation_steps_and_time_from_error: callable = None,
        unitary_error_coefficients: dict[int: float] = None,
        x: float | None = None,
        synthesize_rotation_with: str = "RUS"
) -> dict[str: int]:
    """
    :param target_error: Desired error in the computed energy
    :param gates_per_trotter_step: ResourceGate counts for each trotter step
    :param unitary_gates: Outer unitary gates in the simulation circuit
    :param num_simulation_steps: Number of calls to the simulation circuits
    :param simulation_steps_and_time_from_error: A function that takes target simulation error and returns a tuple of (evolution time (float), simulation steps (int))
    :param unitary_error_coefficients: Dictionary of error coefficients in the error of the simulation circuit, keys being the power of time, values being coeffiecients.
    :param x: Proportion of error allocated to rotation synthesis. If None, then gate counts are returned in Clifford + T + Tdag + r
    :param synthesize_rotation_with: Strategy to synthesize rotation gates
    :return: dict[str: int] gate counts in {clifford + T + Tdg} or {clifford + T + Tdg + Rx + Ry + Rz} if x is None
    """
    measurement_ancillae = 1
    bit_precision_constant = 0.76*math.pi
    if x is None:
        simulation_error = target_error
    else:
        simulation_error = target_error*(1-x)
        rotation_error = target_error*x

    _, Npe = _simulation_time_and_number_of_calls(bit_precision_constant, simulation_error, num_simulation_steps, simulation_steps_and_time_from_error, unitary_error_coefficients)
    res = {
        "ancilla_qubits": measurement_ancillae
    }

    unitary_applications = math.ceil(np.log2(Npe))

    for nlayers, gate_bundle in zip([unitary_applications, Npe], [unitary_gates, gates_per_trotter_step]):
        for gate, count in gate_bundle.items():
            try:
                cgate = gate.symmetric_controlled()
            except NotImplementedError:
                cgate = gate
            for g in gate_labels + depth_labels:
                if g in res.keys():
                    res[g] += nlayers * count * getattr(cgate, g)
                else:
                    res[g] = nlayers * count * getattr(cgate, g)

    if x is None:
        return {x: y for x, y in res.items() if y != 0}

    return synthesize_resource_dict_rotations(res, rotation_error, synthesize_rotation_with)

def quantum_phase_estimation_resources(
        target_error: float,
        gates_per_trotter_step: dict[ResourceGate: int],
        unitary_gates: dict[ResourceGate: int] = {},
        num_simulation_steps: int | None = None,
        simulation_steps_and_time_from_error: callable = None,
        unitary_error_coefficients: dict[int: float] = None,
        x: float | None = None,
        synthesize_rotation_with: str = "RUS"
) -> dict[str: int]:
    """
        :param target_error: Desired error in the computed energy
        :param gates_per_trotter_step: ResourceGate counts for each trotter step
        :param unitary_gates: Outer unitary gates in the simulation circuit
        :param num_simulation_steps: Number of calls to the simulation circuits
        :param simulation_steps_and_time_from_error: A function that takes target simulation error and returns a tuple of (evolution time (float), simulation steps (int))
        :param unitary_error_coefficients: Dictionary of error coefficients in the error of the simulation circuit, keys being the power of time, values being coeffiecients.
        :param x: Proportion of error allocated to rotation synthesis. If None, then gate counts are returned in Clifford + T + Tdag + r
        :param synthesize_rotation_with: Strategy to synthesize rotation gates
        :return: dict[str: int] gate counts in {clifford + T + Tdg} or {clifford + T + Tdg + Rx + Ry + Rz} if x is None
        """
    bit_precision_constant = 0.5*math.pi
    if x is None:
        simulation_error = target_error
    else:
        simulation_error = target_error*(1-x)
        rotation_error = target_error*x

    t, Npe = _simulation_time_and_number_of_calls(bit_precision_constant, simulation_error, num_simulation_steps, simulation_steps_and_time_from_error, unitary_error_coefficients)
    measurement_ancillae = math.ceil(np.log2(math.pi / (np.sqrt(2)*simulation_error*t)))
    res = {
        "ancilla_qubits": measurement_ancillae,
        "h": measurement_ancillae - 1,
        "rz": measurement_ancillae - 1,
    }

    unitary_applications = 1

    for nlayers, gate_bundle in zip([unitary_applications, Npe], [unitary_gates, gates_per_trotter_step]):
        for gate, count in gate_bundle.items():
            try:
                cgate = gate.symmetric_controlled()
            except NotImplementedError:
                cgate = gate
            for g in gate_labels + depth_labels:
                if g in res.keys():
                    res[g] += nlayers * count * getattr(cgate, g)
                else:
                    res[g] = nlayers * count * getattr(cgate, g)

    if x is None:
        return {x: y for x, y in res.items() if y != 0}

    return synthesize_resource_dict_rotations(res, rotation_error, synthesize_rotation_with)

def _simulation_time_and_number_of_calls(
    bit_precision_constant: float,
    target_error: float,
    num_simulation_steps: int | None = None,
    simulation_steps_and_time_from_error: callable = None,
    unitary_error_coefficients: dict[int: float] = None,
) -> tuple[float, int]:
    if simulation_steps_and_time_from_error is not None:
        return simulation_steps_and_time_from_error(target_error)
    if num_simulation_steps:
        return bit_precision_constant/(target_error * num_simulation_steps), num_simulation_steps
    if unitary_error_coefficients:
        def f(t) -> float:
            res = 0
            for power, coeff in unitary_error_coefficients.items():
                res += coeff * t ** power
            return bit_precision_constant/(t*target_error - res)

        min_power = min(unitary_error_coefficients)
        t0 = (target_error / (min_power * unitary_error_coefficients[min_power]))**(1/(min_power - 1))

        tm = (target_error / unitary_error_coefficients[min_power]) ** (1 / (min_power - 1))

        bnds = [(0.00001, tm)]  # These bounds ensure Npe > 0
        res = minimize(f, t0, bounds = bnds)
        if not res.success:
            warnings.warn(
                f"minimize failed to converge with msg:\n{res.message}"
            )

        return res.x[0], math.ceil(res.fun)