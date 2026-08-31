import math

from src.hubbard_models.hexagonal_hubbard import (
    compute_evolution_time_and_number_of_simulation_circuits_for_qpe_kwargs,
    compute_number_of_trotter_steps_kwargs
)

from src.resource_estimates.gate_costs.protocol import ResourceGate
from src.resource_estimates.gate_costs.free_fermionic import FreeFermionicS1Tile, FreeFermionicS2Tile, FreeFermionicS3Tile
from src.resource_estimates.gate_costs.interactions import ShiftedOccupationPair
from src.resource_estimates.gate_costs.hubbard_commutators import (
    SpinSymmetricMixedControlledKappa,
    SpinSymmetricMixedControlledHappa
)
from src.resource_estimates.gate_costs.hamming_weight_phasing import HWPGate
from src.resource_estimates.gate_costs.parallelized_gate import ParallelizedResourceGate

from src.resource_estimates.quantum_phase_estimation import (
    quantum_phase_estimation_resources,
    adaptive_phase_estimation_resources,
)
from src.resource_estimates.hamiltonian_simulation import (
    hamiltonian_simulation_cost
)

import warnings

def hubbard_model_phase_estimation_resources(
        U: float,
        tau: float,
        Lx: int,
        Ly: int,
        target_error: float,
        phase_estimation_algorithm: str = "adaptive",
        time_evolution_algorithm: str = "tile",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        hwp: bool = False
) -> dict[str: int]:
    if time_evolution_algorithm == "qubitization":
        return get_qubitization_gate_counts(U, tau, Lx, Ly, target_error)
    extra_kwargs =  compute_evolution_time_and_number_of_simulation_circuits_for_qpe_kwargs(
        target_error if x is None else target_error * (1 - x),
        U,
        tau,
        Lx,
        Ly,
        time_evolution_algorithm,
        unitary_decomp=False
    )
    trotter_step_gates, unitary_gates = get_qpe_gate_set(time_evolution_algorithm, Lx, Ly, hwp)

    match phase_estimation_algorithm:
        case "adaptive":
            return adaptive_phase_estimation_resources(
                target_error,
                trotter_step_gates,
                unitary_gates,
                **extra_kwargs,
                x=x,
                synthesize_rotation_with=synthesize_rotation_with
            )

        case "QPE":
            return quantum_phase_estimation_resources(
                target_error,
                trotter_step_gates,
                unitary_gates,
                **extra_kwargs,
                x=x,
                synthesize_rotation_with=synthesize_rotation_with
            )
        case _:
            raise ValueError(f"phase_estimation_algorithm \"{phase_estimation_algorithm}\" not recognized/implemented")

def hubbard_model_time_evolution_resources(
        U: float,
        tau: float,
        Lx: int,
        Ly: int,
        t: float,
        target_error: float,
        algorithm: str = "tile",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        hwp: bool = False
) -> dict[str: int]:
    extra_kwargs = compute_number_of_trotter_steps_kwargs(
        t,
        target_error if x is None else target_error * (1 - x),
        U,
        tau,
        Lx,
        Ly,
        algorithm
    )
    trotter_step_gates, unitary_gates = get_time_evolution_gate_set(algorithm, Lx, Ly, hwp)

    return hamiltonian_simulation_cost(
        t,
        target_error,
        trotter_step_gates,
        unitary_gates,
        **extra_kwargs,
        x=x,
        synthesize_rotation_with=synthesize_rotation_with
    )

def get_qpe_gate_set(type: str, Lx: int, Ly: int, hwp: bool) -> tuple[dict[ResourceGate: int], ...]:
    N = 2*Lx*Ly
    match type:
        case "tile":
            if hwp:
                trotter_step_gates = {
                    HWPGate(
                        {FreeFermionicS2Tile(): N // 2},
                        N
                    ): 5,
                    HWPGate({ShiftedOccupationPair(): N}, N): 1
                }
                unitary_gates = {
                    HWPGate({ShiftedOccupationPair(): N}, N): 1
                }
            else:
                trotter_step_gates = {
                    ParallelizedResourceGate(
                        {FreeFermionicS2Tile(): N // 2}
                    ): 5,
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 1
                }
                unitary_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 1
                }
        case "tile 4th":
            if hwp:
                trotter_step_gates = {
                    HWPGate(
                        {FreeFermionicS2Tile(): N // 2},
                        N
                    ): 5*5,
                    HWPGate({ShiftedOccupationPair(): N}, N): 5
                }
                unitary_gates = {
                    HWPGate({ShiftedOccupationPair(): N}, N): 1
                }
            else:
                trotter_step_gates = {
                    ParallelizedResourceGate(
                        {FreeFermionicS2Tile(): N // 2}
                    ): 5*5,
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 5
                }
                unitary_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 1
                }
        case "augmented tile":
            if hwp:
                trotter_step_gates = {
                    HWPGate(
                        {FreeFermionicS1Tile(): N},
                        2*N
                    ): 15,
                    HWPGate({ShiftedOccupationPair(): N}, N): 2,
                    HWPGate({SpinSymmetricMixedControlledHappa(): N // 2}, N): 3,
                }
                unitary_gates = {
                    HWPGate(
                        {FreeFermionicS1Tile(): N},
                        2*N
                    ): 15,
                    HWPGate({SpinSymmetricMixedControlledKappa(): N // 2}, N): 3,
                }
            else:
                trotter_step_gates = {
                    ParallelizedResourceGate(
                        {FreeFermionicS1Tile(): N}
                    ): 15,
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 2,
                    ParallelizedResourceGate({SpinSymmetricMixedControlledHappa(): N}): 3,
                }
                unitary_gates = {
                    ParallelizedResourceGate(
                        {FreeFermionicS1Tile(): N}
                    ): 15,
                    ParallelizedResourceGate({SpinSymmetricMixedControlledKappa(): N}): 3,
                }
        case _:
            raise ValueError(f"Unrecognized time evolution algorithm: {type}")

    return trotter_step_gates, unitary_gates

def get_time_evolution_gate_set(type: str, Lx: int, Ly: int, hwp: bool) -> tuple[dict[ResourceGate: int], ...]:
    res1, res2 = get_qpe_gate_set(type, Lx, Ly, hwp)
    if type == "augmented tile":
        for x in res2.keys():
            res2[x] *= 2
    return res1, res2

def get_qubitization_gate_counts(U: float, tau: float, Lx: int, Ly: int, target_error: float) -> dict[str: int]:
    warnings.warn(
        "For hexagonal hubbard model using qubitization, only T gate counts are available"
    )
    Cs = 40*Lx*Ly - 4
    Cp = 0 # This contains only polylog factors that we neglect, much like square model.
    Cr = 0
    lambd = (3*tau + U/4)*2*Lx*Ly
    return math.ceil(lambd*math.pi/(2*target_error)) * (Cs + 2*Cp + Cr)
