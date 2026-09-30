import math

from src.PPP_model.split_operator_coefficients import (
    second_order_coefficient,
    augmented_coefficients,
    dominant_augmented_coefficients,
    dominant_4th_order_coefficients
)

from src.resource_estimates.gate_costs.protocol import ResourceGate
from src.resource_estimates.gate_costs.free_fermionic import FreeFermionicS1Tile, FreeFermionicS2Tile, FreeFermionicS3Tile
from src.resource_estimates.gate_costs.interactions import ShiftedOccupationPair, AncillaOccupationPair, AncillaFreeOccupationPair
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

def PPP_hamiltonian_simulation_resources(
        U: float,
        v: float,
        Lx: int,
        Ly: int,
        t: float,
        target_error: float,
        algorithm: str = "2nd order",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
) -> dict[str: int]:

    d = 3
    N = 2*Lx*Ly
    match algorithm:
        case "2nd order":
            error_coeffs = {3: second_order_coefficient(U, v, d, N)}
        case "4th order":
            error_coeffs = {5: dominant_4th_order_coefficients(U, v, d, N)}
        case "augmented":
            error_coeffs = augmented_coefficients(U, v, d, N)
        case "dominant augmented":
            error_coeffs = dominant_augmented_coefficients(U, v, d, N)
        case _:
            raise NotImplementedError(f"Algorithm {algorithm} not implemented")

    trotter_step_gates, unitary_gates = PPP_model_gates(2*Lx*Ly, algorithm, double_unitary=True)

    return hamiltonian_simulation_cost(
        t,
        target_error,
        trotter_step_gates,
        unitary_gates,
        error_coefficients=error_coeffs,
        x=x,
        synthesize_rotation_with=synthesize_rotation_with
    )

def PPP_qpe_resources(
        U: float,
        v: float,
        Lx: int,
        Ly: int,
        target_error: float,
        phase_estimation_algorithm: str = "adaptive",
        algorithm: str = "2nd order",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
) -> dict[str: int]:
    d = 3
    N = 2 * Lx * Ly
    match algorithm:
        case "2nd order":
            error_coeffs = {3: second_order_coefficient(U, v, d, N)}
        case "4th order":
            error_coeffs = {5: dominant_4th_order_coefficients(U, v, d, N)}
        case "augmented":
            error_coeffs = augmented_coefficients(U, v, d, N)
        case "dominant augmented":
            error_coeffs = dominant_augmented_coefficients(U, v, d, N)
        case _:
            raise NotImplementedError(f"Algorithm {algorithm} not implemented")

    trotter_step_gates, unitary_gates = PPP_model_gates(2 * Lx * Ly, algorithm, double_unitary=False)

    match phase_estimation_algorithm:
        case "adaptive":
            return adaptive_phase_estimation_resources(
                target_error,
                trotter_step_gates,
                unitary_gates,
                unitary_error_coefficients=error_coeffs,
                x=x,
                synthesize_rotation_with=synthesize_rotation_with
            )
        case "qpe":
            return quantum_phase_estimation_resources(
                target_error,
                trotter_step_gates,
                unitary_gates,
                unitary_error_coefficients=error_coeffs,
                x=x,
                synthesize_rotation_with=synthesize_rotation_with
            )

def PPP_model_gates(
        N: int,
        type: str,
        double_unitary: bool
) -> dict[ResourceGate: int]:
    d = 3
    match type:
        case "2nd order":
            trotter_step_gates = {
                FreeFermionicS2Tile(): 5*N // 2,
                AncillaOccupationPair(): 2*N*(N-1) + N
            }
            unitary_gates = {
                FreeFermionicS2Tile(): 5*N // 2,
            }
        case "4th order":
            trotter_step_gates = {
                FreeFermionicS2Tile(): 25*N // 2,
                AncillaOccupationPair(): 10*N*(N-1) + 5*N
            }
            unitary_gates = {
                FreeFermionicS2Tile(): 5*N // 2,
            }
        case "augmented" | "dominant augmented":
            trotter_step_gates = {
                FreeFermionicS1Tile(): 9*N,
                FreeFermionicS2Tile(): 3*N,
                AncillaOccupationPair(): 2*N*(N-1) + N,
                SpinSymmetricMixedControlledKappa(): 2*N**2*d
            }
            unitary_gates = {
                FreeFermionicS1Tile(): 9*N,
                FreeFermionicS2Tile(): 3*N,
                SpinSymmetricMixedControlledKappa(): 2*N**2*d
            }
    if double_unitary:
        return trotter_step_gates, {x: 2*y for x, y in unitary_gates.items()}
    return trotter_step_gates, unitary_gates