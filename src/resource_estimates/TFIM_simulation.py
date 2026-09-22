import networkx as nx

from src.TFIM.split_operator_error_coefficients import (
    second_order_coefficient,
    fourth_order_coefficient,
    augmented_coefficients
)
from src.TFIM.standard_interaction_models import (
    nearest_neighbor_from_graph,
    hexagonal_power_law,
    square_power_law,
    cubic_power_law,
)
from src.TFIM.standard_interaction_gates import (
    gates_from_graph_2nd,
    gates_from_graph_4th,
    gates_from_graph_augmented,
    gates_from_positions_and_interaction_augmented,
    gates_from_positions_and_interaction_2nd_order,
    gates_from_positions_and_interaction_4th_order
)
from src.resource_estimates.quantum_phase_estimation import (
    quantum_phase_estimation_resources,
    adaptive_phase_estimation_resources,
)
from src.resource_estimates.hamiltonian_simulation import hamiltonian_simulation_cost

def TFIM_phase_estimation_resources(
        target_error: float,
        lattice: str | nx.Graph,
        J: float,
        U: float,
        simulation_type: str,
        phase_estimation_algorithm: str = "qpe",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        alpha: float = None,
        Ls: list[int] = None,
        max_dist: float = None,
        tol: float = 1e-14,
) -> dict[str: int]:
    if isinstance(lattice, nx.Graph):
        interaction, positions = nearest_neighbor_from_graph(lattice, J)
    elif lattice == "honeycomb":
        if not len(Ls) == 2 and (alpha is not None):
            raise ValueError(f"When lattice = honeycomb, please specify Ls = [Lx, Ly] and alpha")
        interaction, positions = hexagonal_power_law(*Ls, J, alpha, max_dist)

    elif lattice == "square":
        if not len(Ls) == 2 and (alpha is not None):
            raise ValueError(f"When lattice = square, please specify Ls = [Lx, Ly] and alpha")
        interaction, positions = square_power_law(*Ls, J, alpha, max_dist)

    elif lattice == "cubic":
        if not len(Ls) == 3 and (alpha is not None):
            raise ValueError(f"When lattice = cubic, please specify Ls = [Lx, Ly, Lz] and alpha")
        interaction, positions = cubic_power_law(*Ls, J, alpha, max_dist)
    else:
        raise NotImplementedError(f"lattice: \"{lattice}\" not recognized / implemented")


    match simulation_type:
        case "2nd order":
            error_coeffs = {3: second_order_coefficient(positions, interaction, U)}
            if isinstance(lattice, nx.Graph):
                gates = gates_from_graph_2nd(G)
            else:
                gates = gates_from_positions_and_interaction_2nd_order(positions, interaction, tol)
        case "4th order":
            error_coeffs = {5: fourth_order_coefficient(positions, interaction, U)}
            if isinstance(lattice, nx.Graph):
                gates = gates_from_graph_4th(G)
            else:
                gates = gates_from_positions_and_interaction_4th_order(positions, interaction, tol)
        case "augmented":
            error_coeffs = augmented_coefficients(positions, interaction, U, unitary_decomp=False)
            if isinstance(lattice, nx.Graph):
                gates = gates_from_graph_augmented(G)
            else:
                gates = gates_from_positions_and_interaction_augmented(positions, interaction, tol)
        case _:
            raise ValueError(f"Simulation type \"{simulation_type}\" not recognized/implemented")

    match phase_estimation_algorithm:
        case "adaptive":
            return adaptive_phase_estimation_resources(
                target_error,
                gates_per_trotter_step = gates[0],
                unitary_gates = gates[1],
                unitary_error_coefficients = error_coeffs,
                x = x,
                synthesize_rotation_with=synthesize_rotation_with
            )
        case "qpe":
            return quantum_phase_estimation_resources(
                target_error,
                gates_per_trotter_step = gates[0],
                unitary_gates = gates[1],
                unitary_error_coefficients = error_coeffs,
                x = x,
                synthesize_rotation_with=synthesize_rotation_with
            )

def TFIM_time_evolution_resources(
        time: float,
        target_error: float,
        lattice: str | nx.Graph,
        J: float,
        U: float,
        simulation_type: str,
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        alpha: float = None,
        Ls: list[int] = None,
        max_dist: float = None,
        tol: float = 1e-14,
) -> dict[str: int]:
    if isinstance(lattice, nx.Graph):
        interaction, positions = nearest_neighbor_from_graph(lattice, J)
    elif lattice == "honeycomb":
        if not len(Ls) == 2 and (alpha is not None):
            raise ValueError(f"When lattice = honeycomb, please specify Ls = [Lx, Ly] and alpha")
        interaction, positions = hexagonal_power_law(*Ls, J, alpha, max_dist)

    elif lattice == "square":
        if not len(Ls) == 2 and (alpha is not None):
            raise ValueError(f"When lattice = square, please specify Ls = [Lx, Ly] and alpha")
        interaction, positions = square_power_law(*Ls, J, alpha, max_dist)

    elif lattice == "cubic":
        if not len(Ls) == 3 and (alpha is not None):
            raise ValueError(f"When lattice = cubic, please specify Ls = [Lx, Ly, Lz] and alpha")
        interaction, positions = cubic_power_law(*Ls, J, alpha, max_dist)
    else:
        raise NotImplementedError(f"lattice: \"{lattice}\" not recognized / implemented")


    match simulation_type:
        case "2nd order":
            error_coeffs = {3: second_order_coefficient(positions, interaction, U)}
            if isinstance(lattice, nx.Graph):
                gates = gates_from_graph_2nd(G)
            else:
                gates = gates_from_positions_and_interaction_2nd_order(positions, interaction, tol)
        case "4th order":
            error_coeffs = {5: fourth_order_coefficient(positions, interaction, U)}
            if isinstance(lattice, nx.Graph):
                gates = gates_from_graph_4th(G)
            else:
                gates = gates_from_positions_and_interaction_4th_order(positions, interaction, tol)
        case "augmented":
            error_coeffs = augmented_coefficients(positions, interaction, U, unitary_decomp=True)
            if isinstance(lattice, nx.Graph):
                gates = gates_from_graph_augmented(G)
            else:
                gates = gates_from_positions_and_interaction_augmented(positions, interaction, tol)
        case _:
            raise ValueError(f"Simulation type \"{simulation_type}\" not recognized/implemented")

    return hamiltonian_simulation_cost(
        time,
        target_error,
        gates[0],
        gates[1],
        error_coefficients=error_coeffs,
        x=x,
        synthesize_rotation_with=synthesize_rotation_with
    )