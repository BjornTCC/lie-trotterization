import math

import networkx as nx
import numpy as np

from scipy.optimize import fsolve, minimize_scalar

from src.hubbard_models._free_fermionic_computations import multiply_graph_by
from src.hubbard_models.free_fermionic_errors import error_between_exp_of_free_fermionic, second_order_error_coefficient_of_free_fermionic_operator
from src.hubbard_models.split_operator_error_coefficients import (
    second_order_split_operator_error_coefficient,
    fourth_order_suzuki_trotter_split_operator_error_coefficient,
    fourth_order_augmented_split_operator_error_coefficients,
)
from src.hubbard_models.s1_tiling_strategy import (
    second_order_s1_tiling,
    augmented_s1_tiling,
    _color_edge_decomp
)

from src.resource_estimates.gate_costs.protocol import ResourceGate
from src.resource_estimates.gate_costs.free_fermionic import FreeFermionicS1Tile
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

def hubbard_model_phase_estimation_resources(
        U: float,
        tau: float,
        hopping_graph: nx.Graph,
        target_error: float,
        phase_estimation_algorithm: str = "adaptive",
        time_evolution_algorithm: str = "trotter",
        decomp: list[nx.Graph] = None,
        decomp_strategy: str = "DSATUR",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        hwp: bool = False
) -> dict[str: int]:
    match time_evolution_algorithm:
        case "trotter":
            return _hubbard_model_phase_estimation_resources_2nd_order(
                U,
                tau,
                hopping_graph,
                target_error,
                phase_estimation_algorithm,
                decomp,
                decomp_strategy,
                x,
                synthesize_rotation_with,
                hwp=hwp
            )

        case "trotter_4th":
            return _hubbard_model_phase_estimation_resources_4th_order_suzuki_trotter(
                U,
                tau,
                hopping_graph,
                target_error,
                phase_estimation_algorithm,
                decomp,
                decomp_strategy,
                x,
                synthesize_rotation_with,
                hwp=hwp
            )
        case "augmented":
            return _hubbard_model_phase_estimation_resources_augmented(
                U,
                tau,
                hopping_graph,
                target_error,
                phase_estimation_algorithm,
                decomp,
                decomp_strategy,
                x,
                synthesize_rotation_with,
                hwp=hwp
            )

def hubbard_model_time_evolution_resources(
        U: float,
        tau: float,
        hopping_graph: nx.Graph,
        t: float,
        target_error: float,
        algorithm: str = "trotter",
        decomp: list[nx.Graph] = None,
        decomp_strategy: str = "DSATUR",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        hwp: bool = False
) -> dict[str: int]:
    match algorithm:
        case "trotter":
            return _hubbard_model_time_evolution_resources_trotter(
                U,
                tau,
                hopping_graph,
                t,
                target_error,
                decomp,
                decomp_strategy,
                x = x,
                synthesize_rotation_with=synthesize_rotation_with,
                hwp=hwp
            )
        case "trotter_4th":
            return _hubbard_model_time_evolution_resources_trotter_4th_suzuki_trotter(
                U,
                tau,
                hopping_graph,
                t,
                target_error,
                decomp,
                decomp_strategy,
                x = x,
                synthesize_rotation_with=synthesize_rotation_with,
                hwp=hwp
            )
        case "augmented":
            return _hubbard_model_time_evolution_resources_augmented(
                U,
                tau,
                hopping_graph,
                t,
                target_error,
                decomp,
                decomp_strategy,
                x = x,
                synthesize_rotation_with=synthesize_rotation_with,
                hwp=hwp
            )

def _hubbard_model_phase_estimation_resources_2nd_order(
        U: float,
        tau: float,
        hopping_graph: nx.Graph,
        target_error: float,
        phase_estimation_algorithm: str = "adaptive",
        decomp: list[nx.Graph] = None,
        decomp_strategy: str = "DSATUR",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        hwp: bool = False,
) -> dict[str: int]:
    d = max(dict(hopping_graph.degree).values())
    N = len(hopping_graph.nodes)
    if nx.is_weighted(hopping_graph):
        hopping_integral = tau * max([hopping_graph[e[0]][e[1]]["weight"] for e in hopping_graph.edges])
    else:
        hopping_integral = tau

    hopping_decomposition = second_order_s1_tiling(hopping_graph, decomp, strategy=decomp_strategy)

    W = (second_order_split_operator_error_coefficient(U, hopping_integral, N, hopping_graph)
         + second_order_error_coefficient_of_free_fermionic_operator(hopping_decomposition)
         )
    Nedges = len(hopping_graph.edges)

    trotter_step_gates, unitary_gates = _gate_step_qpe(N, None, hopping_decomposition, Nedges, "trotter", hwp)
    match phase_estimation_algorithm:
        case "adaptive":
            return adaptive_phase_estimation_resources(
                target_error,
                trotter_step_gates,
                unitary_gates,
                unitary_error_coefficients={3: W},
                x = x,
                synthesize_rotation_with=synthesize_rotation_with
            )

        case "QPE":
            return quantum_phase_estimation_resources(
                target_error,
                trotter_step_gates,
                unitary_gates,
                unitary_error_coefficients={3: W},
                x = x,
                synthesize_rotation_with=synthesize_rotation_with
            )

def _hubbard_model_phase_estimation_resources_4th_order_suzuki_trotter(
        U: float,
        tau: float,
        hopping_graph: nx.Graph,
        target_error: float,
        phase_estimation_algorithm: str = "adaptive",
        decomp: list[nx.Graph] = None,
        decomp_strategy: str = "DSATUR",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        hwp: bool = False,
) -> dict[str: int]:
    d = max(dict(hopping_graph.degree).values())
    N = len(hopping_graph.nodes)
    if nx.is_weighted(hopping_graph):
        hopping_integral = tau * max([hopping_graph[e[0]][e[1]]["weight"] for e in hopping_graph.edges])
    else:
        hopping_integral = tau

    hopping_decomposition = second_order_s1_tiling(hopping_graph, decomp, strategy=decomp_strategy)

    W = fourth_order_suzuki_trotter_split_operator_error_coefficient(U, hopping_integral, N, d)
    Nedges = len(hopping_graph.edges)

    trotter_step_gates, unitary_gates = _gate_step_qpe(N, None, hopping_decomposition, Nedges, "trotter_4th", hwp)
    match phase_estimation_algorithm:
        case "adaptive":
            return adaptive_phase_estimation_resources(
                target_error,
                trotter_step_gates,
                unitary_gates,
                unitary_error_coefficients={5: W},
                x=x,
                synthesize_rotation_with=synthesize_rotation_with
            )

        case "QPE":
            return quantum_phase_estimation_resources(
                target_error,
                trotter_step_gates,
                unitary_gates,
                unitary_error_coefficients={5: W},
                x=x,
                synthesize_rotation_with=synthesize_rotation_with
            )

def _hubbard_model_phase_estimation_resources_augmented(
        U: float,
        tau: float,
        hopping_graph: nx.Graph,
        target_error: float,
        phase_estimation_algorithm: str = "adaptive",
        decomp: list[nx.Graph] = None,
        decomp_strategy: str = "DSATUR",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        hwp: bool = False,
) -> dict[str: int]:
    d = max(dict(hopping_graph.degree).values())
    N = len(hopping_graph.nodes)
    if nx.is_weighted(hopping_graph):
        hopping_integral = tau * max([hopping_graph[e[0]][e[1]]["weight"] for e in hopping_graph.edges])
    else:
        hopping_integral = tau
    match phase_estimation_algorithm:
        case "adaptive":
            bit_precision_constant = 0.76
        case "QPE":
            bit_precision_constant = 0.5
    _func = lambda e: _augmented_s1_evolution_time_and_steps_qpe(
        e,
        U,
        tau,
        hopping_integral,
        d,
        N,
        hopping_graph,
        decomp,
        decomp_strategy,
        unitary_decomp= False,

    )

    hopping_decomposition = augmented_s1_tiling(hopping_graph, decomp, strategy=decomp_strategy)
    Nedges = len(hopping_graph.edges)

    if decomp is None:
        _decomp = _color_edge_decomp(hopping_graph, strategy=decomp_strategy)
    else:
        _decomp = decomp

    trotter_step_gates, unitary_gates = _gate_step_qpe(N, _decomp, hopping_decomposition, Nedges, "augmented", hwp)
    match phase_estimation_algorithm:
        case "adaptive":
            return adaptive_phase_estimation_resources(
                target_error,
                trotter_step_gates,
                unitary_gates,
                simulation_steps_and_time_from_error=_func,
                x=x,
                synthesize_rotation_with=synthesize_rotation_with
            )

        case "QPE":
            return quantum_phase_estimation_resources(
                target_error,
                trotter_step_gates,
                unitary_gates,
                simulation_steps_and_time_from_error=_func,
                x=x,
                synthesize_rotation_with=synthesize_rotation_with
            )

def _hubbard_model_time_evolution_resources_trotter(
        U: float,
        tau: float,
        hopping_graph: nx.Graph,
        t: float,
        target_error: float,
        decomp: list[nx.Graph] = None,
        decomp_strategy: str = "DSATUR",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        hwp: bool = False,
) -> dict[str: int]:
    d = max(dict(hopping_graph.degree).values())
    N = len(hopping_graph.nodes)
    if nx.is_weighted(hopping_graph):
        hopping_integral = tau * max([hopping_graph[e[0]][e[1]]["weight"] for e in hopping_graph.edges])
    else:
        hopping_integral = tau

    hopping_decomposition = second_order_s1_tiling(hopping_graph, decomp, strategy=decomp_strategy)

    W = (second_order_split_operator_error_coefficient(U, hopping_integral, N, hopping_graph)
         + second_order_error_coefficient_of_free_fermionic_operator(hopping_decomposition)
         )
    Nedges = len(hopping_graph.edges)

    trotter_step_gates, outer_gates = _gate_step_evolution(N, None, hopping_decomposition, Nedges, "trotter", hwp)

    return hamiltonian_simulation_cost(
        t,
        target_error,
        trotter_step_gates,
        outer_gates,
        error_coefficients={3:W},
        x=x,
        synthesize_rotation_with=synthesize_rotation_with
    )

def _hubbard_model_time_evolution_resources_trotter_4th_suzuki_trotter(
        U: float,
        tau: float,
        hopping_graph: nx.Graph,
        t: float,
        target_error: float,
        decomp: list[nx.Graph] = None,
        decomp_strategy: str = "DSATUR",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        hwp: bool = False,
) -> dict[str: int]:
    d = max(dict(hopping_graph.degree).values())
    N = len(hopping_graph.nodes)
    if nx.is_weighted(hopping_graph):
        hopping_integral = tau * max([hopping_graph[e[0]][e[1]]["weight"] for e in hopping_graph.edges])
    else:
        hopping_integral = tau

    hopping_decomposition = second_order_s1_tiling(hopping_graph, decomp, strategy=decomp_strategy)

    W = fourth_order_suzuki_trotter_split_operator_error_coefficient(U, hopping_integral, N, d)
    Nedges = len(hopping_graph.edges)
    trotter_step_gates, outer_gates = _gate_step_evolution(N, None, hopping_decomposition, Nedges, "trotter_4th", hwp)
    return hamiltonian_simulation_cost(
        t,
        target_error,
        trotter_step_gates,
        outer_gates,
        error_coefficients={5:W},
        x=x,
        synthesize_rotation_with=synthesize_rotation_with
    )

def _hubbard_model_time_evolution_resources_augmented(
        U: float,
        tau: float,
        hopping_graph: nx.Graph,
        t: float,
        target_error: float,
        decomp: list[nx.Graph] = None,
        decomp_strategy: str = "DSATUR",
        x: float | None = None,
        synthesize_rotation_with: str = "RUS",
        hwp: bool = False,
) -> dict[str: int]:
    d = max(dict(hopping_graph.degree).values())
    N = len(hopping_graph.nodes)
    if nx.is_weighted(hopping_graph):
        hopping_integral = tau * max([hopping_graph[e[0]][e[1]]["weight"] for e in hopping_graph.edges])
    else:
        hopping_integral = tau

    _func = lambda t,e: _augmented_s1_evolution_time_and_steps(
        t,
        e,
        U,
        tau,
        hopping_integral,
        d,
        N,
        hopping_graph,
        decomp,
        decomp_strategy,
    )

    hopping_decomposition = augmented_s1_tiling(hopping_graph, decomp, strategy=decomp_strategy)

    Nedges = len(hopping_graph.edges)
    Num_s1_tiles = 2 * sum([len(G.edges) for G in hopping_decomposition])

    if decomp is None:
        _decomp = _color_edge_decomp(hopping_graph, strategy=decomp_strategy)
    else:
        _decomp = decomp

    trotter_step_gates, outer_gates = _gate_step_evolution(N, _decomp, hopping_decomposition, Nedges, "augmented", hwp)
    return hamiltonian_simulation_cost(
        t,
        target_error,
        trotter_step_gates,
        outer_gates,
        trotter_steps_from_time_and_error=_func,
        x=x,
        synthesize_rotation_with=synthesize_rotation_with
    )


def _augmented_s1_evolution_time_and_steps_qpe(
        eps: float,
        U: float,
        tau: float,
        hopping_integral: float,
        d: int,
        N: int,
        hopping_graph: nx.Graph,
        decomp: list[nx.Graph] = None,
        strategy: str = "DSATUR",
        unitary_decomp: bool = True,
        bit_precision_constant: float = 0.76
) -> tuple[float, int]:
    W5, W6, W7 = fourth_order_augmented_split_operator_error_coefficients(U, hopping_integral, N, d, unitary_decomp=unitary_decomp)
    def f(t: float) -> float:
        _decomp = augmented_s1_tiling(hopping_graph, decomp = decomp, strategy=strategy, time = tau*t)
        target_graph = multiply_graph_by(hopping_graph, tau*t)
        return W5 * t**5 + W6 * t**6 + W7 * t**7 + error_between_exp_of_free_fermionic(_decomp, target_graph, scale = 1j)

    opt_func = lambda t: bit_precision_constant*np.pi/(t*eps - f(t))

    t0 = (eps / (5*W5))**(1/4)
    tm = (eps/ W5)**(1/4)
    bnds = [(0.00001,tm)]# These bounds ensure Npe > 0

    min_res = minimize_scalar(opt_func, bounds = bnds[0])

    if not min_res.success:
        warnings.warn(f"Npe optimizer failed with message: {min_res.message}")

    Npe = min_res.fun

    return  min_res.x, math.ceil(Npe)

def _augmented_s1_evolution_time_and_steps(
        t: float,
        eps: float,
        U: float,
        tau: float,
        hopping_integral: float,
        d: int,
        N: int,
        hopping_graph: nx.Graph,
        decomp: list[nx.Graph] = None,
        strategy: str = "DSATUR",
) -> tuple[float, int]:
    W5, W6, W7 = fourth_order_augmented_split_operator_error_coefficients(U, hopping_integral, N, d)
    def solve_func(n: float |list[float]) -> float:
        if not isinstance(n, float):
            return [solve_func(x) for x in n]
        _decomp = augmented_s1_tiling(hopping_graph, decomp = decomp, strategy=strategy, time = tau*t / n)
        target_graph = multiply_graph_by(hopping_graph, tau*t / n)
        return W5 * t**5 / n**4 + W6 * t**6 / n**4 + W7 * t**7 / n**6 + n*error_between_exp_of_free_fermionic(_decomp, target_graph, scale = 1j) - eps

    n0 = (W5 / eps)**(1/4) * t**(5/4)

    root, _, ier, mesg = fsolve(solve_func, n0, full_output=True)
    if not ier:
        warnings.warn(
            f"fsolve failed to converge with msg:\n{mesg}"
        )

    n = math.ceil(root[0])
    return n

def _gate_step_qpe(
        N: int,
        single_decomp: list[nx.Graph],
        hopping_decomposition: list[nx.Graph],
        Nedges: int,
        type: str,
        hamming_weight_phasing: bool
) -> dict[ResourceGate: int]:
    if hamming_weight_phasing:
        match type:
            case "trotter":
                trotter_step_gates = _hamming_weight_phasing_hopping_terms(hopping_decomposition)
                interaction_gates = HWPGate({ShiftedOccupationPair(): N}, N)
                trotter_step_gates[interaction_gates] = 1
                outer_gates = {HWPGate({
                    ShiftedOccupationPair(): N
                },N) : 1}
                return trotter_step_gates, outer_gates
            case "trotter_4th":
                trotter_step_gates = _hamming_weight_phasing_hopping_terms(hopping_decomposition)
                for x in trotter_step_gates.keys():
                    trotter_step_gates[x] *= 5
                interaction_gates = HWPGate({ShiftedOccupationPair(): N}, N)
                trotter_step_gates[interaction_gates] = 5
                outer_gates = {HWPGate({
                    ShiftedOccupationPair(): N
                },N) : 1}
                return trotter_step_gates, outer_gates
            case "augmented":
                trotter_step_gates = _hamming_weight_phasing_hopping_terms(hopping_decomposition)
                interaction_gates = HWPGate({ShiftedOccupationPair(): N}, N)
                commutator_gates_inner = _hamming_weight_phasing_commutator_term(single_decomp, type = 2)
                trotter_step_gates[interaction_gates] = 1
                trotter_step_gates.update(commutator_gates_inner)
                outer_gates = _hamming_weight_phasing_commutator_term(single_decomp, type = 1)
                return trotter_step_gates, outer_gates
    else:
        match type:
            case "trotter":
                trotter_step_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 1
                }
                for g in hopping_decomposition:
                    new_gate = ParallelizedResourceGate({FreeFermionicS1Tile() : 2*len(g.edges)})
                    trotter_step_gates[new_gate] = 1
                outer_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 1
                }
                return trotter_step_gates, outer_gates
            case "trotter_4th":
                trotter_step_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 5
                }
                for g in hopping_decomposition:
                    new_gate = ParallelizedResourceGate({FreeFermionicS1Tile() : 2*len(g.edges)})
                    trotter_step_gates[new_gate] = 5
                outer_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 1
                }
                return trotter_step_gates, outer_gates
            case "augmented":
                trotter_step_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 2
                }
                for g in hopping_decomposition:
                    new_gate = ParallelizedResourceGate({FreeFermionicS1Tile(): 2 * len(g.edges)})
                    trotter_step_gates[new_gate] = 1

                for g in single_decomp:
                    new_gate = ParallelizedResourceGate({SpinSymmetricMixedControlledHappa(): len(g.edges)})
                    trotter_step_gates[new_gate] = 1

                outer_gates = {
                    ParallelizedResourceGate({SpinSymmetricMixedControlledKappa(): len(g.edges)}) : 1 for g in single_decomp
                }
                return trotter_step_gates, outer_gates

def _gate_step_evolution(
        N: int,
        single_decomp: list[nx.Graph],
        hopping_decomposition: list[nx.Graph],
        Nedges: int,
        type: str,
        hamming_weight_phasing: bool
) -> dict[ResourceGate: int]:
    if hamming_weight_phasing:
        match type:
            case "trotter":
                trotter_step_gates = _hamming_weight_phasing_hopping_terms(hopping_decomposition)
                interaction_gates = HWPGate({ShiftedOccupationPair(): N}, N)
                trotter_step_gates[interaction_gates] = 1
                outer_gates = {HWPGate({
                    ShiftedOccupationPair(): N
                },N) : 1}
                return trotter_step_gates, outer_gates
            case "trotter_4th":
                trotter_step_gates = _hamming_weight_phasing_hopping_terms(hopping_decomposition)
                for x in trotter_step_gates.keys():
                    trotter_step_gates[x] *= 5
                interaction_gates = HWPGate({ShiftedOccupationPair(): N}, N)
                trotter_step_gates[interaction_gates] = 5
                outer_gates = {HWPGate({
                    ShiftedOccupationPair(): N
                },N) : 1}
                return trotter_step_gates, outer_gates
            case "augmented":
                trotter_step_gates = _hamming_weight_phasing_hopping_terms(hopping_decomposition)
                interaction_gates = HWPGate({ShiftedOccupationPair(): N}, N)
                commutator_gates_inner = _hamming_weight_phasing_commutator_term(single_decomp, type = 2)
                trotter_step_gates[interaction_gates] = 1
                trotter_step_gates.update(commutator_gates_inner)
                outer_gates = _hamming_weight_phasing_commutator_term(single_decomp, type = 1)
                for x in outer_gates.keys():
                    outer_gates[x] *= 2
                return trotter_step_gates, outer_gates

    else:
        match type:
            case "trotter":
                trotter_step_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 1
                }
                for g in hopping_decomposition:
                    new_gate = ParallelizedResourceGate({FreeFermionicS1Tile(): 2 * len(g.edges)})
                    trotter_step_gates[new_gate] = 1
                outer_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 1
                }
                return trotter_step_gates, outer_gates
            case "trotter_4th":
                trotter_step_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 5
                }
                for g in hopping_decomposition:
                    new_gate = ParallelizedResourceGate({FreeFermionicS1Tile(): 2 * len(g.edges)})
                    trotter_step_gates[new_gate] = 5
                outer_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 1
                }
                return trotter_step_gates, outer_gates
            case "augmented":
                trotter_step_gates = {
                    ParallelizedResourceGate({ShiftedOccupationPair(): N}): 2
                }
                for g in hopping_decomposition:
                    new_gate = ParallelizedResourceGate({FreeFermionicS1Tile(): 2 * len(g.edges)})
                    trotter_step_gates[new_gate] = 1

                for g in single_decomp:
                    new_gate = ParallelizedResourceGate({SpinSymmetricMixedControlledHappa(): len(g.edges)})
                    trotter_step_gates[new_gate] = 1

                outer_gates = {
                    ParallelizedResourceGate({SpinSymmetricMixedControlledKappa(): len(g.edges)}): 1 for g in
                    single_decomp
                }
                return trotter_step_gates, outer_gates

def _hamming_weight_phasing_hopping_terms(
        hopping_decomposition: list[nx.Graph]
) -> dict[ResourceGate, int]:
    res = {}
    for G in hopping_decomposition:
        if not nx.is_weighted(G):
            m = 2*len(G.edges)
            hwp_gate = HWPGate({
                FreeFermionicS1Tile(): m
            }, 2*m)
            res[hwp_gate] = 1
            continue

        _group_labels = []
        _group_counts = []
        for edge in G.edges:
            weight = G[edge[0]][edge[1]]["weight"]
            index = -1
            for i, c in enumerate(_group_labels):
                if abs(weight - c) < 1e-12:
                    index = i
                    break
            if index > -1:
                _group_counts[i] += 1
            else:
                _group_labels.append(weight)
                _group_counts.append(1)
        for c in _group_counts:
            hwp_gate = HWPGate({
                FreeFermionicS1Tile(): 2*c},
                4*c
            )
            res[hwp_gate] = 1
    return res

def _hamming_weight_phasing_commutator_term(
        hopping_graph_decomp: list[nx.Graph],
        type: int = 1
) -> dict[ResourceGate, int]:
    res = {}
    for G in hopping_graph_decomp:
        if not nx.is_weighted(G):
            m = len(G.edges)
            if type == 1:
                hwp_gate = HWPGate({
                    SpinSymmetricMixedControlledKappa(): m
                }, m)
                res[hwp_gate] = 1
                return res
            if type == 2:
                hwp_gate = HWPGate({
                    SpinSymmetricMixedControlledHappa(): m
                }, m)
                res[hwp_gate] = 1
                return res

        _group_labels = []
        _group_counts = []
        for edge in G.edges:
            weight = G[edge[0]][edge[1]]["weight"]
            index = -1
            for i, c in enumerate(_group_labels):
                if abs(weight - c) < 1e-12:
                    index = i
                    break
            if index > -1:
                _group_counts[i] += 1
            else:
                _group_labels.append(weight)
                _group_counts.append(1)
        for c in _group_counts:
            if type == 1:
                hwp_gate = HWPGate({
                    SpinSymmetricMixedControlledKappa(): c},
                    c
                )
                res[hwp_gate] = 1

            if type == 2:
                hwp_gate = HWPGate({
                    SpinSymmetricMixedControlledHappa(): c},
                    c
                )
                res[hwp_gate] = 1
    return res