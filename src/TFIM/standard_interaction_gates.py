import networkx as nx
import numpy as np

from src.resource_estimates.gate_costs.protocol import ResourceGate
from src.resource_estimates.gate_costs.skew_field_ising_model import (
    XField,
    ZInteraction,
    YInteraction,
    ZYInteraction,
)

def gates_from_graph_2nd(
        G: nx.Graph,
        hwp: bool = False
) -> tuple[dict[ResourceGate: int], dict[ResourceGate: int]]:
    if hwp:
        raise NotImplementedError()
    N = len(G.nodes)
    total_interactions = len(G.edges)
    return {
        ZInteraction(): total_interactions,
        XField(): N
    }, {
        XField(): N,
    }

def gates_from_graph_4th(
        G: nx.Graph,
        hwp: bool = False
) -> tuple[dict[ResourceGate: int], dict[ResourceGate: int]]:
    if hwp:
        raise NotImplementedError()
    N = len(G.nodes)
    total_interactions = len(G.edges)
    return {
        ZInteraction(): 5*total_interactions,
        XField(): 5*N
    }, {
        XField(): N,
    }

def gates_from_graph_augmented(
        G: nx.Graph,
        hwp: bool = False
) -> tuple[dict[ResourceGate: int], dict[ResourceGate: int]]:
    if hwp:
        raise NotImplementedError()
    N = len(G.nodes)
    total_interactions = len(G.edges)
    return {
        ZInteraction(): total_interactions,
        YInteraction(): total_interactions,
        XField(): 2*N
    }, {
        ZInteraction(): total_interactions,
        ZYInteraction(): 2*total_interactions
    }

def gates_from_positions_and_interaction_2nd_order(
        positions: np.ndarray,
        interaction: callable,
        tol: float = 1e-14,
        hwp: bool = False
) -> tuple[dict[ResourceGate: int], dict[ResourceGate: int]]:
    if hwp:
        raise NotImplementedError()
    N = len(positions)
    num_interactions_per_site = sum([1 for x in positions[1:] if abs(interaction(positions[0], x)) > tol])
    total_interactions = N * num_interactions_per_site // 2

    return {
        ZInteraction(): total_interactions,
        XField(): N
    },{
        XField(): N,
    }


def gates_from_positions_and_interaction_4th_order(
        positions: np.ndarray,
        interaction: callable,
        tol: float = 1e-14,
        hwp: bool = False
) -> tuple[dict[ResourceGate: int], dict[ResourceGate: int]]:
    if hwp:
        raise NotImplementedError()
    N = len(positions)
    num_interactions_per_site = sum([1 for x in positions[1:] if abs(interaction(positions[0], x)) > tol])
    total_interactions = N * num_interactions_per_site // 2

    return {
        ZInteraction(): 5*total_interactions,
        XField(): 5*N
    }, {
        XField(): N,
    }


def gates_from_positions_and_interaction_augmented(
        positions: np.ndarray,
        interaction: callable,
        tol: float = 1e-14,
        hwp: bool = False
) -> tuple[dict[ResourceGate: int], dict[ResourceGate: int]]:
    if hwp:
        raise NotImplementedError()
    N = len(positions)
    num_interactions_per_site = sum([1 for x in positions[1:] if abs(interaction(positions[0], x)) > tol])
    total_interactions = N * num_interactions_per_site // 2

    return {
        ZInteraction(): total_interactions,
        YInteraction(): total_interactions,
        XField(): 2*N
    }, {
        ZInteraction(): total_interactions,
        ZYInteraction(): 2*total_interactions
    }