import networkx as nx
import matplotlib.pyplot as plt

def generate_bethe_lattice(generations: int, coordination_number: int) -> nx.graph:
    """
    Generates a Bethe lattice using a recursive tree approach.
    """
    G = nx.Graph()
    node_counter = 0

    def add_branches(parent_id, current_generation, parent_branch):
        nonlocal node_counter

        # Determine how many branches this node needs
        num_neighbors = coordination_number if current_generation == 0 else coordination_number - 1

        for _ in range(num_neighbors):
            if current_generation < generations:
                node_counter += 1
                child_id = node_counter
                G.add_edge(parent_id, child_id)
                # Recursively add branches to this child
                add_branches(child_id, current_generation + 1, parent_id)

    # Step 1: Create the center node (generation 0)
    center_node = 0
    add_branches(center_node, 0, None)

    return G
