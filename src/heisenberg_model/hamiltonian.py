import networkx as nx

from src.bch_formula.unitary_modified_bch_terms import unitary_modified_bch_3_v1

from qiskit.quantum_info import SparsePauliOp

def _fcom(H1: SparsePauliOp, H2: SparsePauliOp) -> SparsePauliOp:
    return (H1 @ H2 - H2 @ H1).simplify().chop(tol = 0.1)


def heisenberg_terms(graph: nx.Graph, Js: list[float] | float, h: float = None) -> tuple[SparsePauliOp, ...]:
    X_terms = []
    Y_terms = []
    Z_terms = []

    if isinstance(Js, float):
        Js = [Js, Js, Js]

    _graph = nx.convert_node_labels_to_integers(graph)

    N = len(_graph.nodes)
    for edge in _graph.edges:
        Xstring_ls = ["X" if i in edge else "I" for i in range(N)]
        Ystring_ls = ["Y" if i in edge else "I" for i in range(N)]
        Zstring_ls = ["Z" if i in edge else "I" for i in range(N)]

        X_terms.append("".join(Xstring_ls))
        Y_terms.append("".join(Ystring_ls))
        Z_terms.append("".join(Zstring_ls))

    X_coeffs = [Js[0] for _ in X_terms]
    Y_coeffs = [Js[1] for _ in Y_terms]
    Z_coeffs = [Js[2] for _ in Z_terms]
    if h is not None:
        Z_terms.extend(
            ["I" * i + "Z" + "I" * (N- i - 1) for i in range(N)]
        )
        Z_coeffs.extend([h for _ in range(N)])

    X = SparsePauliOp(X_terms, X_coeffs)
    Y = SparsePauliOp(Y_terms, Y_coeffs)
    Z = SparsePauliOp(Z_terms, Z_coeffs)
    return X, Y, Z

def heisenberg_hamiltonian_with_correctors(graph: nx.Graph, Js: list[float] | float, h: float = None) -> tuple[SparsePauliOp, ...]:
    X,Y,Z = heisenberg_terms(graph, Js, h)
    HS = [Z,X,Y]
    F, C = unitary_modified_bch_3_v1(HS, _fcom)
    return *HS, F,C