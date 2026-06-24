"""
Functions for generating circuits to simulate hubbard models on the 2D square lattice with fully symmetric coulomb and hopping hamiltonians.

Second order circuits are copied directly from:
    Earl T Campbell 2022 Quantum Sci. Technol. 7 015007

In all calculations we assume a periodic lattice
"""

from openfermion.ops import FermionOperator

def get_fermionic_operator(U: float, tau: float, L: int) -> FermionOperator:
    ...

def _second_order_error_coefficient(U: float, tau: float, L: int) -> float:
    ...

def _fourth_order_suzuki_trotter_error_coefficient(U: float, tau: float, L: int) -> float:
    ...

def _fourth_order_augmented_error_coefficient(U: float, tau: float, L: int) -> float:
    ...