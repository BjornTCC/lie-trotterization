from __future__ import annotations

from copy import deepcopy

import numpy as np

class Kappa:
    def __init__(
            self,
            i: int,
            j: int,
            dimension: int = -1
    ) -> None:

        if i == j:
            raise ValueError(f"Equal rotation indices not allowed: i = j = {i}")

        self.i = i
        self.j = j
        if dimension == -1:
            self.dimension = max(self.i, self.j)
        else:
            self.dimension = dimension
            if self.dimension < max(self.i, self.j):
                raise ValueError(
                    f"dimension {self.dimension} was smaller than a rotation index, (i,j) = ({self.i},{self.j})")

        self._matrix_representation = None

    @property
    def matrix_representation(self) -> np.ndarray:
        if self._matrix_representation is None:
            self._matrix_representation = self._construct_matrix_representation()
        return self._matrix_representation

    def _construct_matrix_representation(self) -> np.ndarray:
        res = np.zeros((self.dimension, self.dimension))

        res[self.i, self.j] = 1
        res[self.j, self.i] = -1
        return res

class Happa:

    def __init__(
            self,
            i: int,
            j: int,
            dimension: int = -1
    ) -> None:

        self.i = i
        self.j = j
        self.is_diagonal = (self.i == self.j)

        if dimension == -1:
            self.dimension = max(self.i, self.j) + 1
        elif dimension <= max(self.i, self.j):
            raise ValueError(
                f"dimension {self.dimension} was smaller than a rotation index, (i,j) = ({self.i},{self.j})")
        else:
            self.dimension = dimension

        self._matrix_representation = None

    def _construct_matrix_representation(self) -> np.ndarray:
        res = np.zeros((self.dimension, self.dimension), dtype=np.complex128)
        if self.is_diagonal:
            res[i, i] = 1j
            return res
        res[i, j] = 1j
        res[j, i] = 1j
        return res

    @property
    def matrix_representation(self) -> np.ndarray:
        if self._matrix_representation is None:
            self._matrix_representation = self._construct_matrix_representation()
        return self._matrix_representation



class PlaneRotation:

    def __init__(
            self,
            i: int,
            j: int,
            theta: float,
            dimension: int = -1
    ) -> None:

        if i == j:
            raise ValueError(f"Equal rotation indices not allowed: i = j = {i}")

        self.i = i
        self.j = j
        self.theta = theta
        self._cos = np.cos(self.theta)
        self._sin = np.sin(self.theta)
        if dimension == -1:
            self.dimension = max(self.i, self.j) + 1
        else:
            self.dimension = dimension
            if self.dimension < max(self.i, self.j) + 1:
                raise ValueError(
                    f"dimension {self.dimension} was smaller than a rotation index, (i,j) = ({self.i},{self.j})")

        self._matrix_representation = None

    @property
    def matrix_representation(self) -> np.ndarray:
        if self._matrix_representation is None:
            self._matrix_representation = self._construct_matrix_representation()
        return self._matrix_representation

    def _construct_matrix_representation(self) -> np.ndarray:
        res = np.identity(self.dimension)

        res[self.i, self.i] = self._cos
        res[self.j, self.j] = self._cos

        res[self.i, self.j] = self._sin
        res[self.j, self.i] = -self._sin
        return res

    def kappa(self) -> Kappa:
        return Kappa(self.i, self.j, self.dimension)

    def right_mult(self, matrix: np.ndarray, inplace: bool = False) -> np.ndarray | None:
        if inplace:
            return self._right_mult(matrix)

        res = deepcopy(matrix)
        self._right_mult(res)
        return res

    def left_mult(self, matrix: np.ndarray, inplace: bool = False) -> np.ndarray | None:
        if inplace:
            return self._left_mult(matrix)

        res = deepcopy(matrix)
        self._left_mult(res)
        return res

    def _right_mult(self, matrix: np.ndarray) -> None:
        a_li = deepcopy(matrix[:, self.i])
        a_lj = deepcopy(matrix[:, self.j])
        for l in range(self.dimension):
            matrix[l, self.i] = self._cos * a_li[l] - self._sin * a_lj[l]
            matrix[l, self.j] = self._cos * a_lj[l] + self._sin * a_li[l]

    def _left_mult(self, matrix: np.ndarray) -> None:
        a_il = deepcopy(matrix[self.i, :])
        a_jl = deepcopy(matrix[self.j, :])
        for l in range(self.dimension):
            matrix[self.i, l] = self._cos * a_il[l] + self._sin * a_jl[l]
            matrix[self.j, l] = self._cos * a_jl[l] - self._sin * a_il[l]

    def inverse(self) -> PlaneRotation:
        return PlaneRotation(self.i, self.j, -self.theta, self.dimension)

class HappaRotation:
    def __init__(
            self,
            i: int,
            j: int,
            theta: float,
            dimension: int = -1
    ) -> None:

        self.i = i
        self.j = j
        self.is_diagonal = (self.i == self.j)
        self.theta = theta

        if dimension == -1:
            self.dimension = max(self.i, self.j) + 1
        elif dimension <= max(self.i, self.j):
            raise ValueError(
                f"dimension {self.dimension} was smaller than a rotation index, (i,j) = ({self.i},{self.j})")
        else:
            self.dimension = dimension

        self._cos = np.cos(self.theta)
        self._sin = np.sin(self.theta)

        self._matrix_representation = None

    def _construct_matrix_representation(self) -> np.ndarray:
        res = np.identity(self.dimension, dtype=np.complex128)
        if self.is_diagonal:
            res[self.i, self.i] = np.exp(self.theta * 1j)
            return res

        res[self.i, self.i] = self._cos
        res[self.j, self.j] = self._cos

        res[self.i, self.j] = 1j * self._sin
        res[self.j, self.i] = 1j * self._sin

        return res

    @property
    def matrix_representation(self) -> np.ndarray:
        if self._matrix_representation is None:
            self._matrix_representation = self._construct_matrix_representation()
        return self._matrix_representation

    def happa(self) -> Happa:
        return Happa(self.i, self.j, self.dimension)

    def inverse(self) -> HappaRotation:
        return HappaRotation(self.i, self.j, -self.theta, self.dimension)

def matrix_is_real_rotation(matrix: np.ndarray, tol: float = 1e-12) -> tuple[bool, int]:
    # check shape

    if len(matrix.shape) != 2:
        return False, 1

    if matrix.shape[0] != matrix.shape[1]:
        return False, 2

    # check determinant

    if abs(np.linalg.det(matrix) - 1.0) > tol:
        return False, 3

    # Check real
    if not np.isreal(matrix).all():
        return False, 4

    # check orthogonality

    matrix_T = np.transpose(matrix)

    identity = np.identity(matrix.shape[0])

    if not np.allclose(matrix_T @ matrix, identity, atol = tol) and np.allclose(matrix @ matrix_T, identity, atol = tol):
        return False, 5

    return True, 0

def matrix_is_special_unitary(matrix: np.ndarray, tol: float = 1e-12) -> tuple[bool, int]:
    # check shape

    if len(matrix.shape) != 2:
        return False, 1

    if matrix.shape[0] != matrix.shape[1]:
        return False, 2

    # check determinant

    if abs(np.linalg.det(matrix) - 1.0) > tol:
        return False, 3

    # check orthogonality

    matrix_T = np.matrix(matrix).H

    identity = np.identity(matrix.shape[0])

    if not np.allclose(matrix_T @ matrix, identity, atol = tol) and np.allclose(matrix @ matrix_T, identity, atol = tol):
        return False, 4

    return True, 0

def check_matrix_is_imag_generator(matrix: np.ndarray, tol: float = 1e-12) -> tuple[bool, int]:
    # check shape

    if len(matrix.shape) != 2:
        return False, 1

    if matrix.shape[0] != matrix.shape[1]:
        return False, 2

    # check determinant

    if abs(np.linalg.trace(matrix)) > tol:
        return False, 3

    # check imaginary:
    if not np.isreal(-1j*matrix).all():
        return False, 4

    # check anti-hermitian

    matrix_T = np.matrix(matrix).H

    if not np.allclose(matrix_T, -matrix, atol = tol):
        return False, 5

    return True, 0

def check_matrix_is_generator(matrix: np.ndarray, tol: float = 1e-12) -> tuple[bool, int]:
    # check shape

    if len(matrix.shape) != 2:
        return False, 1

    if matrix.shape[0] != matrix.shape[1]:
        return False, 2

    # check determinant

    if abs(np.linalg.trace(matrix)) > tol:
        return False, 3

    # check anti-hermitian

    matrix_T = np.matrix(matrix).H

    if not np.allclose(matrix_T, -matrix, atol=tol):
        return False, 4

    return True, 0