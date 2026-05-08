from copy import deepcopy

import numpy as np
import scipy as sp

from rotation_utils import (
    PlaneRotation,
    HappaRotation,
    matrix_is_real_rotation,
    matrix_is_special_unitary,
    check_matrix_is_imag_generator,
    check_matrix_is_generator
)

def decompose_real_rotation_into_plane_rotations(
        rotation: np.ndarray, check: bool = True, tol: float = 1e-12, generator: bool = False
) -> list[PlaneRotation]:

    if generator:
        rotation = sp.linalg.expm(rotation)

    if check:
        check, code = matrix_is_real_rotation(rotation, tol=tol)
        if not check:
            raise ValueError(f"Matrix not in SO(n), code: {code}")

    _matrix = deepcopy(rotation)
    dim = rotation.shape[0]
    res = []

    for i in range(dim):
        for j in range(i + 1, dim):
            if abs(_matrix[i, j]) < tol:
                continue
            theta = -np.atan2(_matrix[i, j], _matrix[i, i])
            rot = PlaneRotation(i, j, theta, dim)
            rot.right_mult(_matrix, inplace=True)
            res.append(PlaneRotation(i, j, -theta, dim))
    return res[::-1]

def decompose_imaginary_rotation_into_happa_rotations(
        rotation: np.ndarray, check: bool = True, tol: float = 1e-12, generator: bool = True
) -> list[HappaRotation]:

    if check and not generator:
        check, code = matrix_is_special_unitary(rotation, tol=tol)
        if not check:
            raise ValueError(f"Matrix not in SU(n), code: {code}")

    if check and generator:
        check, code = check_matrix_is_imag_generator(rotation, tol=tol)
        if not check:
            raise ValueError(f"Matrix not in SU(n), code: {code}")

    if not generator:
        generator = sp.linalg.logm(rotation)
    else:
        generator = rotation

    generator_r = -1j * generator

    D, V = sp.linalg.eigh(generator_r)

    if np.linalg.det(V) < 0:
        V *= -1
    V = np.real_if_close(V)

    rot1 = decompose_real_rotation_into_plane_rotations(V, check=False, generator=False)
    diag = [HappaRotation(i,i,d, dimension=rotation.shape[0]) for i,d in enumerate(D)]
    return rot1 + diag + [R.inverse() for R in rot1[::-1]]

def decompose_unitary_into_rotations(
    rotation: np.ndarray, check: bool = True, tol: float = 1e-12, generator: bool = True
) -> list[HappaRotation | PlaneRotation]:

    if check and not generator:
        check, code = matrix_is_special_unitary(rotation, tol=tol)
        if not check:
            raise ValueError(f"Matrix not in SU(n), code: {code}")

    if check and generator:
        check, code = check_matrix_is_generator(rotation, tol=tol)
        if not check:
            raise ValueError(f"Matrix not in SU(n), code: {code}")

    if not generator:
        generator = sp.linalg.logm(rotation)
    else:
        generator = rotation

    raise ValueError("Not fully implemented")

if __name__ == "__main__":
    # Check code
    def compare_decomposition(
            original: np.ndarray,
            decomp: list[any],
            generator: bool = False,
            verbose: bool = False
    ) -> float:
        if generator:
            comp_matrix = sp.linalg.expm(original)
        else:
            comp_matrix = original

        res_matrix = np.identity(original.shape[0])

        for i,rot in enumerate(decomp):
            if verbose:
                print(f"\nRotation {i+1} /{len(decomp)}")
                print(rot.i, rot.j, rot.theta)
                print(rot.matrix_representation)
            res_matrix = res_matrix @ rot.matrix_representation

        if verbose:
            print("\nOriginal:")
            print(comp_matrix)
            print("Result:")
            print(res_matrix)

        return np.linalg.norm(comp_matrix - res_matrix)

    d = 100

    herm_mat = np.random.rand(d,d)
    herm_mat = 1j*(np.matrix(herm_mat) + np.matrix(herm_mat).H - 2*np.trace(herm_mat) * np.identity(d)/d)

    decomp = decompose_imaginary_rotation_into_happa_rotations(herm_mat, check=False)

    print(compare_decomposition(herm_mat, decomp, generator=True, verbose=False))