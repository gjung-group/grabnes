"""Unit tests for sparse-to-dense TAPW projection mathematics."""

import numpy as np
import scipy.sparse as sp


def hermitian_sparse_matrix(size, density, seed):
    rng = np.random.default_rng(seed)
    matrix = sp.random(
        size,
        size,
        density=density,
        format="csr",
        dtype=complex,
        random_state=rng,
        data_rvs=lambda count: rng.normal(size=count) + 1j * rng.normal(size=count),
    )
    matrix = (matrix + matrix.conj().T) / 2
    matrix.setdiag(rng.uniform(-2, 2, size))
    return matrix.tocsr()


def test_sparse_two_step_transform_matches_direct_transform():
    rng = np.random.default_rng(123)
    hamiltonian = hermitian_sparse_matrix(size=50, density=0.15, seed=42)
    basis = rng.normal(size=(50, 20)) + 1j * rng.normal(size=(50, 20))

    intermediate = hamiltonian.dot(basis)
    projected = basis.conj().T @ intermediate
    direct = basis.conj().T @ hamiltonian.toarray() @ basis

    assert intermediate.shape == (50, 20)
    assert projected.shape == (20, 20)
    np.testing.assert_allclose(projected, direct, rtol=1e-13, atol=1e-12)
    np.testing.assert_allclose(projected, projected.conj().T, atol=1e-12)
