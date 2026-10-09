"""Checks for projection and reconstruction with an orthonormal TAPW basis."""

import numpy as np

from test_transform import hermitian_sparse_matrix


def orthonormal_basis(rows, columns, seed):
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(rows, columns)) + 1j * rng.normal(size=(rows, columns))
    basis, _ = np.linalg.qr(raw)
    return basis[:, :columns]


def test_orthonormal_basis_and_projected_reconstruction():
    size = 30
    projected_size = 15
    basis = orthonormal_basis(size, projected_size, seed=42)
    hamiltonian = hermitian_sparse_matrix(size=size, density=0.2, seed=17).toarray()

    np.testing.assert_allclose(
        basis.conj().T @ basis, np.eye(projected_size), atol=1e-12
    )

    projected = basis.conj().T @ hamiltonian @ basis
    reconstructed = basis @ projected @ basis.conj().T
    projected_original = (
        basis @ (basis.conj().T @ hamiltonian @ basis) @ basis.conj().T
    )

    np.testing.assert_allclose(projected, projected.conj().T, atol=1e-12)
    np.testing.assert_allclose(reconstructed, projected_original, atol=1e-12)
