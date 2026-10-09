"""Checks for the array conventions used by the Fortran TAPW projection."""

import numpy as np
import scipy.sparse as sp


def build_tapw_basis(x, y, labels, gx, gy, label_count):
    """Build X with Fortran ordering: column = G index * labels + label."""
    orbital_count = len(x)
    counts = np.bincount(labels, minlength=label_count)
    if np.any(counts == 0):
        raise ValueError("each TAPW label must contain at least one orbital")

    basis = np.zeros((orbital_count, len(gx) * label_count), dtype=complex)
    for g_index, (g_x, g_y) in enumerate(zip(gx, gy)):
        phase = np.exp(1j * (g_x * x + g_y * y))
        for label_index in range(label_count):
            column = g_index * label_count + label_index
            mask = labels == label_index
            basis[mask, column] = phase[mask] / np.sqrt(counts[label_index])
    return basis


def test_fortran_column_ordering_and_normalization():
    rng = np.random.default_rng(123)
    orbital_count = 100
    g_count = 25
    label_count = 2
    labels = np.tile(np.arange(label_count), orbital_count // label_count)
    rng.shuffle(labels)

    basis = build_tapw_basis(
        rng.uniform(-10, 10, orbital_count),
        rng.uniform(-10, 10, orbital_count),
        labels,
        rng.uniform(-2, 2, g_count),
        rng.uniform(-2, 2, g_count),
        label_count,
    )

    assert basis.shape == (orbital_count, g_count * label_count)
    np.testing.assert_allclose(np.linalg.norm(basis, axis=0), 1.0, atol=1e-13)
    for column in range(basis.shape[1]):
        expected_label = column % label_count
        assert np.count_nonzero(basis[labels != expected_label, column]) == 0


def test_two_step_projection_matches_dense_expression():
    rng = np.random.default_rng(42)
    orbital_count = 40
    projected_count = 12
    dense = rng.normal(size=(orbital_count, orbital_count))
    dense = (dense + dense.T) / 2
    hamiltonian = sp.csr_matrix(dense)
    basis = rng.normal(size=(orbital_count, projected_count)) + 1j * rng.normal(
        size=(orbital_count, projected_count)
    )

    two_step = basis.conj().T @ hamiltonian.dot(basis)
    direct = basis.conj().T @ hamiltonian.toarray() @ basis

    np.testing.assert_allclose(two_step, direct, rtol=1e-13, atol=1e-12)
    np.testing.assert_allclose(two_step, two_step.conj().T, atol=1e-12)
