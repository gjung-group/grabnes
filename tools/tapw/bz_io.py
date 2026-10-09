"""Parsers shared by the TAPW Brillouin-zone visualization tools."""

from pathlib import Path

import numpy as np


def _numbers(line):
    return [float(value) for value in line.split()]


def parse_bz_file(path):
    """Return vertices, K point, and reciprocal vectors from a debug file."""
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    vertices = []
    k_point = None
    reciprocal = []

    for index, line in enumerate(lines):
        marker = line.strip()
        if marker == "# BZ vertices (x, y):":
            vertices = [_numbers(value)[:2] for value in lines[index + 1 : index + 7]]
        elif marker == "# K-point:":
            k_point = _numbers(lines[index + 1])[:2]
        elif marker == "# Reciprocal lattice vectors b1, b2:":
            reciprocal = [_numbers(value)[:2] for value in lines[index + 1 : index + 3]]

    if len(vertices) != 6:
        raise ValueError(f"{path}: expected six Brillouin-zone vertices")
    if len(reciprocal) != 2:
        raise ValueError(f"{path}: expected two reciprocal lattice vectors")

    return (
        np.asarray(vertices),
        None if k_point is None else np.asarray(k_point),
        np.asarray(reciprocal[0]),
        np.asarray(reciprocal[1]),
    )


def load_table(path, *, required=False):
    """Load a numeric debug table, returning None for an optional missing file."""
    path = Path(path)
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return None
    return np.atleast_2d(np.loadtxt(path))


def closed_polygon(vertices):
    """Append the first vertex so Matplotlib draws a closed polygon."""
    return np.vstack((vertices, vertices[0]))
