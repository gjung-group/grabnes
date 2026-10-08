# Hamiltonian verification

`verify_tables.py` compares the tight-binding Hamiltonian assembled by GRABNES
with one built independently in Python (NumPy required).

1. Run a band calculation with `WriteDataFiles .true.` in a scratch directory.
   The solver then writes its neighbor list (`v`), the lattice translation and
   displacement of every entry (`neighCell.dat`, `neighD.dat`), the hopping
   values (`<prefix>.s.mag`), and the positions, cell, and on-site energies.
2. Pass the directory and the intralayer matrix elements of the model, in eV
   and one per neighbor shell (as many as `TB.NeighLevels`):

```sh
# example 03 (TB.NeighLevels 1)
python3 tools/hamiltonian/verify_tables.py /path/to/run --intralayer=-2.9888

# default F2G2 intralayer model (TB.NeighLevels 5)
python3 tools/hamiltonian/verify_tables.py /path/to/run \
    --intralayer=-2.9888,0.2354,-0.1877,0,0.0633

# pristine graphene, example 01 (TB.Hopping 3.1, nonBulkSmall, no interlayer search)
python3 tools/hamiltonian/verify_tables.py /path/to/run --g0 3.1 --intralayer=-3.1 \
    --interlayer-cutoff 1.0 --periodic-z

# original Koshino model (KoshinoIntralayer .true., TB.NeighLevels 5): two-center
# form within the layers as well, with the Moon-Koshino value of vpppi0
python3 tools/hamiltonian/verify_tables.py /path/to/run --koshino-intralayer 5 --vpppi0 2.7

# hBN monolayer: on-site energies by species number
python3 tools/hamiltonian/verify_tables.py /path/to/run --g0 3.1 --intralayer=-3.0294 \
    --onsite 3:3.09,4:-1.89 --interlayer-cutoff 1.0 --periodic-z

# no model available: geometry, translations, reverse partners and Hermiticity only,
# with a listing of the stored elements by species pair and distance
python3 tools/hamiltonian/verify_tables.py /path/to/run --structure-only 1
```

The second Hamiltonian uses only the atomic positions and these parameters:
every pair inside the search radii is enumerated exhaustively over lattice
translations and its matrix element is computed from the shell value or from
the two-center interlayer formula. The script then reports, with residuals,

- whether each stored displacement equals `r_m + R - r_i`;
- missing, extra, duplicate, and unpaired entries of the solver's list;
- the largest difference between stored and independent matrix elements;
- Hermiticity and agreement of `H(k)` at Gamma, K, M, generic, and seeded
  random k-points;
- invariance under reciprocal lattice vectors and under moving atoms to other
  unit cells (a unitary change of basis);
- the difference between the independent eigenvalues and the band file.

It exits with status 1 if a residual exceeds its limit (`--tol-matrix`,
default 1e-9 eV; `--tol-bands`, default 1e-6 eV, the resolution of the band
file). Options exist for the interlayer parameters, the band path
(`--path`), z-periodic systems, and the random seed; see `--help`.

What agreement means: the solver assembles the Hamiltonian that the stated
pair-selection rules and parameters define. It is not a statement about the
physical adequacy of those rules or parameters. The Bloch convention and the
results for the public examples are described in
[`docs/development/cluster-build-and-validation.md`](../../docs/development/cluster-build-and-validation.md).
The regression suite runs the script for three cases
(`tests/regression/run_examples.sh`).
