# GRABNES test harness

This directory contains a small regression harness for the canonical solver in
[`../lanczosKuboCode`](../lanczosKuboCode/). It holds **no copy of the solver
and no copy of the examples**: every run compiles the current
`lanczosKuboCode/` sources out of tree, runs the inputs in
[`../examples`](../examples/), and compares the results with
`examples/*/reference/`.

## Commands

Both scripts can be started from any directory.

```sh
./grabnes_testrun/smoke_test.sh      # build + pristine graphene bands (about 1.5 minutes)
./grabnes_testrun/run_examples.sh    # build + all four public examples + solver checks
```

The exit status is `0` when every check passes, `1` when the build, a run, or
a comparison fails, and `2` when a prerequisite is missing (compiler, GNU
Make, Python 3, input file).

Useful variations:

```sh
# Checked build: -O0 -g -fcheck=all -fbacktrace
./grabnes_testrun/run_examples.sh --debug

# Another build configuration, for example Intel Fortran + Intel MPI + MKL
./grabnes_testrun/run_examples.sh --make-sys grabnes_testrun/config/intel.make.sys

# Keep the build so that the next run only recompiles what changed
./grabnes_testrun/run_examples.sh --work-dir /path/to/scratch/grabnes-test

# Test an executable that is already built
GRABNES_BIN=/path/to/grabnes ./grabnes_testrun/run_examples.sh

# Start the executable through an MPI launcher or the batch system
GRABNES_LAUNCHER="mpirun -np 1" ./grabnes_testrun/run_examples.sh
```

Run `./grabnes_testrun/smoke_test.sh --help` for the complete list.

The examples use exact diagonalization and must be started with **one MPI
process**; see the limitations in
[`docs/development/cluster-build-and-validation.md`](../docs/development/cluster-build-and-validation.md).

## Requirements

- GNU Make, a Fortran compiler with an MPI wrapper, BLAS/LAPACK, and ARPACK,
  as required by the solver itself;
- Python 3 (standard library only) for the numerical comparison;
- a POSIX shell.

## Build configuration

The harness uses, in this order: `--debug` or `--make-sys FILE`, the
`GRABNES_MAKE_SYS` environment variable, your local
`lanczosKuboCode/make.sys` if it exists, and otherwise the tracked
`lanczosKuboCode/make.sys.example`.

| File | Purpose |
| --- | --- |
| `config/gfortran.debug.make.sys` | GNU Fortran with runtime checks (`--debug`) |
| `config/intel.make.sys` | Intel Fortran (`ifort`), Intel MPI, MKL, OpenMP |
| `config/macos-homebrew.make.sys` | Homebrew GNU toolchain on macOS; historical, see `BUILD_REPORT.md` |

## Where files go

Nothing is written into `lanczosKuboCode/` or `examples/`. Each run creates a
work directory, by default a new one under `$TMPDIR` (or `/tmp`):

```text
<work directory>/
├── build/      object files, modules, generated sources
├── bin/        the grabnes executable
├── build.log   complete compiler output
└── run/<check>/     copied input, job.out, job.err, and all solver output
```

A temporary work directory is removed when every check passes and kept, with
its path printed, when something fails. A directory given with `--work-dir`
(or `GRABNES_TEST_DIR`) is always kept.

## What is checked

For every example the harness requires that the solver exits with status 0,
reports `0 errors, 0 warnings`, and writes its output file; `compare_output.py`
then checks that

- result and reference have the same number of lines and of values per line;
- every value agrees within `atol + rtol * |reference|`;
- band files are internally consistent (header against data size, ascending
  eigenvalues at each k-point) and DOS files have an increasing energy grid and
  no negative density.

| Output | atol | rtol | Reason |
| --- | --- | --- | --- |
| `generate.bands` | 2e-6 | 0 | written with six decimals; allows one last-digit rounding flip |
| `generate.diag.DOS` | 1e-12 | 1e-9 | written with 17 significant digits; limited by BLAS/LAPACK and math-library rounding |

`run_examples.sh` then runs checks that do not depend on stored solver output:

| Check | Requirement |
| --- | --- |
| `dirac_point_degeneracy` | The four Dirac states of example 03 at the moire K point form two pairs degenerate within 1e-6 eV (a symmetry the historical neighbor search violated by 1e-4 eV). |
| `hamiltonian_example03`, `hamiltonian_f2g2`, `hamiltonian_small_cell` | The solver is run with `WriteDataFiles .true.` on example 03, on its F2G2 variant (`TB.NeighLevels 5`), and on a 28-atom cell smaller than the interlayer search radius. [`tools/hamiltonian/verify_tables.py`](../tools/hamiltonian/) compares the tables it writes with a Hamiltonian enumerated independently from the atomic positions: no missing, extra, duplicate, or unpaired entry; matrix elements, Hermiticity, and eigenvalues at Gamma, K, M, generic and random k-points within 1e-9 eV; bands within 1e-6 eV of the six-decimal band file. Needs NumPy; skipped without it. |
| `neighlevels_0`, `neighlevels_9` | `TB.NeighLevels` outside 1 to 8 must be refused with its error message. |
| `legacy_cutatnn3` | The deprecated `Neigh.CutAtNN3 .true.` with `TB.NeighLevels 5` must give the same band file as `TB.NeighLevels 3`. |
| `kubo_graphene_dos` | The stochastic Kubo (Lanczos recursion) DOS of a 180000-atom graphene cell with a fixed seed must match the exact broadened DOS in `kubo_graphene_dos/reference/exact_dos.dat` (analytic, written by `make_reference.py`): maximum deviation below 0.03, rms below 0.006. About 2 s and 0.5 GB. |
| `two_mpi_processes` | Started with `mpirun -np 2`, the solver must refuse to run with its single-process error. Skipped when no `mpirun` is found; set `GRABNES_MPIRUN` to use another launcher. |

The Hamiltonian checks use tolerances set by floating-point precision
(observed: 4e-16 eV with GNU Fortran, 2e-14 eV with Intel Fortran) and by the
six decimals of the band file (observed: 5.0e-7 eV). The Kubo check is
statistical because the random numbers depend on the compiler; its limits are
about twice the deviation observed with GNU and Intel Fortran.

Byte identity is reported for information. `--exact` additionally requires it
for the band files, which have been byte-identical for every compiler tested
so far. If a comparison fails, find out why before changing a tolerance or a
reference file.

The header of every run records the Git commit, the build configuration, the
compiler version and flags, and the md5 sum of the tested executable.

## Python tests

The TAPW linear-algebra tests are independent of the Fortran solver:

```sh
python3 -m pytest tests/tapw
```
