# GRABNES test harness

This directory contains a small regression harness for the canonical solver in
[`../../lanczosKuboCode`](../../lanczosKuboCode/). It holds **no copy of the solver
and no copy of the examples**: every run compiles the current
`lanczosKuboCode/` sources out of tree, runs the inputs in
[`../../examples`](../../examples/), and compares the results with
`examples/*/reference/`.

## Commands

Both scripts can be started from any directory.

```sh
./tests/regression/smoke_test.sh      # build + pristine graphene bands (about 1.5 minutes)
./tests/regression/run_examples.sh    # build + all four public examples + solver checks
```

The exit status is `0` when every check passes, `1` when the build, a run, or
a comparison fails, and `2` when a prerequisite is missing (compiler, GNU
Make, Python 3, input file).

Useful variations:

```sh
# Checked build: -O0 -g -fcheck=all -fbacktrace
./tests/regression/run_examples.sh --debug

# Another build configuration, for example Intel Fortran + Intel MPI + MKL
./tests/regression/run_examples.sh --make-sys tests/regression/config/intel.make.sys

# Keep the build so that the next run only recompiles what changed
./tests/regression/run_examples.sh --work-dir /path/to/scratch/grabnes-test

# Test an executable that is already built
GRABNES_BIN=/path/to/grabnes ./tests/regression/run_examples.sh

# Start the executable through an MPI launcher or the batch system
GRABNES_LAUNCHER="mpirun -np 1" ./tests/regression/run_examples.sh
```

Run `./tests/regression/smoke_test.sh --help` for the complete list.

The examples use exact diagonalization and must be started with **one MPI
process**; see the limitations in
[`docs/development/cluster-build-and-validation.md`](../../docs/development/cluster-build-and-validation.md).

## Survey of the model switches

`model_survey.py` runs one or more small cases (tens to hundreds of atoms) for
every model switch of the solver; the cases are listed in `model_cases.py`. The
solver only assembles the Hamiltonian, or diagonalises it at two k-points for
terms that act on the spin. For every case the script records how the run ended

| Status | Meaning |
| --- | --- |
| `ok` | ended normally; the tables are finite and Hermitian |
| `INERT` | as `ok`, but identical to the base input: the switch had no effect in this configuration |
| `NONHERM` | the two directions of a bond differ, or a reverse entry is missing |
| `REFUSED` | the solver stopped with its own error message |
| `CRASH`, `NAN` | abnormal end, or non-finite values |

and a fingerprint of the Hamiltonian (or of the eigenvalues) that does not
depend on the order of the neighbor list. `run_examples.sh` compares both with
`model_survey_reference.json`, so a change in any model is detected. A case
that is `REFUSED`, `INERT`, or `NONHERM` in the reference documents the present
state of that switch; it is not a statement that the switch is correct.

```sh
# table for one family of cases, with a build that checks array bounds
python3 tests/regression/model_survey.py --bin /path/to/grabnes --only 'xyz4*' --report survey.md
# after an intended change of a model: inspect the differences, then record them
python3 tests/regression/model_survey.py --bin /path/to/grabnes --check
python3 tests/regression/model_survey.py --bin /path/to/grabnes --update-reference
```

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
| `hamiltonian_example03`, `hamiltonian_f2g2`, `hamiltonian_small_cell` | The solver is run with `WriteDataFiles .true.` on example 03, on its F2G2 variant (`TB.NeighLevels 5`), and on a 28-atom cell smaller than the interlayer search radius. [`tools/hamiltonian/verify_tables.py`](../../tools/hamiltonian/) compares the tables it writes with a Hamiltonian enumerated independently from the atomic positions: no missing, extra, duplicate, or unpaired entry; matrix elements, Hermiticity, and eigenvalues at Gamma, K, M, generic and random k-points within 1e-9 eV; bands within 1e-6 eV of the six-decimal band file. Needs NumPy; skipped without it. |
| `koshino_intralayer_default` | With `KoshinoIntralayer .true.` and no `vpppi0` in the input, the tables must match the independent two-center model with 2.7 eV, and the log must report 2.7 eV; the F2G2 run must report 3.5 eV. This protects the model-dependent convention for `vpppi0`. |
| `hbn_monolayer` | hBN monolayer against the independent model (t = 3.0294 eV, on-site energies 3.09 and -1.89 eV). |
| `twisted_bulk` | A twisted cell periodic along z: the stored tables must be complete, symmetric, Hermitian, and reproduce the bands (structure only; no model). |
| `supercell_neighlist` | The legacy `NeighList` routine must give the same band file as the default search for a 32-atom graphene cell. |
| `legacy_fastnn`, `legacy_notrectangle` | `Neigh.fastNN` and `Neigh.fastNNnotsquareNotRectangle` must be refused. |
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
