# Linux cluster build and numerical validation

This page records how the canonical solver (`lanczosKuboCode/`) was built and
validated on a Linux cluster in October 2026, which defects had to be fixed
first, and what remains unverified. It replaces the Apple Silicon report of
September 2026 (`grabnes_testrun/BUILD_REPORT.md`) as the reference for the
current source tree.

Everything below was measured; nothing is carried over from the earlier report.

## Baseline

| Item | Value |
| --- | --- |
| Baseline commit | `941e8598b4ace297fb694f947244cc586556f3cd` (`main` at the time) |
| Working branch | `nicolas/development` |
| Date | 2026-10-08 |
| System | Fedora Linux 36, x86_64, kernel 6.6, cluster login node |
| GNU Make | 4.3 |
| GNU toolchain | GNU Fortran 12.2.1 through the Intel MPI `mpif90` wrapper |
| Intel toolchain | `ifort` 2021.6.0 through `mpiifort` |
| MPI | Intel MPI 2021.12 (both toolchains) |
| BLAS/LAPACK | GNU build: reference LAPACK/BLAS 3.10.1. Intel build: MKL 2022.1 (sequential) |
| ARPACK | system ARPACK 3.8.0 (both builds) |
| Python | 3.10.11, NumPy 1.24.4, SciPy 1.10.1, pytest 6.2.5 |

The unmodified baseline did **not** compile with `make.sys.example`. The
failures, in the order they appear, were: the 100-character file-name buffers
of the code generators (any checkout path longer than that), source lines
beyond 132 columns, the `TIMER` preprocessor symbol being lost between the
nested Makefiles, `.eq.` applied to logicals, and OpenMP `REDUCTION` clauses
on `POINTER` arrays. A baseline executable therefore does not exist, and
"before/after" comparisons below are made against the reference data.

## Build

```sh
cd lanczosKuboCode
cp make.sys.example make.sys     # adjust if needed; make.sys is not tracked
make
```

`make.sys.example` worked unchanged on the cluster described above. The build
can also be started from another directory (`make -C lanczosKuboCode`) and can
be placed outside the source tree:

```sh
make -C lanczosKuboCode MAKE_SYS=/path/to/make.sys \
     BUILD_DIR=/scratch/me/grabnes/build BIN_DIR=/scratch/me/grabnes/bin
```

| Configuration | Flags | Result |
| --- | --- | --- |
| GNU, optimized (`make.sys.example`) | `-O3 -g -DMPI -DPOINTER_SIZE=8 -ffree-line-length-none -fallow-argument-mismatch`, link `-fopenmp`, `-larpack -llapack -lblas` | builds, about 70 s on one core |
| GNU, checked (`grabnes_testrun/config/gfortran.debug.make.sys`) | as above with `-O0 -g -fcheck=all -fbacktrace` | builds, about 25 s |
| Intel (`grabnes_testrun/config/intel.make.sys`) | `-O2 -g -traceback -DMPI -qopenmp -DPOINTER_SIZE=8`, `-larpack -qmkl=sequential` | builds, about 150 s, no diagnostics |
| `ifx` | | not tested |

GNU Fortran emits 35 warnings, all of legacy character: 29 argument rank/type
mismatches in the generated MPI wrappers (accepted through
`-fallow-argument-mismatch`), four uses of deleted `DO` termination forms, and
two legacy I/O list commas.

### GNU Fortran and OpenMP

GNU Fortran 12 rejects the OpenMP directives of `Src/diag.F90` because six
`REDUCTION` clauses name `POINTER` arrays (`DOS`, `Ake1Loc`, `Ake2Loc`).
`make.sys.example` therefore does **not** pass `-fopenmp` to the compiler, and
the GNU build is serial within a process; the OpenMP runtime is still linked
because the sources call its query functions. Intel Fortran accepts the
directives, and the Intel build was run with 1 and 4 threads. Making these
loops portable is a source change to the DOS and spectral-function routines
and was deliberately not done here.

## Tests

```sh
./grabnes_testrun/smoke_test.sh        # build + example 01
./grabnes_testrun/run_examples.sh      # build + examples 01-04
./grabnes_testrun/run_examples.sh --debug
./grabnes_testrun/run_examples.sh --make-sys grabnes_testrun/config/intel.make.sys
python3 -m pytest tests/tapw
```

Each example needs between 0.5 and 1.5 s on one core, so the suite can be run
on a login node; the build dominates the run time. On clusters where
executables must be started through the scheduler, set `GRABNES_LAUNCHER`
(for example `srun -n 1`).

### Results

Largest absolute deviation from the reference data (`id.` = byte-identical):

| Example | Output | GNU `-O3` | GNU `-O0` checked | Intel `-O2`, MKL, 1 and 4 threads | Tolerance |
| --- | --- | --- | --- | --- | --- |
| `01_graphene_bands` | 152 values | 0 (id.) | 0 (id.) | 0 (id.) | atol 2e-6 |
| `02_graphene_dos` | 602 values | 2.8e-16 | 2.8e-16 | 5.3e-15 | atol 1e-12, rtol 1e-9 |
| `03_twisted_bilayer_bands` | 2241 values | 0 (id.) | 0 (id.) | 0 (id.) | atol 2e-6 |
| `04_twisted_bilayer_dos` | 602 values | 1.6e-13 | 2.0e-13 | 1.8e-13 | atol 1e-12, rtol 1e-9 |

The DOS reference values reach 0.47 (example 02) and 21 (example 04), so the
largest relative deviations are 5e-15 and 9e-14. All runs reported
`0 errors, 0 warnings`. The GNU `-O3` executable gives the same results when
started with `mpirun -np 1`. The four TAPW pytest cases pass.

The checked GNU build was additionally run under Valgrind (memcheck) for all
four examples. After the fixes below, the only remaining reports are inside
MPI initialization and finalization; none comes from GRABNES source lines.

## Defects fixed in the canonical solver

Compatibility fixes taken from the former test copy:

1. `.eq.` between logicals replaced by `.eqv.` (`calc.F90`, ten places in
   `diag.F90`).
2. Unreachable block calling `move_alloc(Gx, Gx)` removed from
   `generate_G_list_from_rcell` (`diag.F90`); the loop always fills the arrays
   completely, so the block could never run.
3. `MIO/mem_temp.f90`: the size of an array is taken before it is deallocated.
4. `MIO/mem_gen.f90`, `MIO/MPI/type_gen.f90`: file-name buffers of 1024
   characters.
5. `ham.F90`: `size(rcvList)` is only evaluated when the MPI receive list is
   associated. The same unguarded expression in `calc.F90` (Kubo arrays) was
   fixed as well, and the halo-list pointers of `neigh.F90` are now
   initialized to `NULL()`; they were undefined before, which makes even
   `associated()` invalid. With one MPI process GNU Fortran returned a size of
   1 for the undefined pointer and Intel Fortran 0, so the GNU build carried a
   spurious extra halo element.

Defects that only appeared once the canonical solver was built with
optimization and run on all examples:

6. `MIO/parser.f90`: `ParserProcess` did not define its `intent(out)` error
   status on success. With `-O3` every integer input parameter failed with
   "Unrecognized error (32767)".
7. `MIO/sys.F90`: `MPI_Abort` was called without arguments. Every fatal error
   ended in a segmentation fault that also discarded the error message.
8. `ham.F90`, `HamHopping`: the local flag `helicalTwistedMBM` was tested in
   nine places but never read from the input in this routine. When it
   happened to be true, the second- and third-neighbor intralayer hoppings
   were taken from the uninitialized variables `t2KA`, `t2KB`, and `t3K`. The
   flag is now read with the other layer flags.
9. `ham.F90`, `HamInit`: the flag `cutAtNN3`, which selects the intralayer
   neighbor cutoff when `TB.NeighLevels` is larger than 2, was never assigned.
   It is now the input parameter `Neigh.CutAtNN3` (default `.false.`, the
   fifth-neighbor F2G2 cutoff).
10. `neigh.F90`: the five `fastNNnotsquare*` routines and `fastNN` wrote
    beyond the neighbor arrays when an atom had more neighbors than `maxnn`.
    They now stop with an explanatory error (`NeighCheckCount`).
11. `neigh.F90`: the same routines compared an uninitialized string with
    `'BoronNitride'`; they now read `TypeOfSystem` first, as `NeighList` does.

Build system:

12. The Makefile derived the source root from `$PWD`; it now derives it from
    its own location, and `MAKE_SYS`, `BUILD_DIR`, and `BIN_DIR` can be
    overridden.
13. `Src/Makefile` re-read `make.sys`, which discarded the `-DTIMER` that the
    top-level Makefile adds; the sources do not compile without it. The
    configuration is now read once.
14. The rule that regenerated `version.info` from `.version` was removed. It
    dates from Subversion keyword expansion, rewrote two tracked files
    depending on their time stamps, and produced a malformed version string.
    `make clean` removes only `bin/grabnes` instead of every file in the
    binary directory.

## Reference data of examples 03 and 04

Items 8 and 9 mean that the twisted-bilayer results used to depend on
uninitialized memory. On this cluster the unfixed code either failed in the
diagonalization (checked build) or crashed in the neighbor search (optimized
build).

The existing reference files were not regenerated. They are reproduced, band
file byte for byte, when the intralayer model is stated explicitly:

```text
Neigh.CutAtNN3 .true.
SingleLayert2KSL 0.0
SingleLayert3K 0.0
```

that is, by a third-shell neighbor search with the second- and third-neighbor
hoppings set to zero. This is what the September 2026 run effectively computed:
the uninitialized flag was true and the uninitialized hoppings were zero. These
three lines were added to the inputs of examples 03 and 04, so the examples
are now deterministic and their reference data unchanged.

They consequently describe **nearest-neighbor intralayer hopping with the
Koshino interlayer model**, not the F2G2 intralayer model that the solver uses
by default. With the default model (`TB.NeighLevels 5`, the three lines
removed) the bands differ by more than 2 eV; for example the Dirac point at
K moves from 0 to about -0.33 eV. Whether the examples should be switched to
the default model, with new reference data, is an open decision.

## Known limitations and unverified functionality

- **Several MPI processes.** The diagonalization examples crash when started
  with `mpirun -np 2` (in the input layer and in the neighbor search). Only
  single-process runs are validated.
- **GNU OpenMP.** See above; only the Intel build was run with threads.
- **Kubo / Lanczos transport.** A small Kubo DOS calculation completes, but
  two consecutive runs of the same executable differ by more than 10 % of the
  peak DOS: the default random-phase state is not reproducible from run to run. No
  regression test exists, and the change to `calc.F90` (item 5) could only be
  checked to be indistinguishable within that noise.
- **Other uninitialized variables.** GNU Fortran flags about 120 further
  "may be used uninitialized" locations, almost all in `ham.F90`. The compiler
  did not flag the two variables that actually broke the examples, so these
  warnings are neither complete nor all real. Only the code paths of the four
  examples were checked with Valgrind. A build with poisoned initial values
  (`-finit-*`) cannot be used as a detector yet, because the input library
  relies on zero-initialized module variables.
- **Not exercised at all:** TAPW, spin-orbit coupling, Berry curvature and
  Chern numbers, spectral functions, sparse (ARPACK) diagonalization, magnetic
  field, self-consistent Hubbard terms, semiclassical orbits, hBN and
  multilayer systems, and every neighbor routine other than
  `fastNNnotsquare` and `fastNNnotsquareBulkSmall`.
- **`ifx`, macOS, Open MPI, MPICH, OpenBLAS:** not tested with the current
  sources.

## Suggested next regression tests

Ordered by the amount of existing functionality each one protects:

1. **Analytic graphene checks**: nearest-neighbor bands against
   `±t|f(k)|` at Gamma, M, and K, and the DOS sum rule (integral equal to the
   number of bands). These do not depend on stored reference files.
2. **F2G2 twisted-bilayer example** with the default intralayer model, once
   its reference data are agreed.
3. **Symmetry invariants of the twisted-bilayer bands**: equality of the
   spectra at symmetry-related k-points and under a change of the unit-cell
   origin.
4. **Sparse against dense diagonalization** for the same small cell: every
   ARPACK eigenvalue must have a dense partner and the counts must match in
   the chosen window.
5. **TAPW against full tight binding** on a small commensurate cell, as a
   two-way comparison: each TAPW band has a tight-binding partner within a
   stated tolerance, and no tight-binding band in the window is missing.
6. **Spin-orbit baselines** from `tests/soc/README.md` (tests 0 to 3 and 9):
   spin degeneracy without SOC, a `2*lambda_I` gap, and Zeeman splitting, as
   committed inputs with numerical assertions.
7. **Chern numbers of a model with a known answer** (Haldane model, C = ±1,
   and C = 0 with the Haldane phase off) before any work on TAPW Berry
   curvature; the sum over all bands must vanish.
8. **Deterministic Kubo DOS**: after adding a reproducible seed, compare the
   Lanczos DOS with the exact-diagonalization DOS of example 02 within the
   stochastic error, and check the norm conservation of the time evolution.
9. **MPI consistency**: the same result with 1 and 2 processes, once the
   multi-process crashes are understood.

Every new feature flag should come with a run that shows the four examples
unchanged when the flag is off.
