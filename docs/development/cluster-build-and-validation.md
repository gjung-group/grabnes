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
./grabnes_testrun/run_examples.sh      # build + examples 01-04 + solver checks
./grabnes_testrun/run_examples.sh --debug
./grabnes_testrun/run_examples.sh --make-sys grabnes_testrun/config/intel.make.sys
python3 -m pytest tests/tapw
```

Each example needs between 0.5 and 1.5 s on one core and the Kubo check about
3 s and 0.5 GB, so the suite can be run on a login node; the build dominates
the run time. On clusters where
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

`run_examples.sh` also runs three solver checks, which pass with all three
builds:

| Check | What it requires | Measured |
| --- | --- | --- |
| `hamiltonian_tables` | example-03 bands rebuilt from the solver's tables within 2e-5 eV and from an independent model within 3e-4 eV (needs NumPy) | 3.9e-6 and 1.0e-4 eV |
| `kubo_graphene_dos` | recursion DOS of 180000-atom graphene, fixed seed, against the exact DOS: maximum deviation below 0.03 and rms below 0.006 | GNU 0.0082 / 0.0027, Intel 0.0115 / 0.0033 |
| `two_mpi_processes` | `mpirun -np 2` is refused with the single-process error | refused, status 1 |

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
   spurious extra halo element. With a fixed seed the Kubo DOS and recursion
   coefficients are byte-identical with and without that element.

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

Found during the follow-up investigation:

12. `ham.F90`: with `WriteDataFiles .true.` the hopping table
    (`<prefix>.s.mag`) was written to a file opened with an 80-character
    record, which stops GNU Fortran with "End of record". The file now gets a
    record long enough for all neighbors of an atom, so the solver's own
    tables can be used for verification.
13. `parallel.F90`: runs with more than one MPI process are refused with an
    error instead of overrunning arrays (see the MPI section).

Build system:

14. The Makefile derived the source root from `$PWD`; it now derives it from
    its own location, and `MAKE_SYS`, `BUILD_DIR`, and `BIN_DIR` can be
    overridden.
15. `Src/Makefile` re-read `make.sys`, which discarded the `-DTIMER` that the
    top-level Makefile adds; the sources do not compile without it. The
    configuration is now read once.
16. The rule that regenerated `version.info` from `.version` was removed. It
    dates from Subversion keyword expansion, rewrote two tracked files
    depending on their time stamps, and produced a malformed version string.
    `make clean` removes only `bin/grabnes` instead of every file in the
    binary directory.

## The Hamiltonian of examples 03 and 04

Items 8 and 9 mean that the twisted-bilayer results used to depend on
uninitialized memory. On this cluster the unfixed code either failed in the
diagonalization (checked build) or crashed in the neighbor search (optimized
build). The existing reference files were not regenerated. They are
reproduced, band file byte for byte, when the inputs state the intralayer
model explicitly, which they now do:

```text
Neigh.CutAtNN3 .true.
SingleLayert2KSL 0.0
SingleLayert3K 0.0
```

**These two examples are regression tests of one specific Hamiltonian. They
are not independent evidence that this Hamiltonian is the physically
appropriate one.** The active terms were read from the tables the solver
writes with `WriteDataFiles .true.` and checked with
`tools/hamiltonian/verify_tables.py`:

| Term | Value |
| --- | --- |
| Geometry | `(m,n) = (3,2)`, 76 atoms, twist 13.1736 degrees, cell length `a*sqrt(19)` = 10.72289 A, all bonds 1.42028 A, layers at z = 18.33 and 21.67 A (3.34 A apart) |
| Energy unit | `g0 = 12.14 - 3.72 a` = 2.98880 eV for a = 2.46 A (the default of `TB.Hopping`) |
| On-site energies | 0 |
| Intralayer, 1st shell (1.42 A, 3 neighbors) | H = -2.98880 eV |
| Intralayer, 2nd and 3rd shells (2.46 and 2.84 A) | in the neighbor list, H = 0 |
| Interlayer | two-center form `Vpi exp(-(d-a_cc)/delta)(1-(dz/d)^2) + Vsigma exp(-(d-d0)/delta)(dz/d)^2` with Vpi = -3.5 eV (`vpppi0`), Vsigma = +0.48 eV (`vppsigma0`), delta = 0.184 a = 0.4526 A (`BLdelta`), d0 = 3.34 A; the tabulated values agree with this formula to 2e-6 eV; largest element 0.48 eV |
| Interlayer range | pairs with in-plane distance below `1.1 * a_cc * Neigh.LayerDistFactor` = 9.69 A, on average 83.0 per atom |

Two remarks on the parameters: the interlayer formula uses `vpppi0` = 3.5 eV,
not the 2.7 eV of the original Moon-Koshino parametrization, and the
intralayer nearest-neighbor value is `g0`, not that same `vpppi0`. Both are
the solver's defaults and should be confirmed against the intended reference
before these examples are presented as "the Koshino model".

### Independent reconstruction

`verify_tables.py` rebuilds H(k) in Python and compares its eigenvalues with
the bands written by the solver (29 k-points, 76 bands):

| Reconstruction | Largest deviation |
| --- | --- |
| From the solver's tables, with the lattice image the solver assigns to each neighbor | 3.9e-6 eV |
| From the solver's tables, with the lattice image given by the displacement vector | 1.0e-4 eV |
| From scratch on the solver's atomic positions, complete and Hermitian neighbor set | 1.0e-4 eV |

The first line validates the diagonalization, the band output, and the
reading of the tables (its accuracy is limited by the five decimals of the
displacement files). The other two expose approximations of the neighbor
search `fastNNnotsquare` that amount to 1e-4 eV in this cell:

1. **Lattice image chosen by distance only.** After the search, the periodic
   image of each neighbor is identified by scanning lattice vectors for one
   whose *distance* matches within 0.1 A; the last match wins. In a cell of
   10.7 A with a 9.7 A interlayer range, two images of the same atom can lie at
   nearly the same distance, and 321 of 7223 entries get the wrong image and
   hence the wrong Bloch phase (matrix elements up to 1.4e-4 eV). This is the
   whole 1e-4 eV deviation.
2. **Truncated interlayer set.** The search visits a block of 5 x 5 bins of
   3 A around each atom, so pairs beyond 6 to 9 A (depending on where the atom
   sits in its bin) are missed: 6311 of the 8572 pairs inside the nominal
   cutoff are found (73.6 %). The missing elements are below 2e-5 eV and their
   effect on these bands is below 1e-6 eV.
3. **Non-Hermitian tables.** For the same reason 659 entries have no reverse
   partner. `ZHEEV` only reads the lower triangle, so the result is the
   spectrum of a Hermitian matrix, but which of the two unequal elements is
   used depends on the atom numbering.

None of this was changed: correcting item 1 moves the bands of examples 03
and 04 by up to 1e-4 eV, which is more than the 2e-6 eV tolerance of the
unchanged reference files. In large moire cells the images are unambiguous
and item 1 does not occur; items 2 and 3 are present whenever
`Neigh.LayerDistFactor` asks for more than about 6 A.

### Assessment of the default F2G2 intralayer model

With the three lines removed and `TB.NeighLevels 5` the solver uses its
default intralayer model. The tables then contain:

| Shell | Distance | Neighbors | H (eV) | Input parameter (stored as -H) |
| --- | --- | --- | --- | --- |
| 1 | 1.42 A | 3 | -2.98880 | `TB.Hopping` (`g0`) |
| 2 | 2.46 A | 6 | +0.23540 | `SingleLayert2KSL` = -0.2354 |
| 3 | 2.84 A | 3 | -0.18770 | `SingleLayert3K` = 0.1877 |
| 4 | 3.76 A | 6 | 0 | `SingleLayert4K` = 0 |
| 5 | 4.26 A | 6 | +0.06330 | `SingleLayert5KSL` = -0.0633 |
| 6, 7 | 4.92, 5.12 A | 6, 6 | 0 | inside the search radius, no term assigned |

The interlayer part is identical to the one above. Findings:

- **Implementation.** The independent model with these shells reproduces the
  solver's bands to 8.9e-5 eV, the same lattice-image deviation as before. The
  F2G2 terms are therefore assembled as tabulated.
- **Energy reference.** The Dirac point of a monolayer with these parameters
  is at `E_D = -3 t2 + 6 t5` = -3(0.2354) + 6(0.0633) = **-0.3264 eV**
  (second- and fifth-neighbor shells connect atoms of the same sublattice and
  shift both bands rigidly at K). The solver gives -0.3302 eV for the four
  Dirac states at the moire K point; the remaining 4 meV is the interlayer
  coupling. The offset is a property of the parameter set with zero on-site
  energy, not an error; the solver does not move the Dirac point to zero.
- **Dirac velocity.** `hbar v = (sqrt(3)/2) a |t1 - 2 t3|` = 5.568 eV A, that
  is v = 8.46e5 m/s, the LDA-like value this kind of fit is meant to give
  (checked numerically on the monolayer dispersion).
- **Size of the difference.** For a monolayer the band edges at Gamma move
  from -8.97/+8.97 eV (nearest neighbor) to -7.74/+11.32 eV, because the
  same-sublattice shells add `6 t2 + 6 t5` = +1.79 eV at Gamma and -0.33 eV
  at K. Differences above 2 eV between the two variants of example 03 are
  therefore expected.

What was **not** verified is the origin of the numbers themselves: the F2G2
values and the rule `g0 = 12.14 - 3.72 a` were not compared with the
publications they come from. With that caveat the default F2G2 calculation is
internally consistent and correctly implemented to 1e-4 eV. New F2G2 reference
data should be generated only after (a) the parameter values are confirmed and
(b) the lattice-image assignment is corrected, since a reference accurate to
1e-4 eV would otherwise freeze that defect in.

## `TB.NeighLevels` and `Neigh.CutAtNN3`

The intralayer part of the neighbor search is controlled in `HamInit`, the
hopping values in `HamHopping`; the two do not share information.

- The **search radius** depends on `TB.NeighLevels` only through three cases:
  1 (first shell), 2 (up to 3.10 A, which already includes the third shell),
  and "3 or more". In the last case the radius is the same for 3, 4, 5, ...
  and is chosen by `Neigh.CutAtNN3` (third shell, 3.58 A) or else by
  `F2G2Model` (5.36 A, seven shells, or 7.15 A when `F2G2Model` is false).
- The **array size** is `sum(numN(1:TB.NeighLevels)) + Neigh.LayerNeighbors`
  with `numN = [3,6,3,6,6]`. It assumes exactly `TB.NeighLevels` shells and has
  no entry beyond five.
- The **hopping terms** are assigned by distance to whatever neighbors were
  found, independently of `TB.NeighLevels`.

Measured on example 03 (`Neigh.LayerNeighbors 100`):

| `TB.NeighLevels` | `Neigh.CutAtNN3` | In-plane neighbors found / allocated | Non-zero intralayer terms | Outcome |
| --- | --- | --- | --- | --- |
| 1 | either | 3 / 3 | t1 | consistent |
| 2 | either | 12 / 9 | t1, t2 | F2G2 cut after t2; fits only thanks to the interlayer allowance |
| 3, 4, 5 | `.true.` | 12 / 12, 18, 24 | t1, t2, t3 | F2G2 cut after t3 |
| 3, 4 | `.false.` | 36 / 12, 18 | | stops: more neighbors than the arrays hold (before the fix: memory corruption) |
| 5 | `.false.` | 36 / 24 | t1, t2, t3, t5 | full F2G2; 123 of 124 list entries used |
| 6 or more | either | | | reads beyond `numN`; undefined array size |

So the two controls overlap, contradictory settings are possible, and a
truncated search radius silently truncates the Hamiltonian model. The array
size is right only for `TB.NeighLevels` 1 and for 3 with `Neigh.CutAtNN3`;
otherwise it relies on unused interlayer slots.

**Recommendation** (not implemented; it changes the input interface):

1. Keep one table of shell radii and cumulative neighbor counts, and derive
   both the search radius (midway between shell n and shell n+1) and the array
   size from it.
2. Let `TB.NeighLevels` be the only user control and mean "number of
   intralayer shells". `Neigh.CutAtNN3 .true.` is then `TB.NeighLevels 3` and
   can be removed.
3. Have each intralayer model declare the outermost shell it needs (F2G2: 5,
   the eight-neighbor set: 8) and stop with a message when `TB.NeighLevels`
   is smaller, instead of dropping terms; reject values beyond the table.
4. Until then, the provisional default `Neigh.CutAtNN3 .false.` together with
   the overflow check gives either the full model or an error for
   `TB.NeighLevels` 3 to 5. The remaining silent case is `TB.NeighLevels 2`.

## MPI

The solver cannot currently run on more than one MPI process, and this is not
specific to the examples:

- `ParallelDiv` (`parallel.F90`) computes the number of domains and then
  overrides it with `nDiv = 1`. This has been so since the first public
  commit. As a result the first process owns all atoms and every other process
  owns none (`inode1 = nAt + 1`), while the neighbor routines, the Hamiltonian
  setup, and the Kubo initialization still loop over all atoms on every
  process and write outside the empty per-process arrays. With a checked build
  the first such write is reported (`NList` in the neighbor search); an
  optimized build crashes later in an unrelated place.
- The exact-diagonalization path (`diag.F90`) contains no MPI at all; it is
  parallelized with OpenMP only.
- Removing the override in a scratch copy distributes the atoms as intended
  (3660 + 3540 for a 7200-atom cell on two processes), but the run then does
  not get past the neighbor-list generation within two minutes, where one
  process needs a few seconds. Restoring MPI domain decomposition is a
  separate task.

The solver now stops immediately with "must be run with a single MPI process"
when started on several processes. The regression suite checks exactly that
with `mpirun -np 2`. Parallel runs must use OpenMP threads, which at present
means the Intel build.

## Kubo density of states and reproducibility

Repeated Kubo runs differ because the random-phase initial state is seeded
from the system clock by default (`KuboInitWF`). This is intended stochastic
behavior, not a defect:

- `setSeed .true.` with `seedValue <integer>` fixes the seed; two runs are
  then byte-identical. These two parameters were undocumented.
- Without a fixed seed the scatter follows the expected `1/sqrt(N)` law of a
  stochastic trace: the mean relative standard deviation over four runs was
  8.5 % for 7200 atoms and 2.4 % for 115200 atoms (expected ratio 4, found
  3.5). The earlier observation of more than 10 % was a 7200-atom test.
- The result is correct: for 180000 atoms the recursion DOS agrees with the
  exact DOS of nearest-neighbor graphene, broadened by a Lorentzian of half
  width `Epsilon`, with an rms deviation of 0.0027 to 0.0033 (peak 0.417) for
  three seeds and both compilers, equal to the statistical error. Controls: a
  hopping wrong by 2 % gives 0.0093 and a broadening in the wrong energy unit
  0.012. A 20 % error in the broadening is *not* resolved (0.0033).
- The DOS file is in units of the hopping: energies in `g0`, DOS in states
  per atom and per `g0`.

The regression suite runs this comparison (`kubo_graphene_dos`). Because the
random numbers depend on the compiler, it is a statistical test with limits
of about twice the observed deviation, not a digit-by-digit comparison.

Only the DOS was examined. The time evolution, the diffusion coefficient, and
the conductivity remain untested, and the other initialization routines
(layer- and species-resolved DOS) use a different, fixed seed parameter
(`SeedSet`).

## Known limitations and unverified functionality

- **Several MPI processes.** Not supported; the solver stops with an error
  (see the MPI section).
- **GNU OpenMP.** See above; only the Intel build was run with threads, so the
  GNU build has no parallel mode at all at present.
- **Neighbor search in small cells.** Wrong lattice images and a truncated,
  non-symmetric interlayer neighbor set in `fastNNnotsquare`, worth 1e-4 eV in
  examples 03 and 04 (see above). The other `fastNN*` routines use the same
  scheme and were not examined.
- **Kubo transport.** The DOS is validated on graphene; time evolution and
  conductivity are not.
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
- **Published parameter values.** No hopping parameter was compared with its
  literature source.
- **`ifx`, macOS, Open MPI, MPICH, OpenBLAS:** not tested with the current
  sources.

## Suggested next regression tests

Ordered by the amount of existing functionality each one protects:

1. **Analytic graphene checks**: nearest-neighbor bands against
   `±t|f(k)|` at Gamma, M, and K, and the DOS sum rule (integral equal to the
   number of bands). These do not depend on stored reference files.
2. **F2G2 twisted-bilayer example** with the default intralayer model, once
   the parameter values are confirmed and the lattice-image assignment of the
   neighbor search is corrected; `tools/hamiltonian/verify_tables.py` then has
   to agree to 1e-6 eV instead of 1e-4 eV.
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
8. **Kubo time evolution**: with a fixed seed, norm conservation and the
   ballistic mean-square spreading of clean graphene, whose slope is fixed by
   the band velocity. (The Kubo DOS is covered by `kubo_graphene_dos`.)
9. **MPI consistency**: the same result with 1 and 2 processes, once domain
   decomposition is restored.

Every new feature flag should come with a run that shows the four examples
unchanged when the flag is off.
