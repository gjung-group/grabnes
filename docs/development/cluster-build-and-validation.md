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

State after the neighbor-search correction (commit series ending with the
documentation commit of 2026-10-08). Largest absolute deviation from the
reference data (`id.` = byte-identical):

| Example | Output | GNU `-O3` | GNU `-O0` checked | Intel `-O2`, MKL, 1 and 4 threads | Tolerance | Reference |
| --- | --- | --- | --- | --- | --- | --- |
| `01_graphene_bands` | 152 values | 0 (id.) | 0 (id.) | 0 (id.) | atol 2e-6 | unchanged |
| `02_graphene_dos` | 602 values | 4.7e-16 | 4.7e-16 | 5.3e-15 | atol 1e-12, rtol 1e-9 | unchanged |
| `03_twisted_bilayer_bands` | 2241 values | 0 (id.) | 0 (id.) | 0 (id.) | atol 2e-6 | **replaced**, see below |
| `04_twisted_bilayer_dos` | 602 values | 0 (id.) | 7.1e-15 | 1.9e-13 | atol 1e-12, rtol 1e-9 | **replaced**, see below |

All runs reported `0 errors, 0 warnings`. The four TAPW pytest cases pass.

`run_examples.sh` also runs checks that do not rely on stored solver output.
All pass with the three builds:

| Check | What it requires | Measured |
| --- | --- | --- |
| `dirac_point_degeneracy` | the four Dirac states of example 03 at K form two pairs degenerate within 1e-6 eV | splitting 0 at six decimals |
| `hamiltonian_example03`, `hamiltonian_f2g2`, `hamiltonian_small_cell` | solver tables against the independent model of `verify_tables.py`: no missing, extra, duplicate, or unpaired entry; matrix elements, Hermiticity, and eigenvalues within 1e-9 eV; model bands against the band file within 1e-6 eV (needs NumPy) | matrix elements 4e-16 eV (GNU), 2e-14 eV (Intel); bands 5.0e-7 eV |
| `neighlevels_0`, `neighlevels_9` | `TB.NeighLevels` outside 1 to 8 is refused with its error message | refused |
| `legacy_cutatnn3` | `TB.NeighLevels 5` + `Neigh.CutAtNN3 .true.` gives byte-identical bands to `TB.NeighLevels 3` | identical |
| `kubo_graphene_dos` | recursion DOS of 180000-atom graphene, fixed seed, against the exact DOS: maximum deviation below 0.03 and rms below 0.006 | GNU 0.0082 / 0.0027, Intel 0.0115 / 0.0033 |
| `two_mpi_processes` | `mpirun -np 2` is refused with the single-process error | refused, status 1 |

The checked GNU build was run under Valgrind (memcheck) for the four examples
and, with the new neighbor search, for the 28-atom F2G2 case. The only
reports are inside MPI initialization and finalization; none comes from
GRABNES source lines.

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
   It was first turned into the input parameter `Neigh.CutAtNN3` and is now
   deprecated in favor of `TB.NeighLevels` (see "Neighbor shells").
10. `neigh.F90`: the five `fastNNnotsquare*` routines and `fastNN` wrote
    beyond the neighbor arrays when an atom had more neighbors than `maxnn`.
    They now stop with an explanatory error (`NeighCheckCount`).
11. `neigh.F90`: the same routines compared an uninitialized string with
    `'BoronNitride'`; they now read `TypeOfSystem` first, as `NeighList` does.

Found during the follow-up investigations:

12. `ham.F90`: with `WriteDataFiles .true.` the hopping table
    (`<prefix>.s.mag`) was written to a file opened with an 80-character
    record, which stops GNU Fortran with "End of record". The file now gets a
    record long enough for all neighbors of an atom, so the solver's own
    tables can be used for verification.
13. `parallel.F90`: runs with more than one MPI process are refused with an
    error instead of overrunning arrays (see the MPI section).
14. `neigh.F90`: wrong lattice translations, incomplete and one-sided neighbor
    lists (see "Neighbor search").
15. `ham.F90`, `tbpar.f90`: inconsistent shell controls and array sizes (see
    "Neighbor shells").
16. `ham.F90`: seventh-shell hopping assigned the sixth-shell value in the
    bilayer intralayer branch.

Build system:

17. The Makefile derived the source root from `$PWD`; it now derives it from
    its own location, and `MAKE_SYS`, `BUILD_DIR`, and `BIN_DIR` can be
    overridden.
18. `Src/Makefile` re-read `make.sys`, which discarded the `-DTIMER` that the
    top-level Makefile adds; the sources do not compile without it. The
    configuration is now read once.
19. The rule that regenerated `version.info` from `.version` was removed. It
    dates from Subversion keyword expansion, rewrote two tracked files
    depending on their time stamps, and produced a malformed version string.
    `make clean` removes only `bin/grabnes` instead of every file in the
    binary directory.

## Neighbor search

### Defects of the original search

The default search (`fastNNnotsquare`, and its `Bulk`/`BulkSmall` variants)
binned the atoms and their eight adjacent periodic images and then, in a
second step, guessed the lattice translation of every neighbor. Measured on
the 76-atom twisted bilayer of examples 03 and 04:

| Defect | Evidence | Effect on the bands |
| --- | --- | --- |
| Lattice translation inferred from a *distance* match within 0.1 A (last match in a scan wins) | 321 of 7223 entries carried the translation of another image and hence a wrong Bloch phase | 1.03e-4 eV |
| Only 5 x 5 bins of 3 A visited; only adjacent images generated | 6311 of the 8572 pairs inside the interlayer radius found (73.6 %); the missing matrix elements reach 2e-5 eV | 1.3e-5 eV |
| List not symmetric (consequence of the previous item) | 659 entries without the reverse entry; `ZHEEV` reads one triangle, so which of two unequal elements was used depended on the atom numbering | 2.1e-6 eV |
| Fixed bin capacity and guessed array size (`sum of shell counts + Neigh.LayerNeighbors`) | overflow, see the earlier defect list | |

The numbers in the last column are the largest eigenvalue changes along the
band path when the defects are removed one after the other in a Python
reconstruction of the historical tables.

### Corrected algorithm

`NeighSearchLayered` (`neigh.F90`) is now the single search behind the three
routines:

- For every atom `i` and every periodic image `m + u1 A1 + u2 A2 (+ u3 A3)`
  it forms `d = (r_m - r_i) + T`. A pair is kept when `|dz| < 1.5 A` and
  `dx^2 + dy^2` is below the squared intralayer radius, or when
  `1.5 A < |dz| < 4.5 A` and it is below the squared interlayer radius
  (`4.5 A < |dz| < 7.5 A` as well with `addSecondLayerInteractions`). These
  are the criteria of the original routines.
- The stored translation `neighCell = (u1, u2, u3)` is the one that built the
  image. Nothing is inferred afterwards.
- As many images are generated as the radii require, from the distance
  between lattice lines and the spread of the fractional coordinates, so
  cells smaller than the radii are handled (28-atom cell of 6.5 A with a
  9.7 A radius tested).
- Bins are at least one radius wide and are linked lists without capacity.
- Because `d` is evaluated as `(r_m - r_i) + T`, the entries `(i -> m, T)` and
  `(m -> i, -T)` see exactly opposite vectors and are accepted or rejected
  together: the list is symmetric by construction.
- The neighbor arrays are allocated after counting. `Neigh.LayerNeighbors`
  therefore only switches the interlayer search on (non-zero) or off.
- The search is faster than before (0.6 s instead of 1.8 s for 180000 atoms).

`fastNNnotsquareSmall` (unreachable), `fastNNnotsquareNotRectangle`, `fastNN`,
and the MPI-era `NeighList` were not changed and have not been examined.

### Bloch convention

`DiagHam` builds, with `R = ucell * neighCell(:,j,i)` and `m = NList(j,i)`,

```text
H(m,i)(k) = sum_j  -hopp(j,i) * exp(-i k.R) * g0
```

Only lattice vectors enter the phases. Consequences, all tested numerically:

- Hermiticity requires for every entry `(i, m, R)` the entry `(m, i, -R)` with
  the conjugate value. This now holds for the stored tables themselves, not
  only after `ZHEEV` has ignored one triangle.
- `H(k + G) = H(k)` element by element for reciprocal lattice vectors `G`.
- Moving atom `a` by a lattice vector `T_a` changes the translations, not the
  displacements: `H'(k) = D^+ H(k) D` with `D = diag(exp(-i k.T_a))`, a unitary
  change of basis. Matrix elements differ; eigenvalues do not.
- The convention with the full displacement in the phase, `exp(-i k.d)`, has
  the same eigenvalues but is not periodic in `k` element by element.

### Independent validation

`tools/hamiltonian/verify_tables.py` builds a second Hamiltonian from the
atomic positions and the model parameters alone: it enumerates every pair
`(i, m, R)` inside the search radii exhaustively and evaluates the matrix
element from the formulas. The solver's neighbor list, translations, and
hopping values are not used for it. With `WriteDataFiles .true.` the solver
writes its own tables, including the translation and the full-precision
displacement of every entry (`neighCell.dat`, `neighD.dat`), and the two are
compared entry by entry in real space, which is independent of any Bloch
gauge. Results with the checked GNU build (the optimized GNU and the Intel
builds give the same within 2e-14 eV):

| Case | Atoms | Entries | Missing / extra / duplicate / unpaired | max abs(H_solver - H_model), entries | Hermiticity of H(k) | max abs(dE) at 9 k-points | Model bands vs. band file |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Example 03, nearest-neighbor intralayer | 76 | 8800 | 0 / 0 / 0 / 0 | 4.4e-16 eV | 3.1e-17 eV | 2.7e-14 eV | 5.0e-7 eV |
| Same cell, F2G2 intralayer (5 shells) | 76 | 10396 | 0 / 0 / 0 / 0 | 4.4e-16 eV | 3.1e-17 eV | 2.0e-14 eV | 5.0e-7 eV |
| `(m,n) = (2,1)`, 21.79 degrees, F2G2 | 28 | 3796 | 0 / 0 / 0 / 0 | 4.4e-16 eV | 1.1e-16 eV | 1.4e-14 eV | 5.0e-7 eV |
| `(m,n) = (5,4)`, 7.34 degrees, nearest neighbor | 244 | 28192 | 0 / 0 / 0 / 0 | 8.3e-16 eV | below 1e-17 eV | 3.6e-14 eV | 5.0e-7 eV |
| Same cell as example 03, eight-shell model (`F2G2Model .false.`) | 76 | 11536 | 0 / 0 / 0 / 0 | 4.4e-16 eV | 3.1e-17 eV | 2.8e-14 eV | 5.0e-7 eV |
| Pristine graphene (example 01, `nonBulkSmall`) | 2 | 6 | 0 / 0 / 0 / 0 | 0 | 8.9e-16 eV | 5.6e-17 eV | 5.0e-7 eV |

The nine k-points are Gamma, K, M, two fixed generic points, and four points
from a seeded random generator. The last column is limited by the six
decimals of the band file (5e-7 eV is exactly half a unit of the last digit);
it was also evaluated on a solver band path through five generic k-points,
with the same result. In addition, for every case: stored displacements equal
`r_m + R - r_i` to 4e-15 A, `H(k+G) - H(k)` is below 2e-14 eV, moving atoms by
random lattice vectors leaves the number of pairs and the eigenvalues
unchanged (2e-14 eV) and reproduces `D^+ H D` to 1e-14 eV, and the
displacement-phase convention gives the same eigenvalues (2e-14 eV).

A check that needs no second implementation: the four Dirac states of example
03 at the moire K point must form two degenerate pairs. The historical data
had -0.000514, -0.000419, -0.000291, -0.000197 eV; the corrected solver gives
-0.000411 twice and -0.000299 twice.

The independent model shares with the solver the atomic positions (whose
honeycomb geometry, twist angle, and cell were checked separately), the
pair-selection criteria, and the parameter values. Agreement therefore
establishes that the solver assembles the Hamiltonian those criteria and
parameters define; it says nothing about whether they are the right physics.

## Neighbor shells: `TB.NeighLevels`

`TB.NeighLevels` is the number of intralayer neighbor shells, from 1 to 8. One
table in the `neigh` module holds the shell radii of the honeycomb lattice
(1, sqrt(3), 2, sqrt(7), 3, sqrt(12), sqrt(13), 4 bond lengths, with 3, 6, 3,
6, 6, 6, 6, 3 sites), and `NeighShellCutoff2` places the search radius halfway
between the last requested shell and the next one.

| Quantity | Determined by |
| --- | --- |
| Intralayer search radius | `TB.NeighLevels` through the shell table |
| Interlayer search radius | `1.1 a_cc * Neigh.LayerDistFactor` if `Neigh.LayerNeighbors` is non-zero; 1 A otherwise (vertical pairs only) |
| Array sizes | the search itself (largest neighbor count found) |
| Hopping value of a pair | its distance, in `HamHopping`, as before |

Rules and backward compatibility:

- Values below 1 or above 8 are refused with an error (previously: array
  overruns).
- `Neigh.CutAtNN3 .true.` is deprecated. It still limits the search to three
  shells, exactly as before, and prints a notice; `TB.NeighLevels 5` with
  `Neigh.CutAtNN3 .true.` gives byte-identical bands to `TB.NeighLevels 3`.
- If fewer shells are searched than the intralayer model reaches (5 for F2G2,
  8 with `F2G2Model .false.`), the solver prints which terms are left out.
  Truncating a model this way remains allowed, because that is what
  `TB.NeighLevels` 1 to 3 always meant; it is no longer silent.
- Inputs that ran before give the same Hamiltonian, apart from the
  neighbor-search corrections: the earlier radii for `TB.NeighLevels` 2 and 5
  also collected the third, respectively the sixth and seventh, shell with zero
  hopping, and those zero entries are gone. `TB.NeighLevels` 3 or 4 without
  `Neigh.CutAtNN3`, which used to overflow, now simply means 3 or 4 shells.
- Every value from 1 to 8 was verified against the independent model for the
  F2G2 parameters, and 5 and 8 for the eight-shell parameters.

Limitation: the shell radii are taken from the graphene lattice parameter.
The halfway rule tolerates a few percent of strain, except between the sixth
and seventh shells, which are only 4 % apart; `TB.NeighLevels 6` is marginal
for lattices that deviate from it by 2 % (hBN).

The eight-shell check exposed one more defect: in the bilayer intralayer
branch of `HamHopping` the sixth-shell distance window (4.43 to 5.41 A) also
contained the seventh shell at 5.12 A, which received `t6` instead of `t7`.
The window is now 4.77 to 5.07 A as in the other branches. The default F2G2
model has no sixth- or seventh-shell term and was not affected.

## Examples 03 and 04: Hamiltonian, parameters, reference data

**These examples are regression tests of one specific Hamiltonian. Agreement
with the independent model shows that it is assembled as specified, not that
it is the physically appropriate model.**

| Term | Value |
| --- | --- |
| Geometry | `(m,n) = (3,2)`, 76 atoms, twist 13.1736 degrees, cell length `a*sqrt(19)` = 10.72289 A, bonds 1.42028 A, layers 3.34 A apart |
| Energy unit | `g0 = 12.14 - 3.72 a` = 2.98880 eV for a = 2.46 A (default of `TB.Hopping`) |
| On-site energies | 0 |
| Intralayer (`TB.NeighLevels 1`) | first shell only, H = -2.98880 eV |
| Interlayer | `Vpi exp(-(d-a_cc)/delta)(1-(dz/d)^2) + Vsigma exp(-(d-d0)/delta)(dz/d)^2` with `vpppi0` = 3.5 eV (Vpi = -3.5 eV), `vppsigma0` = 0.48 eV, `BLdelta` = 0.184 a = 0.4526 A, d0 = 3.34 A |
| Interlayer range | in-plane distance below 9.686 A: 112.8 neighbors per atom, elements from 0.48 eV down to 1.4e-9 eV |

The inputs now contain `TB.NeighLevels 1` instead of the three lines used
before (`Neigh.CutAtNN3 .true.`, `SingleLayert2KSL 0.0`, `SingleLayert3K 0.0`);
the Hamiltonian is the same.

### Parameter provenance

| Parameter | Status |
| --- | --- |
| Two-center form, `vppsigma0` = 0.48 eV, `BLdelta` = 0.184 a | Moon and Koshino, Phys. Rev. B 85, 195458 (2012), checked against the preprint (arXiv:1202.4365): "V0ppπ ≈ −2.7 eV, V0ppσ ≈ 0.48 eV", decay length "0.184a" |
| `vpppi0` = 3.5 eV | **Deliberate phenomenological calibration**, not the Moon-Koshino value of 2.7 eV. It was chosen to obtain a more realistic Dirac velocity; see the next section for where it acts. Do not change it to 2.7 eV |
| Interlayer distance 3.34 A | set in the example inputs; Moon and Koshino quote about 3.35 A; the solver default is 3.22 A |
| `g0 = 12.14 - 3.72 a` (2.9888 eV) | origin not identified; not found in the two Jung-MacDonald papers consulted. Unverified |
| F2G2 single-layer values (`SingleLayert2KSL` -0.2354, `SingleLayert3K` 0.1877, `SingleLayert5KSL` -0.0633) | not found in the text of Jung and MacDonald, Phys. Rev. B 87, 195450 (2013), whose five-parameter set starts from t1 = -3.00236 eV. Unverified. They do reproduce the Dirac-velocity coefficient quoted for the LDA bands in Jung and MacDonald, arXiv:1309.5429 (5.567 eV A, 8.45e5 m/s): the solver gives 5.568 eV A |
| Eight-shell values used with `F2G2Model .false.` | not found in the same text. Unverified |

No parameter value was changed.

### Dirac velocity and the role of `vpppi0`

In Moon and Koshino's model one parameter, V0ppπ, sets both the intralayer
hopping and the pi part of the interlayer hopping, and the monolayer velocity
is `sqrt(3) a |V0ppπ| / (2 hbar)`. In GRABNES `vpppi0` enters the intralayer
hopping only with `KoshinoIntralayer .true.`; by default, and in the public
examples, the intralayer nearest-neighbor hopping is `g0`. All velocities
below are `hbar v = dE/dk` in eV A (and v in 1e5 m/s).

Bare monolayer values (analytic or from the monolayer dispersion):

| Intralayer model | hbar v | v |
| --- | --- | --- |
| Nearest neighbor, t = `g0` = 2.9888 eV | 6.367 | 9.67 |
| Nearest neighbor, t = 2.7 eV | 5.752 | 8.74 |
| Nearest neighbor, t = 3.5 eV | 7.457 | 11.33 |
| Two-center form on five shells, `vpppi0` = 2.7 eV | 5.220 | 7.93 |
| Two-center form on five shells, `vpppi0` = 3.5 eV | 6.767 | 10.28 |
| F2G2, `(sqrt(3)/2) a abs(t1 - 2 t3)` | 5.568 | 8.46 |

Twisted bilayer at 13.17 degrees, from solver bands on a line through the
moire K point: half the difference between the mean of the two upper and the
mean of the two lower Dirac bands at `+q` and `-q`, divided by `q`, for
`q` = 0.0039 and 0.0098 1/A and extrapolated quadratically to `q = 0`:

| Solver input | Dirac energy (eV) | hbar v | v | Relative to the bare monolayer |
| --- | --- | --- | --- | --- |
| Example 03 (`TB.NeighLevels 1`, `vpppi0` 3.5) | -0.00036 | 6.2327 | 9.47 | -2.1 % |
| Example 03 with `vpppi0 2.7` | -0.00035 | 6.2339 | 9.47 | -2.1 % |
| `KoshinoIntralayer .true.`, 5 shells, `vpppi0` 3.5 | +1.01282 | 6.6468 | 10.10 | -1.8 % |
| `KoshinoIntralayer .true.`, 5 shells, `vpppi0 2.7` | +0.77975 | 5.0649 | 7.70 | -3.0 % |
| F2G2 intralayer (`TB.NeighLevels 5`), `vpppi0` 3.5 | -0.33023 | 5.4274 | 8.25 | -2.5 % |

Reading:

- In the example configuration the value of `vpppi0` is irrelevant for the
  velocity (0.02 %), because it only scales the pi part of the interlayer
  hopping, which is small at interlayer distances. The velocity there is set
  by `TB.Hopping`.
- With `KoshinoIntralayer .true.`, 3.5 eV instead of 2.7 eV raises the
  velocity by 31 %, to 1.01e6 m/s in this bilayer. This is the configuration
  in which the calibration does what it is meant to do.
- In all cases the bilayer velocity is 2 to 3 % below the bare monolayer
  value, the expected size of the interlayer renormalization at this angle.
- No experimental target value is asserted here; these numbers only state
  what the implemented models give.

With the interlayer coupling switched off (`deactivateInterlayer .true.`) the
solver reproduces the monolayer expressions directly: Dirac energy -0.326400
eV for F2G2 (analytic -0.326400) and 0 for nearest neighbor, velocities 5.5683
and 6.3682 eV A (analytic 5.5676 and 6.3674; the difference is the
extrapolation error).

### The F2G2 intralayer model

With `TB.NeighLevels 5` the matrix elements are -2.98880 eV (shell 1),
+0.23540 (2), -0.18770 (3), 0 (4), and +0.06330 eV (5); as input parameters
they are entered with the opposite sign. The solver agrees with the
independent model to 4e-16 eV (table above).

- The Dirac point of the decoupled layers lies at
  `-3 H2 + 6 H5 = 3*SingleLayert2KSL - 6*SingleLayert5KSL` = -0.3264 eV: the
  second and fifth shells connect sites of the same sublattice (six sites
  each) and their Bloch sums at K are -3 and +6. This is an energy offset of
  the parameter set with zero on-site energy, confirmed by the decoupled run
  above; the interlayer coupling adds -3.8 meV. The solver does not recenter
  energies.
- Its bands differ from the nearest-neighbor model by more than 2 eV at the
  band edges, as they must: at Gamma the same shells add
  `6 H2 + 6 H5` = +1.79 eV.

No F2G2 reference data were added, because the provenance of the parameter
values is not established. The model is covered by the independent-model
check, which needs no stored reference.

### Reference data replaced

The reference files of examples 03 and 04 were regenerated with the corrected
neighbor search after the validation above. The historical files remain in
the Git history (last present in commit `b360f13`).

| | Historical | Corrected | Largest change | rms change |
| --- | --- | --- | --- | --- |
| `03/reference/bands.dat` (md5) | `7d74a17f...` | `0a7c9313...` | 1.03e-4 eV (Dirac states at K) | 1.7e-5 eV |
| `04/reference/dos.dat` (md5) | `bd2d0d7a...` | `70fab82d...` | 1.44e-4 of 15.7 (relative 9e-6) | 2.2e-5 |

1712 of the 2204 band energies change in the sixth decimal or more; the mean
change is 9e-9 eV and the integrated DOS changes by 2e-7 relative. Causes, in
order of size: lattice translations (1.03e-4 eV), completed neighbor list
(1.3e-5 eV), symmetric list (2.1e-6 eV). Shell selection contributes nothing:
the Hamiltonian terms are the same. The corrected files are byte-identical
between the optimized GNU, checked GNU, and Intel builds for the bands, and
agree to 2e-13 for the DOS. Examples 01 and 02 are unaffected (byte-identical
bands; DOS within 5e-16).

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
- **Other neighbor routines.** `fastNNnotsquareNotRectangle`, `fastNN`, and
  `NeighList` were not corrected or examined; `fastNNnotsquareSmall` is
  unreachable. The `ReadDataFiles` path only reconstructs translations for
  first-shell neighbors unless `readNeighborDetails` is set.
- **Hopping models other than those listed.** The independent check covers
  the two-center interlayer form and the F2G2 and eight-shell intralayer
  values of a twisted bilayer, and nearest-neighbor graphene. hBN,
  multilayers, strain, relaxation, magnetic field, and the many other branches
  of `HamHopping` are not covered.
- **`KoshinoIntralayer`.** Used for the velocity table; its Dirac energy and
  velocity are consistent with the two-center form, but it was not compared
  entry by entry with the independent model.
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
- **Parameter provenance.** Only the Moon-Koshino interlayer form and values
  were confirmed against a source; see the provenance table.
- **`ifx`, macOS, Open MPI, MPICH, OpenBLAS:** not tested with the current
  sources.

## Suggested next regression tests

Ordered by the amount of existing functionality each one protects:

1. **Analytic graphene checks**: nearest-neighbor bands against
   `±t|f(k)|` at Gamma, M, and K, and the DOS sum rule (integral equal to the
   number of bands). These do not depend on stored reference files.
2. **F2G2 twisted-bilayer example** with stored reference data, once the
   provenance of the parameter values is confirmed. The implementation itself
   is already covered by `hamiltonian_f2g2`.
3. **Further symmetry invariants of the twisted-bilayer bands**: equality of
   the spectra at symmetry-related k-points. (The Dirac-point degeneracy and
   the invariance under unit-cell images are covered.)
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
