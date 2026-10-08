# GRABNES pre-announcement collaboration log

This document records repository changes made with Codex assistance while the
public GRABNES repository is being prepared for a broader announcement. The
repository is already publicly accessible, but the work described here should
be treated as pre-release consolidation until collaborators have reviewed and
committed it.

Last updated: 2026-10-08

## Solver reconciliation, Linux validation, and test harness (2026-10-08)

Work on branch `nicolas/development`, starting from commit `941e859`. Full
details, measured results, and limitations are in
[`cluster-build-and-validation.md`](cluster-build-and-validation.md).

- **One solver tree.** The compatibility fixes that existed only in
  `grabnes_testrun/lanczosKuboCode/` were reviewed individually and applied to
  `lanczosKuboCode/` (logical `.eqv.`, removal of the dead `move_alloc` block,
  the MPI receive-list guard, memory accounting, generator file-name buffers).
  The duplicated solver tree and the duplicated graphene example were then
  removed from `grabnes_testrun/`.
- **Further defects fixed.** Building the canonical tree with optimization and
  running all four examples exposed defects the laptop test had not reached:
  an undefined `intent(out)` status in the input parser, `MPI_Abort` called
  without arguments (fatal errors ended in a segmentation fault without their
  message), two flags used without ever being assigned (`helicalTwistedMBM` in
  `HamHopping`, `cutAtNN3` in `HamInit`), unchecked writes beyond the neighbor
  arrays, an uninitialized string in the neighbor routines, and an unguarded
  `size()` of the receive list in the Kubo setup.
- **New input parameter** `Neigh.CutAtNN3` (default `.false.`) replaces the
  never-assigned `cutAtNN3`.
- **Build system.** The Makefile no longer depends on `$PWD`, supports
  `MAKE_SYS`, `BUILD_DIR`, and `BIN_DIR` for out-of-tree builds, passes
  `-DTIMER` to all sub-builds as intended, and no longer rewrites
  `.version`/`version.info`. `make.sys.example` now contains the flags GNU
  Fortran needs and does not enable OpenMP directives, which GNU Fortran
  rejects in `Src/diag.F90`.
- **Examples 03 and 04.** Their results used to depend on uninitialized
  memory. The reference files are unchanged; the inputs now state explicitly
  the model those files correspond to (nearest-neighbor intralayer hopping).
  Switching the examples to the default F2G2 intralayer model would require
  new reference data and is left as an open decision.
- **Harness.** `grabnes_testrun/` is now a small harness that builds the
  canonical sources out of tree and compares results numerically:
  `smoke_test.sh` (example 01) and `run_examples.sh` (all four), with checked
  GNU and Intel configurations. The example launchers no longer look for a
  test-copy executable.
- **Validation on Linux x86_64.** All four examples reproduce their reference
  data with GNU Fortran 12.2.1 (`-O3` and `-O0 -fcheck=all`) and with Intel
  `ifort` 2021.6 + MKL (1 and 4 OpenMP threads): band files byte-identical,
  DOS files within 2e-13. The four TAPW pytest cases pass.
- **Not validated:** OpenMP with GNU Fortran, Kubo time evolution and
  conductivity, and all TAPW, SOC, Berry-curvature, and semiclassical
  functionality. Runs with more than one MPI process are not supported.

### Follow-up investigation (same day)

- **Examples 03 and 04** are now documented as regression tests of one
  specific Hamiltonian, with its active terms listed. A new tool,
  `tools/hamiltonian/verify_tables.py`, rebuilds that Hamiltonian
  independently; it agreed with the solver to 1e-4 eV and traced the
  remainder to the neighbor search (corrected in the next entry).
- **Default F2G2 model** assessed: correctly assembled to 1e-4 eV; its Dirac
  point at -0.33 eV follows from the parameters (`-3 t2 + 6 t5`). No F2G2
  reference data were added; the parameter values still need to be checked
  against their source.
- **`TB.NeighLevels` / `Neigh.CutAtNN3`** analyzed; the controls overlap and
  can contradict each other. A single-control design is recommended in the
  validation document; the interface was not changed.
- **MPI.** Multi-process runs never worked in this source: the domain
  decomposition is switched off in `ParallelDiv`. The solver now refuses to
  start on more than one process instead of corrupting memory.
- **Kubo DOS.** Run-to-run differences are the expected noise of a
  clock-seeded random-phase state; `setSeed`/`seedValue` make runs
  reproducible. The recursion DOS of graphene agrees with the exact result
  within its statistical error.
- **`WriteDataFiles .true.`** no longer stops at the hopping table.
- **Harness.** `run_examples.sh` gained three checks: `hamiltonian_tables`,
  `kubo_graphene_dos`, and `two_mpi_processes`.

### Neighbor-search correction and Hamiltonian validation (same day)

- **Neighbor search rewritten.** The default search and its two bulk variants
  now share one exhaustive routine. It stores the lattice translation that
  built each periodic image instead of guessing it from a distance, generates
  as many images as the search radii need, produces a symmetric list, and
  sizes the arrays from the result. In the 76-atom twisted bilayer the old
  search attached wrong Bloch phases to 321 of 7223 entries (1.0e-4 eV in the
  bands), missed 26 % of the interlayer pairs (1.3e-5 eV), and left 659
  entries without reverse partner (2e-6 eV).
- **`TB.NeighLevels`** is now the single control for the number of intralayer
  shells (1 to 8). `Neigh.CutAtNN3` is deprecated but keeps its meaning;
  values outside the range are refused; `Neigh.LayerNeighbors` only switches
  the interlayer search on.
- **Seventh-shell hopping.** In the bilayer intralayer branch the seventh
  shell received the sixth-shell value (non-default eight-shell models only).
- **Independent validation.** `tools/hamiltonian/verify_tables.py` now builds
  a Hamiltonian from positions and parameters alone and compares it with the
  solver's tables entry by entry: agreement to 4e-16 eV for six cases (three
  twisted cells, F2G2 and eight-shell models, graphene), with Hermiticity,
  reciprocal-lattice periodicity, and invariance under unit-cell images
  checked.
- **Reference data of examples 03 and 04 replaced** after that validation;
  the historical files remain in the Git history (commit `b360f13`). Largest
  change 1.03e-4 eV. Their inputs now use `TB.NeighLevels 1`.
- **Parameters.** `vpppi0` = 3.5 eV is documented as a deliberate calibration
  relative to the Moon-Koshino value of 2.7 eV, with the measured Dirac
  velocities for both; the Moon-Koshino form and remaining values were checked
  against the preprint. The F2G2 values and the rule for `g0` could not be
  traced to a source and are marked unverified. No value was changed and no
  F2G2 reference data were added.
- **Harness.** New checks: Dirac-point degeneracy, three independent
  Hamiltonian comparisons, shell-control errors, and the legacy
  `Neigh.CutAtNN3` equivalence.

### Parameter convention, wider validation, and release readiness (same day)

- **`vpppi0` is model dependent.** It is one variable that sets the intralayer
  hopping with `KoshinoIntralayer` and the pi part of the interlayer hopping in
  every model. Its default was 3.5 eV for all models; it is now 2.7 eV (Moon
  and Koshino) with `KoshinoIntralayer .true.` and 3.5 eV otherwise, and the
  value is printed. Both combinations are in the regression suite. Measured:
  with the Koshino intralayer model 3.5 eV instead of 2.7 eV raises the Dirac
  velocity by 31 to 48 % (13.2 to 3.9 degrees); with the F2G2-type models the
  choice changes the velocity by 0.03 to 0.3 % only, so the velocity
  calibration attributed to 3.5 eV is not reproduced by this code and needs
  clarification by the authors.
- **More undefined behavior removed:** aliased defaults in `TB.Hopping` and
  `TB.BNHopping` (zero B-N hopping with GNU `-O3`), hoppings converted through
  single precision, and an unset factor in the Koshino intralayer hopping
  (wrong and non-Hermitian matrix with Intel Fortran).
- **Legacy neighbor routines audited.** `NeighList` agrees with the default
  search on graphene supercells; `Neigh.fastNN` and
  `Neigh.fastNNnotsquareNotRectangle` overrun their arrays and are now refused.
- **Other systems.** hBN monolayer validated against an analytic model;
  twisted bulk, trilayer-cell, `GBNtwoLayers`, and encapsulated cells pass
  structural checks only; `Graphene_Over_BN` is inconsistent with the inputs
  tried and stays unvalidated.
- **Performance.** Setup (neighbor search and hopping assembly) scales
  linearly up to 11.5 million atoms and 69 million entries; the search is 6 to
  15 % of it.
- **Release readiness.** New page `functionality-status.md`; a GitHub Actions
  workflow added but not yet run; no license file (and conflicting license
  statements) and no citation file remain release blockers.

### Release preparation (2026-10-09)

- **Scope.** Physics development is frozen. The four examples are the
  acceptance criterion for the first release: from a clean clone the
  documented build succeeds, all four run, reproduce their reference data
  (bands byte-identical, DOS within 5e-16), and their plotting scripts work.
- **Effective graphene/hBN model.** Its on-site and bond terms agree with the
  published expressions, but its stored hopping table is not Hermitian. The
  solver now reports the asymmetry; `Hopping.Symmetrize` is an optional
  remedy, off by default. No further change is planned here: the correction
  is maintained outside this repository.
- **`GBNtwoLayers`** is left untouched; the single statement of this series
  inside its branch was restored to the original. `Graphene_Over_BN` is
  documented as legacy.
- **Repository cleanup.** 45 tracked macOS executables, 38 duplicated
  notebook checkpoints, and a `.DS_Store` removed; `.gitignore` extended.
- **Documents added:** new root README, `release-readiness.md`,
  `licensing-audit.md`, `citation-checklist.md`, `software-paper-outline.md`,
  and an index of the development documentation.
- **License and citation (same day).** GRABNES is distributed under
  `GPL-3.0-or-later`: `LICENSE`, `COPYING.LESSER`, and
  `THIRD_PARTY_LICENSES.md` added, `pyproject.toml` made consistent, the
  commented-out Numerical Recipes routines removed from `kubo.F90`.
  `CITATION.cff` names the authors Jeil Jung, Rafael Martinez-Gordillo, and
  Nicolas Leconte; Nicolas Leconte is the current developer and maintainer.
- **CI.** The workflow now runs on pushes to the development and main
  branches and on pull requests to main.

Sections 1, 5, and 6 below describe the state before this work; where they
mention the test copy or a `grabnes_testrun` executable, this section
supersedes them.

## Documentation and test organization (working tree)

- Expanded the root README into the public quick-start and repository map.
- Corrected the Sphinx and Doxygen project identity to GRABNES.
- Moved user, theory, development, test, and debugging material out of the
  repository root into `docs/`, `tests/`, and `tools/`.
- Replaced overlapping historical SOC/TAPW plans with one source-reviewed
  implementation-status page.
- Converted the standalone TAPW mathematical checks into deterministic pytest
  tests and consolidated shared Brillouin-zone parsing code.

Validation on this laptop: all four TAPW pytest cases pass; both TAPW plotting
commands generated PNG files from a synthetic debug fixture; Python syntax and
the Git whitespace check pass. A warnings-as-errors Sphinx build of all guide
content also passes. The full API build still requires Doxygen, which is not
installed on this laptop but is installed by the Read the Docs configuration.

These changes are included in the documentation-cleanup update for
collaborator review.

## Current Git state

- The consolidation and initial examples were committed and pushed to `main`
  as commit `94964340c3322751f352e44337715912ca13505f`.
- The follow-up public update adds the `grabnes_testrun/` harness described
  below for collaborator testing in commits `e57b597` and `47cef60`.
- Commit `68dd8b9` completes the four-example suite described below.
- The tag `pre-code-consolidation` points to commit `6133ab7`, the repository
  state before the solver directories were consolidated.
- Deleted files remain recoverable from that tag and from Git history.

## 1. Public examples scaffold

Added an `examples/` entry point for new users:

```text
examples/
├── README.md
├── 01_graphene_bands/
├── 02_graphene_dos/
├── 03_twisted_bilayer_bands/
└── 04_twisted_bilayer_dos/
```

Each directory now contains a minimal `Gendata.in`, location-independent
`run.sh`, Python/Matplotlib `plot.py`, focused README, and known-good reference
data. The launchers automatically find either the canonical executable or the
temporary verified `grabnes_testrun` executable and still accept an explicit
`GRABNES_BIN` override.

The examples are:

- pristine graphene bands along `K - Gamma - M - K'`;
- pristine graphene DOS from exact diagonalization on a `30 x 30` k-grid;
- bands of the 76-atom `(m,n)=(3,2)` commensurate twisted bilayer; and
- DOS of the same twisted bilayer on an `8 x 8` k-grid.

The small 13.17-degree twisted cell is intentional: it demonstrates the full
moire-cell workflow while remaining a laptop-scale test. The DOS examples use
deterministic exact diagonalization rather than the much larger stochastic
Kubo calculations from the 2023 tutorial. All four calculations completed with
`0 errors, 0 warnings`, and a second end-to-end run reproduced every reference
dataset numerically.

## 2. Solver directory consolidation

The repository previously contained three competing solver trees:

- `lanczosKuboCode/`
- `lanczosKuboCode_jiaqi/`
- `lanczosKuboCode_jinwoo/`

Investigation of the private repository history showed that the `jinwoo` tree
contains the active development line through September 2026, including recent
TAPW, Berry-link, Wilson-loop, orbital-moment, and semiclassical-orbit work.
The `jiaqi` source tree was also present byte-for-byte inside
`lanczosKuboCode_jinwoo/SrcJiaqi/`.

The latest `jinwoo` implementation was therefore promoted to the canonical
`lanczosKuboCode/` path. The `_jiaqi` and `_jinwoo` directories were removed.
The embedded duplicate `SrcJiaqi/` tree was also removed.

All repository documentation references were changed from
`lanczosKuboCode_jinwoo` to `lanczosKuboCode`, including:

- `Doxyfile`
- `DOCUMENTATION_SETUP.md`
- `SOC_IMPLEMENTATION_SUMMARY.md`
- `TAPW_CHERN_OPTIMIZATION_GUIDE.md`

The generated Chern-file header in `Src/diag.F90` now identifies the producer
as `GRABNES` instead of referring to the former personal directory name.

## 3. Repository cleanup

Removed material that should not be maintained in the public source tree:

- compiled executables and object, module, and archive files;
- tracked `build/` and `bin/` output;
- generated Sphinx documentation under `docs/_build/`;
- notebook checkpoints and Python bytecode;
- editor swap files and `.DS_Store`;
- backup and temporary Makefiles;
- personal `_prathap` source snapshots;
- the obsolete `diag.F90.org.reduck96` snapshot; and
- other generated sparse-matrix binaries and module files.

The three solver directories previously occupied approximately 81 MB in the
working tree. The consolidated canonical directory is approximately 6.6 MB.
This does not shrink existing Git history; it only cleans the checked-out tree
and future commits.

Expanded `lanczosKuboCode/.gitignore` so that local build products,
documentation output, notebook state, editor files, and backup files are not
accidentally committed again.

## 4. Build and documentation improvements

Replaced the placeholder solver README with current information covering:

- GRABNES's purpose;
- required compilers and numerical libraries;
- configuration and compilation;
- executable location and invocation;
- cleanup; and
- the link to the first public example.

Added `lanczosKuboCode/make.sys.example`, a portable GNU-oriented starting
configuration using `mpif90`, `gfortran`, OpenMP, ARPACK, LAPACK, and BLAS.
Site-specific configurations remain available under `lanczosKuboCode/Sys/`.

Made the nested clean rules tolerate partially created `MIO/` and `math/`
build directories that do not yet contain Makefiles.

Moved the `action` argument declaration in `Src/MIO/MPI/mpitime.F90` outside
the `TIMER` preprocessor guard. Without this change, GNU Fortran rejected the
subroutine when `TIMER` was not defined because `implicit none` left `action`
undeclared.

Trailing whitespace was removed mechanically from the promoted Fortran and
Makefile sources. No intended program logic was changed by that formatting
cleanup.

## 5. Verification performed

The following checks completed successfully:

- exactly one `lanczosKuboCode*` directory remains;
- no repository text references the removed `_jiaqi` or `_jinwoo` paths;
- no compiled objects, modules, archives, binaries, caches, swap files, or
  listed backup patterns remain in the canonical tree;
- `git diff --cached --check` reports no whitespace errors;
- the example shell and Python files pass syntax checks; and
- the graphene reference data match the Twistronics 2023 result exactly.

The consolidated source build was attempted with the local Homebrew GNU
Fortran and Open MPI installation. Compilation reached the generated MPI helper
program after the `mpitime.F90` fix, but the host linker failed with:

```text
ld: library 'crt1.o' not found
```

This was the result of the original build attempt with an obsolete Intel
Homebrew compiler. A later isolated native Apple Silicon build is documented
below.

## 6. Isolated native build and smoke-test harness (superseded on 2026-10-08)

Added `grabnes_testrun/` temporarily inside the public checkout so
collaborators can reproduce the current Apple Silicon build while keeping all
compiler products under one ignored test subtree. Its local `.gitignore`
excludes `build/`, `bin/`, and generated example results from Git history.

The harness contains a source snapshot with the portability fixes discovered
during testing, an Apple Silicon `make.sys`, the graphene example, a detailed
`BUILD_REPORT.md`, and a top-level `smoke_test.sh`. Run it with:

```sh
cd grabnes_testrun
./smoke_test.sh
```

The test uses native Homebrew GCC 16.2, Open MPI 5.0.11, ARPACK 3.9.1_1, and
OpenBLAS 0.3.34. A clean checked build completed on this laptop, the program
reported `0 errors, 0 warnings`, and the generated graphene bands matched the
reference file byte-for-byte. Compiler warnings remain in the generated
legacy MPI wrappers and are documented in the build report.

The test-copy compatibility changes include modern logical operators, safe
handling of an unassociated single-process MPI list, corrected memory
accounting during deallocation, support for deeply nested source paths, and
selection of the neighbor routine intended for very small cells. These have
not yet been promoted to the canonical solver and should be reviewed first.

The harness is checkout-location independent: its scripts resolve their own
directories, the MIO generators accept long nested paths, and `make.sys`
discovers Homebrew library prefixes instead of hard-coding one installation
path. The required test `make.sys` is explicitly included, while compiler
products and generated results remain ignored.

## 7. Recommended collaborator review before announcement

1. Review the staged consolidation diff, especially the promoted solver source.
2. Build on a clean Linux environment with MPI, ARPACK, LAPACK, and BLAS.
3. Run `examples/01_graphene_bands/run.sh` and compare the resulting bands with
   `reference/bands.dat`.
4. Decide whether the two sparse-diagonalization test programs should be
   integrated into a maintained test suite or removed.
5. Add the planned graphene DOS and twisted-bilayer examples.
6. Add a repository license, `CITATION.cff`, software DOI, and manuscript
   references before the broader announcement.
7. Replace the current internal-group top-level README with a public-facing
   project landing page.

## Recovery

To inspect the repository before consolidation without changing the current
working tree:

```sh
git show pre-code-consolidation
git ls-tree -r --name-only pre-code-consolidation
```

Collaborators should avoid resetting the current work blindly. Review the
staged changes first, then commit them as one consolidation commit or split
them into examples, cleanup, and solver-promotion commits as appropriate.
