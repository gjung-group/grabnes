# Release-readiness checklist

State of the repository for a first announced release, October 2026.
Each item is marked **Complete**, **Maintainer decision**, **Technical work**,
or **Optional**. The acceptance criterion for this release is that the four
example calculations build, run, and reproduce their reference data;
capabilities outside that set are released as research code.

## Source code

| Item | Status | Note |
| --- | --- | --- |
| One canonical solver tree | Complete | `lanczosKuboCode/` |
| Compiles with GNU Fortran 12.2 and Intel `ifort` 2021.6 | Complete | |
| Undefined behavior found in the tested paths removed | Complete | listed in `cluster-build-and-validation.md` |
| Neighbor search: complete, symmetric, exact lattice translations | Complete | |
| Unsupported run modes refused with a message | Complete | several MPI processes; two legacy neighbor switches; shell counts outside 1 to 8 |
| Physics-code freeze for this release | Complete | no further source changes planned |
| Remaining "may be used uninitialized" compiler warnings (about 120, mostly `ham.F90`) | Optional | not on the paths of the four examples (checked with Valgrind) |
| Dead code (`fastNNnotsquareSmall`, `NeighListOld`, large commented blocks) | Optional | |

## Build and portability

| Item | Status | Note |
| --- | --- | --- |
| Documented build works from a clean clone | Complete | `cp make.sys.example make.sys; make` |
| Build independent of the working directory; out-of-tree builds | Complete | |
| No machine-specific path in the build files that are used | Complete | |
| Historical site configurations under `lanczosKuboCode/Sys/` contain absolute paths | Optional | kept as starting points; could be removed |
| Other compilers and platforms (`ifx`, macOS, MPICH) | Optional | untested; not claimed. GNU Fortran 11.4 with Open MPI 4.1.2 on Ubuntu 22.04 is exercised by the CI workflow |
| OpenMP with GNU Fortran | Complete | |

## Numerical validation

| Item | Status | Note |
| --- | --- | --- |
| Four examples reproduce their references (GNU and Intel) | Complete | bands byte-identical; DOS within 2e-13 |
| Hamiltonian checked against an independent implementation | Complete | graphene, hBN monolayer, three twisted cells, four intralayer models |
| Kubo recursion DOS against the exact result | Complete | graphene, 180000 atoms |
| Model-dependent default of `vpppi0` | Complete | |
| Validation of capabilities outside the four examples | Optional | not required for this release; status in `functionality-status.md` |

## Examples

| Item | Status | Note |
| --- | --- | --- |
| Four examples run with their launchers | Complete | |
| Plotting scripts run | Complete | Matplotlib |
| READMEs state the Hamiltonian and parameters actually used | Complete | |
| No private or cluster-specific path | Complete | |
| Further examples (graphene/hBN, transport) | Maintainer decision | deliberately postponed |

## Documentation

| Item | Status | Note |
| --- | --- | --- |
| Root README: purpose, capabilities, build, quick start, tests, limits | Complete | |
| Validation report and functionality status | Complete | |
| Input keyword reference | Technical work | only the keywords of the examples are documented; the solver reads several hundred |
| Sphinx/Doxygen build of the documentation | Optional | configuration present; not rebuilt in this pass |
| Consolidating the development log into a release changelog | Optional | |

## Continuous integration

| Item | Status | Note |
| --- | --- | --- |
| Workflow `.github/workflows/regression.yml` | Complete | builds with GNU Fortran and Open MPI on Ubuntu 22.04, runs the four examples and the regression checks against the reference data, then the TAPW tests; read-only permissions |
| Result of the workflow on GitHub | Complete | passes on `nicolas/development` (first green run: commit `445656e`, GNU Fortran 11.4, Open MPI 4.1.2, Ubuntu 22.04): build, four examples, regression checks, plotting scripts, TAPW tests |

## Licensing, copyright, third-party code

| Item | Status | Note |
| --- | --- | --- |
| License chosen and added | Complete | `GPL-3.0-or-later`; `LICENSE`, `COPYING.LESSER`, `THIRD_PARTY_LICENSES.md`; `pyproject.toml` consistent |
| Third-party notices preserved | Complete | GPL (MIO, math), LGPL (four routines), BSD (`fracToCart.py`) |
| Commented-out Numerical Recipes routines | Complete | removed from `kubo.F90` (comments only) |
| Written agreement of the co-authors with the license for the sources without notice | Maintainer decision | see "Open points" in `licensing-audit.md` |
| Copyright headers in the solver sources | Optional | |

## Citation

| Item | Status | Note |
| --- | --- | --- |
| `CITATION.cff` | Complete | three authors, maintainer as contact, schema-valid |
| Version, release date, DOI, ORCID identifiers, affiliations | Maintainer decision | to be added with the first release; not invented |

## Repository cleanliness

| Item | Status | Note |
| --- | --- | --- |
| Tracked executables removed | Complete | 45 `a.out` files |
| Duplicated notebook checkpoints and `.DS_Store` removed; `.gitignore` extended | Complete | |
| `LAMMPSNotebooks/` and `PyBinding/` | Complete | removed from the repository (about 500 MB of research notebooks and data that the solver and the tests did not use); they remain in the Git history |
| Absolute paths in three tutorial notebooks | Optional | outputs of earlier runs |
| Size of the Git history | Optional | unchanged; rewriting history was not attempted |

## Known unsupported functionality (acceptable for the release)

- More than one MPI process.
- `Neigh.fastNN`, `Neigh.fastNNnotsquareNotRectangle`.
- `Graphene_Over_BN`: legacy geometry generator; with it the hBN layer is
  coupled through in-plane hopping constants. Use `GBNtwoLayers` with a
  structure file, or the effective moire-potential model, for graphene/hBN.

## Scientific limitations (to be stated, not blockers)

- The supported examples use a nearest-neighbor intralayer model; they are
  regression cases, not recommendations of a parameter set.
- Literature provenance of the F2G2 intralayer values and of the rule for
  `g0` is unconfirmed.
- Effective graphene/hBN model: the bond terms are evaluated at one end of
  each bond, so the stored hopping table is not Hermitian (up to 1.7 meV in a
  55 x 55 cell). The solver warns about the asymmetry;
  `MoireOffDiagMidpoint .true.` evaluates the term at the bond midpoint and
  makes the table Hermitian. It is off by default and covers the plain
  effective model only; the other position-dependent bond terms
  (`tBGOffDiag`, `GBNOffDiag`) have no such option yet.
- `GBNtwoLayers` is not validated in this repository.
- Strained or strongly corrugated structures can have pairs of the last
  intralayer shell beyond the rigid-lattice search radius. The solver warns;
  `Neigh.IntralayerRadius` sets a larger radius.
- Kubo transport, TAPW, spin-orbit terms, and Berry curvature are retained as
  research functionality without examples.

## Planned developments

Recorded so that they are not lost; none is required for the first release.

- **Spin terms in exact tight-binding calculations.** `ZeemanTerm`,
  `PseudoZeemanTerm`, `IsingSOCterm`, and `RashbaSOCterm` are implemented in the
  TAPW path only. Without `useTAPW` the first three have no effect (the solver
  prints a warning) and Rashba is refused. They are to be added to the exact
  diagonalization as well.
- **Self-consistent Hubbard calculations (`scf.F90`, `SpinPolarized` with
  `EnableSCF`).** The module was never completed or tested; a spin-polarized run
  with the default `EnableSCF .true.` ends with a segmentation fault. Finishing
  and validating it is an important item for future development.
- **Models refused until repaired** (the solver stops with a message):
  `TypeOfBL Jeil`, `RandomStrain`, `printBubble`, `realisticBubbles` without its
  list of centres, and every combination of switches that gives non-finite
  matrix elements (found so far: `TypeOfBL BLKaxiras` on a twisted bilayer,
  `singleLayerXYZ`, `MoireTrilayer` on the effective model).
- **Hermitian evaluation of the remaining position-dependent bond terms**
  (`tBGOffDiag`, `GBNOffDiag`, the non-plain variants of `MoireOffDiag`), in the
  way `MoireOffDiagMidpoint` does it for the plain effective model.
- **Hopping part of the PIA spin-orbit term.** `ApplyPIAHopping` exists but is
  not called and was never tested, so `PIASOCterm` acts through its on-site
  part only. To be wired in and validated.
- **Geometry names that are tested in the source but not accepted:**
  `TwistedBilayerBasedOnMoireCellRectangular`, `BLtoSLYoungju`, `Hybrid`.
- **Keys read with different defaults in different places.** Most are the
  parameter sets of different models read under one name (`SingleLayert2KSL`,
  `BilayertAB1`, `CAA`, `PhiAA`, ...), which is harmless but undocumented. A
  few are inconsistencies whose resolution changes a default and is therefore
  left to a decision:
  `SuperCell` (1 everywhere, but 60 where `ham.F90` computes the moire length
  for `periodicStrain` and the strained-moire terms);
  `CellSize` (50; 55 for the position of `MoireBilayerElectricFieldInvert`;
  1 in two spectral-function routines);
  `InterlayerDistance` (3.22; 3.35 in two branches of `ham.F90`);
  `MagField.Integer` (0; 1 where `ham.F90` reads it for `FrankMagneticField`);
  `Neigh.LayerNeighbors` (0, 1 or 2 depending on the routine);
  `Epsilon` (0.01; 0.001 in `DiagHamChern`);
  `bubbleSigmaR` (1.0 and 1.42); `vpppi0` (2.7 and 3.5);
  `MoirePotCabG` (0.002235 and 0.001987).
- **MPI domain decomposition**, disabled at present.

## Minimum actions before the announcement

1. Obtain the co-authors' agreement with the license (a written record is
   sufficient).
2. Merge the development branch into `main` once the CI run is green.
3. Tag the first release, archive it to obtain a DOI, and add version, date,
   and DOI to `CITATION.cff`.

Maintainer: Nicolas Leconte.
