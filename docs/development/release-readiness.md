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
| Other compilers and platforms (`ifx`, macOS, Open MPI, MPICH) | Optional | untested; not claimed |
| OpenMP with GNU Fortran | Optional | six `REDUCTION` clauses on pointer arrays |

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
| Further examples (graphene/hBN, semiclassics, transport) | Maintainer decision | deliberately postponed |

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
| Workflow file `.github/workflows/regression.yml` | Complete | syntax checked; uses only repository files and Ubuntu packages |
| First execution on GitHub | Technical work | never run; it could not be reproduced locally (no usable container runtime). Expect to adjust package names or the Open MPI launcher on the first run |

## Licensing, copyright, third-party code

| Item | Status | Note |
| --- | --- | --- |
| Audit of notices and third-party code | Complete | `licensing-audit.md` |
| Choice of license; `LICENSE` file; consistent `pyproject.toml` | Maintainer decision | **release blocker**; needs agreement with the copyright holder of the GPL-marked libraries |
| Copyright holder and headers of the solver sources | Maintainer decision | |
| Commented-out Numerical Recipes routines in `kubo.F90` | Maintainer decision | removal recommended; touches a frozen source file |

## Citation

| Item | Status | Note |
| --- | --- | --- |
| Checklist and template for `CITATION.cff` | Complete | `citation-checklist.md` |
| Author list, version, release date, DOI | Maintainer decision | **release blocker** for a citable release |

## Repository cleanliness

| Item | Status | Note |
| --- | --- | --- |
| Tracked executables removed | Complete | 45 `a.out` files |
| Duplicated notebook checkpoints and `.DS_Store` removed; `.gitignore` extended | Complete | |
| `LAMMPSNotebooks/` (about 380 MB) and `PyBinding/` (about 120 MB) in the software repository | Maintainer decision | research material, not used by the solver or the tests; a separate repository would make the code repository about 20 times smaller. Not moved here |
| Two notebook checkpoints without a live notebook | Maintainer decision | `PyBinding/twistedBilayerGraphene/.ipynb_checkpoints/` |
| Absolute paths in three tutorial notebooks | Optional | outputs of earlier runs |
| Size of the Git history | Optional | unchanged; rewriting history was not attempted |

## Known unsupported functionality (acceptable for the release)

- More than one MPI process.
- OpenMP with GNU Fortran.
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
  55 x 55 cell). The solver reports the asymmetry; `Hopping.Symmetrize` is an
  optional diagnostic remedy and is off by default. A correction exists
  outside this repository and will be transferred by the maintainers.
- `GBNtwoLayers` was neither modified nor validated in this repository.
- Kubo transport, TAPW, spin-orbit terms, Berry curvature, and the
  semiclassical module are retained as research functionality without
  examples.

## Minimum actions before the announcement

1. Choose and add the license (with the copyright holder of the GPL-marked
   code), and make `pyproject.toml` consistent.
2. Supply the author list and create `CITATION.cff`.
3. Push the branch, let the CI workflow run once, and fix what it reports.
4. Decide whether `LAMMPSNotebooks/` and `PyBinding/` ship with the code.
5. Tag the release.
