# Functionality status

This page classifies what the public GRABNES solver does by the evidence that
exists for it in this repository, as of October 2026. "Validated" never means
more than what the listed test establishes. Details and numbers are in
[`cluster-build-and-validation.md`](cluster-build-and-validation.md).

## Validated

Compiles, runs, and has a numerical regression test or an independent check in
`tests/regression/run_examples.sh`.

| Capability | Evidence | Limits |
| --- | --- | --- |
| Tight-binding Hamiltonian of layered carbon systems built by the default neighbor search | Entry-by-entry agreement (4e-16 eV) with an independently enumerated Hamiltonian for pristine graphene, hBN monolayer, and twisted bilayers of 28, 76, and 244 atoms; Hermiticity; reciprocal-lattice periodicity; invariance under unit-cell images | Models covered: nearest neighbor, F2G2 (5 shells), eight-shell, Koshino intralayer; two-center interlayer form |
| Exact-diagonalization bands (`Diag.Calc`, `Calculate.Bands`) | Examples 01 and 03; independent eigenvalues within the six-decimal output (5e-7 eV); Dirac-point degeneracy | Dense `ZHEEV`; single MPI process |
| Exact-diagonalization DOS (`Calculate.DOS`) | Examples 02 and 04 (stored references, agreement 1e-13 between compilers) | Reference-based only; no analytic check yet |
| Kubo recursion DOS (`Kubo.Calc`, `Calculate.OnlyDOS`) | 180000-atom graphene against the exact broadened DOS, within the statistical error | Graphene only; stochastic unless `setSeed .true.` |
| Neighbor-shell control `TB.NeighLevels` 1 to 8 | All eight values against the independent model; out-of-range values refused | Shell radii assume the graphene lattice parameter |
| Model-dependent default of `vpppi0` | 2.7 eV with `KoshinoIntralayer`, 3.5 eV otherwise, each checked against the independent model and in the log | |
| Build with GNU Fortran 12 and Intel `ifort` 2021.6 | All checks pass with optimized, checked, and Intel builds | One MPI process |
| Magnetic field (`MagField.Integer`) | Landau levels of a 1152-atom graphene cell at 137 T: zero mode exact and twofold, E_1 to E_3 within 0.16 %, 0.31 %, 0.47 % of the Dirac formula, independent of k | One field, graphene |
| hBN on-site energies | Band edges of monolayer hBN at K equal the two on-site energies | |
| Sparse diagonalization (`sparseDiagSolver`, ARPACK) | The 20 levels nearest zero energy of a 364-atom twisted bilayer equal the dense result at three k-points | Without shift-invert |
| Every model switch | `tests/regression/model_survey.py`: more than 400 small cases; each must end as recorded and reproduce the fingerprint of its Hamiltonian with the GNU and Intel builds | The survey records the present state of a switch (working, without effect, refused, non-Hermitian), not its physical correctness; 51 cases depend on the compiler and are listed as such |

## Partially validated

Runs and passes consistency checks, but lacks an independent numerical
validation of its physics.

| Capability | What was checked | What is missing |
| --- | --- | --- |
| Systems periodic along z (`Bulk`, `BulkSmall`) | Twisted bulk cell: complete, symmetric neighbor list with exact translations; Hermitian tables that reproduce the solver's bands | Hopping values against a model; a physical reference (graphite) |
| `TrilayerBasedOnMoireCell`, `MoireEncapsulatedBilayerBasedOnMoireCell` | Same structural checks on 76- and 114-atom cells | Inputs were guessed: with them the trilayer came out with zero interlayer coupling. Canonical inputs and expected hoppings are needed |
| Effective moire-potential model of graphene on hBN (`MoirePotential`, `MoireJeil`, `MoireOffDiag`; Jung et al., Phys. Rev. B 89, 205414) | With the inputs of the tutorial notebook, the on-site energies (1e-16 eV) and the bond terms (4e-16 eV) equal the published expressions with the published coefficients, on 11 x 11 and 55 x 55 cells | The bond term is evaluated at the atom a bond starts from, so the stored table is not Hermitian (up to 1.7 meV for 55 x 55). The solver warns about it; `MoireOffDiagMidpoint .true.` evaluates the term at the bond midpoint, which makes the table Hermitian, and is off by default. Not in the regression suite |
| `NeighList` (selected with `Neigh.fastNNnotsquare .false.`) | Bands identical to the default search for graphene supercells, shells 1 to 3 | Refuses twisted bilayers; not examined beyond single-layer cells |
| OpenMP (GNU and Intel builds) | Same results with 1 and 4 threads for every check | Speed-up not measured |
| Large-system setup | Neighbor search and Hamiltonian assembly timed up to 11.5 million atoms; linear scaling | Assembly only; no large calculation was validated |

## Experimental or unvalidated

Present in the source; not established by any test here. Do not rely on these
without your own validation.

| Capability | Status |
| --- | --- |
| Kubo time evolution, diffusion, conductivity | Never tested |
| TAPW (reduced plane-wave basis) | Only the linear-algebra conventions are tested in Python (`tests/tapw`); the Fortran path is not run |
| Spin-orbit coupling terms | Validation plan only (`tests/soc/README.md`) |
| Berry curvature, Chern numbers | Known unresolved issue (`known-issues/berry-curvature.md`) |
| Sparse (ARPACK) diagonalization | Not run |
| Spectral functions, PDOS, 3D bands | Not run |
| Magnetic field, Haldane term, Hubbard self-consistency | Not run |
| `GBNtwoLayers` (explicit graphene/hBN heterostructures, normally with a structure read through `ReadXYZ`) | Under active development by the maintainers; left unmodified and not validated in this repository |
| `Graphene_Over_BN`, `BilayerGraphene_Over_BN` | **Legacy.** Geometry generators kept for backward compatibility. With them the hBN layer is coupled to graphene through the in-plane C-B and C-N hopping constants, which is not a physical interlayer model. For graphene/hBN use `GBNtwoLayers` or the effective moire-potential model |
| Ribbons, hybrid, and other `ReadXYZ` systems | Not run |
| Strained, relaxed, or lattice-mismatched structures | The intralayer search radius is that of the rigid lattice. The solver warns when pairs lie close to it; `Neigh.IntralayerRadius` (in Angstrom) sets a larger one. `TB.NeighLevels 6` is marginal for a 2 % lattice mismatch |
| `MayouIntralayer` | Runs through the same code as `KoshinoIntralayer`; its intended `vpppi0` was not specified and stays at 3.5 eV |

## Not supported

| Item | Behavior |
| --- | --- |
| More than one MPI process | Refused with an error (domain decomposition is disabled in the source) |
| Shift-invert sparse solver (`Diag.SparseUseShift`, `Diag.SparseCount`) without Intel MKL | Stops with a message in builds compiled with `-DNO_MKL` |
| `Neigh.fastNN`, `Neigh.fastNNnotsquareNotRectangle` | Refused with an error (their binning overruns its arrays) |
| `TB.NeighLevels` outside 1 to 8 | Refused with an error |
