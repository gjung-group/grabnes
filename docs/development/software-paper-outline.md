# Outline for a software paper on GRABNES

Working outline for a possible article in *Computer Physics Communications*
(Computer Programs in Physics). It is a plan, not a manuscript: it lists what
each section would contain and which material already exists. No such paper
has been written or submitted; nothing here is a publication claim.

Proposed software authors, as in `CITATION.cff`: Jeil Jung, Rafael
Martinez-Gordillo, Nicolas Leconte. The author list of the article, the
author roles, and any contribution statement are for the authors to decide
and are not anticipated here. Citations and results are left open.

Legend for the state of the supporting material:

- **[ready]** validated and documented in this repository;
- **[implemented]** present in the source, can be described as functionality
  but should not be presented with results until validated;
- **[future]** better reserved for later work.

## Program summary (CPC front matter)

To be filled in: program title, repository and archive link, licensing
provisions (GPL-3.0-or-later, see `licensing-audit.md`), programming language (Fortran 90
and later, Python for tests and tools), external libraries (MPI, BLAS/LAPACK,
ARPACK), nature of the problem, solution method, restrictions (one MPI
process; OpenMP with one compiler family).

## 1. Scientific motivation

- Moire and layered two-dimensional materials: why atomistic real-space
  tight-binding is needed next to continuum models (large cells, lattice
  relaxation, lack of translational symmetry, disorder, magnetic field).
- Scale of the problem: thousands to millions of orbitals; where exact
  diagonalization ends and recursion methods begin.
- Position relative to existing codes. *To be written by the authors; no
  comparison with other packages has been made in this repository.*

## 2. Scope and capabilities

- One paragraph per capability, each tagged with its state, following
  `functionality-status.md`.
- **[ready]** graphene and twisted-bilayer graphene bands and DOS; Kubo
  recursion DOS.
- **[implemented]** hBN and graphene/hBN systems, multilayers, bulk, TAPW,
  spin-orbit terms, Berry curvature and Chern numbers, spectral functions,
  sparse diagonalization, magnetic field, semiclassical analysis.
- A table that separates "demonstrated in this paper" from "available in the
  code".

## 3. Hamiltonian construction and parameterizations

- Single-orbital pz tight-binding Hamiltonian; internal energy unit `g0`;
  sign convention of the stored hoppings. **[ready]**
- Intralayer models: nearest neighbor; F2G2-type model on five neighbor
  shells; eight-shell model; two-center (Koshino) form. Table of the values
  used and of the shells they act on. **[ready as implementation; the
  literature provenance of the F2G2 values and of the rule for `g0` is still
  to be confirmed by the authors]**
- Interlayer model: distance-dependent two-center form; role of `vpppi0` and
  its model-dependent default (2.7 eV with the Koshino intralayer model,
  3.5 eV otherwise); in the F2G2-type models the Dirac velocity is set by the
  intralayer parameters. **[ready]**
- Effective moire-potential model of graphene on hBN (on-site and bond
  terms). **[implemented; on-site and bond terms agree with the published
  expressions, but the assembled bond terms are not Hermitian in this version
  and the model is not part of the supported set]**
- Explicit graphene/hBN heterostructures and multilayers. **[implemented]**
- Bloch convention: lattice-vector phases; consequences for periodicity in
  reciprocal space and for the choice of unit-cell images. **[ready]**

## 4. Atomic geometries and moire systems

- Commensurate twisted cells from integer indices `(m,n)`; number of atoms
  and twist angle; supercells; layers and stacking. **[ready]**
- Structures read from files (relaxed geometries). **[implemented]**
- Neighbor search: exhaustive enumeration over periodic images with exact
  lattice translations, linked-cell binning, shell table and
  `TB.NeighLevels`. **[ready]**

## 5. Numerical methods

- Exact diagonalization of the Bloch Hamiltonian (LAPACK). **[ready]**
- Lanczos recursion and continued-fraction DOS with random-phase states;
  statistical error and its scaling with system size; seeding. **[ready for
  the DOS]**
- Time evolution and Kubo conductivity. **[implemented; future validation]**
- Sparse diagonalization with ARPACK; plane-wave reduction (TAPW).
  **[implemented]**

## 6. Electronic-structure calculations

- Bands of graphene and of twisted-bilayer graphene; Dirac velocities for
  the intralayer models and their renormalization with twist angle (measured
  at 13.17, 6.01, and 3.89 degrees). **[ready]**
- Band-structure-derived quantities: Berry curvature, Chern numbers,
  semiclassical quantization. **[implemented; future]**

## 7. Density of states and transport

- DOS by exact diagonalization; DOS by recursion compared with the exact
  broadened DOS of graphene. **[ready]**
- Transport. **[future]**: needs norm-conservation and ballistic-spreading
  checks before any result is shown.

## 8. Validation and benchmark examples

- The four supported examples and their reference data. **[ready]**
- Independent validation: entry-by-entry comparison of the Hamiltonian with a
  second implementation; Hermiticity; invariances; analytic limits (decoupled
  layers, hBN two-band formula, Dirac-point degeneracy). **[ready]**
- What the validation does not establish (physical adequacy of parameters;
  models outside the tested set). To be stated explicitly.

## 9. Computational performance

- Setup cost (neighbor search and Hamiltonian assembly): linear scaling
  measured from 3e3 to 1.2e7 atoms; memory per stored hopping. **[ready]**
- Cost of the solvers (diagonalization, recursion steps), thread scaling with
  OpenMP. **[not measured yet; required for this section]**

## 10. Software architecture and reproducibility

- Source layout: MIO library (input, memory, timers), math library, solver
  modules. Input format. Build system and out-of-tree builds. **[ready]**
- Regression suite, tolerances and their justification, compilers tested,
  continuous integration. **[ready; the CI workflow has not yet run]**
- Recording of commit, compiler, flags, and executable checksum in every test
  run. **[ready]**

## 11. Applications and limitations

- Representative applications. *To be selected by the authors from published
  work; none is asserted here.*
- Limitations: single MPI process; OpenMP with one compiler family;
  stochastic recursion; validated model set; undocumented research keywords.

## Material still needed before writing

1. Affiliations and identifiers of the authors, and an archived release with
   a DOI (see `citation-checklist.md`).
2. Confirmed literature sources for every parameter set quoted.
3. Solver timings and thread scaling.
4. A decision on which implemented capabilities are described, and for each
   one the validation that would allow results to be shown.
