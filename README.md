# GRABNES

GRABNES (GRAphene and Boron Nitride Electronic Structure) is a Fortran code
for real-space tight-binding calculations of the electronic structure and
quantum transport of graphene, hexagonal boron nitride, and the layered and
moire systems built from them.

It constructs atomistic Hamiltonians for single layers, commensurate twisted
bilayers, and multilayers, and solves them either by exact diagonalization
(band structures, densities of states) or by Lanczos recursion on very large
cells (densities of states of systems with hundreds of thousands to millions
of atoms).

> **Status.** GRABNES is under active development. There is no numbered
> release yet, so please cite the commit you used. The examples in
> [`examples/`](examples/) are tested on every change and reproduce their
> reference data; other functionality is research code of varying maturity
> (see [What is supported](#what-is-supported)).

## Capabilities

Demonstrated by the supported examples and covered by regression tests:

- Tight-binding Hamiltonians of pristine graphene, of commensurate
  twisted-bilayer graphene (down to the magic angle) and of alternating-twist
  trilayer graphene read from a structure file, with nearest-neighbor or
  longer-range intralayer hopping and a distance-dependent two-center
  interlayer hopping.
- The effective moire-potential model of graphene aligned or twisted on hBN.
- Band structures along a k-path and densities of states by exact
  diagonalization, and bands of one valley in the plane-wave reduction (TAPW).
- Density of states of large cells, up to millions of atoms, by Lanczos
  recursion with a random-phase state (the "Kubo" solver), with and without a
  perpendicular magnetic field; reproducible with a fixed seed.

Implemented in the source, but not part of the supported set (use only after
your own validation):

- Kubo time evolution, diffusion, and conductivity.
- Atomistic hexagonal boron nitride and graphene/hBN systems, and multilayer,
  encapsulated, and bulk geometries.
- TAPW beyond the inputs of the examples (both valleys, relaxed cells),
  spin-orbit terms, Berry curvature and Chern numbers, spectral functions,
  sparse diagonalization.

[`docs/development/functionality-status.md`](docs/development/functionality-status.md)
states for each capability what evidence exists.

## Requirements

- A Fortran compiler. Tested: GNU Fortran 11.4 and 12.2, and Intel `ifort`
  2021.6, on Linux x86_64.
- An MPI library with a Fortran compiler wrapper (`mpif90`, `mpiifort`, ...).
  The code is built with MPI but currently runs on **one MPI process**.
- BLAS and LAPACK (reference implementation or MKL).
- ARPACK.
- GNU Make.
- Python 3 for the tests; NumPy and Matplotlib for some checks and for the
  plotting scripts.

Other platforms and compilers (macOS, `ifx`, other MPI libraries) have not
been tested with the current sources.

## Build

```sh
cd lanczosKuboCode
cp make.sys.example make.sys     # adjust compiler and library settings if needed
make
```

The executable is `lanczosKuboCode/bin/grabnes`. `make.sys` is your local,
untracked configuration; `make.sys.example` is a working GNU Fortran setup and
`tests/regression/config/intel.make.sys` an Intel one. Details, including
out-of-tree builds, are in [`lanczosKuboCode/README.md`](lanczosKuboCode/README.md).

Both the GNU and the Intel builds are threaded with OpenMP. A build without
Intel MKL is compiled with `-DNO_MKL` (as in `make.sys.example`); the
shift-invert sparse solver, which needs MKL PARDISO, is then unavailable.

## Quick start

```sh
cd examples/01_graphene_bands
./run.sh              # runs grabnes on Gendata.in, writes generate.bands
python3 plot.py       # writes graphene_bands.png
```

A calculation is controlled by one input file of `keyword value` lines
(`Gendata.in` in the examples) and is started as `grabnes Gendata.in`.

## Supported examples

| Example | Calculation | Size |
| --- | --- | --- |
| [`01_graphene_bands`](examples/01_graphene_bands/) | Bands of pristine graphene | 2 atoms |
| [`02_graphene_dos`](examples/02_graphene_dos/) | DOS of pristine graphene, exact diagonalization on a 30 x 30 k-grid | 2 atoms |
| [`03_twisted_bilayer_bands`](examples/03_twisted_bilayer_bands/) | Bands of commensurate twisted-bilayer graphene at 13.17 degrees | 76 atoms |
| [`04_twisted_bilayer_dos`](examples/04_twisted_bilayer_dos/) | DOS of the same bilayer, 8 x 8 k-grid | 76 atoms |

Each directory holds the input, a launcher, a plotting script, reference data,
and a README that states the Hamiltonian and its parameters. Each runs in a few
seconds on one core. See the [examples guide](examples/README.md).

### Moire examples

Five larger examples on rigid lattices. They take from 20 seconds to 5 minutes on
4 cores and have the same layout as the four above.

| Example | Calculation | Size |
| --- | --- | --- |
| [`05_graphene_hbn_bands`](examples/05_graphene_hbn_bands/) | Bands of graphene aligned with hBN (effective moire model), one valley | 6050 atoms |
| [`06_graphene_hbn_landau_levels`](examples/06_graphene_hbn_landau_levels/) | DOS of the same system at 33 magnetic fields up to 26 T (recursion method) | 387,200 atoms |
| [`07_graphene_hbn_twist_series`](examples/07_graphene_hbn_twist_series/) | DOS of graphene on hBN at twist angles 0, 0.61, 0.92 and 1.87 degrees | about 3.5 million atoms |
| [`08_magic_angle_twisted_bilayer`](examples/08_magic_angle_twisted_bilayer/) | Flat bands of twisted-bilayer graphene at 1.08 degrees | 11,164 atoms |
| [`09_alternating_twist_trilayer`](examples/09_alternating_twist_trilayer/) | Bands of alternating-twist trilayer graphene at 1.54 degrees, structure read from a file | 8322 atoms |

The band examples (05, 08, 09) use the plane-wave reduction (TAPW) for one
valley; for these three inputs it was compared with exact diagonalization in
both directions (every level of either calculation has a partner in the other
within 0.2 meV, in a window of 0.3 eV). Example 05 is part of the standard test
suite; `tests/regression/run_long_examples.sh` reruns 06 to 09 against their
reference data.

## Testing and reproducibility

```sh
./tests/regression/smoke_test.sh      # build + example 01
./tests/regression/run_examples.sh    # build + examples 01 to 05 + further checks
python3 -m pytest tests/tapw         # TAPW linear-algebra conventions (Python only)
```

The same commands run automatically on GitHub for every push to the
development and main branches (`.github/workflows/regression.yml`).

The scripts build the solver in a temporary directory, never in the source
tree, run the example inputs, and compare the results numerically with the
reference data (bands within 2e-6 eV, DOS within 1e-9 relative). They also
compare the assembled Hamiltonian, entry by entry, with one constructed
independently in Python. [`tests/regression/README.md`](tests/regression/README.md)
describes every check; measured results for both compilers are in
[`docs/development/cluster-build-and-validation.md`](docs/development/cluster-build-and-validation.md).

## What is supported

| Level | Meaning | Applies to |
| --- | --- | --- |
| Supported | Runs, is documented, and is reproduced by the test suite | The nine examples above (06 to 09 by `run_long_examples.sh`) |
| Independently validated | Additionally checked against an analytic result or a second implementation | Hamiltonian assembly and bands of graphene, hBN monolayer, and twisted bilayers; Kubo DOS of graphene; TAPW bands of examples 05, 08 and 09 against exact diagonalization |
| Research / experimental | Present in the source without sufficient validation for general use | Everything else, including transport, spin-orbit, Berry curvature, atomistic graphene/hBN, and the models of the examples used with other inputs |

## Known limitations

- One MPI process only; a run with more is refused. Parallelism is by OpenMP.
- The Kubo solver is stochastic: results vary from run to run unless
  `Kubo.SetSeed .true.` and `Kubo.SeedValue` are given.
- The bond term of the effective moire-potential model of graphene on hBN
  (`Moire.OffDiag`) gives a hopping table that is not exactly Hermitian; the
  solver reports the size of the asymmetry. `Moire.OffDiagMidpoint .true.`
  evaluates the term at the bond midpoint and makes it Hermitian. See the
  validation document.
- Many input keys exist only for research workflows. Nearly all have a one-line
  description in [`docs/user-guide/input-keys.md`](docs/user-guide/input-keys.md),
  but only the keys of the examples are tested; the example inputs and READMEs
  are the reference for the supported ones.

## Documentation

- [`examples/README.md`](examples/README.md): the supported calculations.
- [`lanczosKuboCode/README.md`](lanczosKuboCode/README.md): building and running.
- [`tests/regression/README.md`](tests/regression/README.md): the regression suite.
- [`docs/development/`](docs/development/): validation report, functionality
  status, release-readiness checklist, licensing audit, and development log;
  [`docs/development/README.md`](docs/development/README.md) is the index.
- [`docs/user-guide/input-file.md`](docs/user-guide/input-file.md): format of the input file, names of the
  keys (the names used before October 2026 keep working), and the messages the solver prints about it.
- [`docs/user-guide/input-keys.md`](docs/user-guide/input-keys.md): every input key with its type,
  default, and the routine that reads it (generated from the source).
- [`docs/user-guide/stacking-displacements.md`](docs/user-guide/stacking-displacements.md): the stacking
  displacements of the moiré terms (convention, and computing them from the rigid and relaxed structures).
- [`docs/user-guide/`](docs/user-guide/) and [`docs/theory/`](docs/theory/):
  notes on TAPW and spin-orbit input (research functionality).
- [`tools/`](tools/): analysis scripts, including the independent Hamiltonian
  check.
- `usefulGeneralScripts/`: research scripts of the group that are not part of
  the supported interface.

## Authors and maintenance

The authors of GRABNES are Jeil Jung, Rafael Martinez-Gordillo, and Nicolas
Leconte. The code grew out of a tight-binding and Kubo transport program
written by Rafael Martinez-Gordillo, whose input/output and math libraries it
still uses, and was developed further in the group of Jeil Jung.

Nicolas Leconte is the current developer and the maintainer of this
repository; questions and reports should be addressed to him through the
issue tracker.

## Citation

Please cite GRABNES by its authors, the repository, and the exact version or
Git commit you used. The metadata are in [`CITATION.cff`](CITATION.cff), which
GitHub displays under "Cite this repository":

> J. Jung, R. Martinez-Gordillo, and N. Leconte, *GRABNES: Graphene and Boron
> Nitride Electronic Structure*, https://github.com/gjung-group/grabnes,
> commit `<hash>`.

A software paper and an archived release with a DOI do not exist yet; this
section will name them when they do.

### Publications to cite with the code

Please also cite the publications behind the methods and models you use.

Methods of the code (real-space lattice models, Lanczos recursion):

- N. Leconte, A. Ferreira, and J. Jung, *Efficient multiscale lattice
  simulations of strained and disordered graphene*, Semicond. Semimet. **95**,
  35-99 (2016), [doi:10.1016/bs.semsem.2016.04.002](https://doi.org/10.1016/bs.semsem.2016.04.002).

Twisted and multilayer graphene:

- N. Leconte, S. Javvaji, J. An, A. Samudrala, and J. Jung, *Relaxation effects
  in twisted bilayer graphene: A multiscale approach*, Phys. Rev. B **106**,
  115410 (2022), [doi:10.1103/PhysRevB.106.115410](https://doi.org/10.1103/PhysRevB.106.115410)
  (twisted-bilayer model; example 08).
- N. Leconte, Y. Park, J. An, A. Samudrala, and J. Jung, *Electronic structure
  of lattice relaxed alternating twist tNG-multilayer graphene: from few layers
  to bulk AT-graphite*, 2D Mater. **9**, 044002 (2022),
  [doi:10.1088/2053-1583/ac8a00](https://doi.org/10.1088/2053-1583/ac8a00)
  (alternating-twist multilayers; example 09).
- K. Shin, N. Leconte, J. Jung, and H. Min, *Supermoire reconstruction and
  domain-resolved valley topology in helical multilayer graphene*, 2D Mater.
  **13**, 045021 (2026), [doi:10.1088/2053-1583/aeb140](https://doi.org/10.1088/2053-1583/aeb140)
  (helical multilayers).
- P. Moon and M. Koshino, *Energy spectrum and quantum Hall effect in twisted
  bilayer graphene*, Phys. Rev. B **85**, 195458 (2012),
  [doi:10.1103/PhysRevB.85.195458](https://doi.org/10.1103/PhysRevB.85.195458)
  (two-center interlayer hopping; examples 03 and 04).

Graphene on hexagonal boron nitride:

- J. Jung, A. Raoux, Z. Qiao, and A. H. MacDonald, *Ab initio theory of moire
  superlattice bands in layered two-dimensional materials*, Phys. Rev. B **89**,
  205414 (2014), [doi:10.1103/PhysRevB.89.205414](https://doi.org/10.1103/PhysRevB.89.205414)
  (effective moire-potential model; examples 05 to 07).
- N. Leconte and J. Jung, *Commensurate and incommensurate double moire
  interference in graphene encapsulated by hexagonal boron nitride*, 2D Mater.
  **7**, 031005 (2020), [doi:10.1088/2053-1583/ab891a](https://doi.org/10.1088/2053-1583/ab891a)
  (graphene encapsulated by hBN, double moire).
- J. An, N. Leconte, S. Javvaji, Y. Park, and J. Jung, *Atomistic substrate
  relaxation effects in the band gaps of graphene on hexagonal boron nitride*,
  [arXiv:2602.18893](https://arxiv.org/abs/2602.18893) (2026) (atomistic
  graphene/hBN with a relaxed substrate).
- W. Shi, S. Kahn, N. Leconte, T. Taniguchi, K. Watanabe, M. Crommie, J. Jung,
  and A. Zettl, *High-order fractal quantum oscillations in graphene/BN
  superlattices in the extreme doping limit*, Phys. Rev. Lett. **130**, 186204
  (2023), [doi:10.1103/PhysRevLett.130.186204](https://doi.org/10.1103/PhysRevLett.130.186204)
  (graphene/hBN in a magnetic field).
- Z. Wang, P. M. Perez-Piskunow, C. P. Y. Wong, M. Holwill, H. Sun, W. Fu,
  J. Hu, Y. Zheng, X. Zhang, T. Taniguchi, K. Watanabe, A. Ariando, G. Yu,
  L. Lin, K. E. J. Goh, S. Roche, J. Jung, K. S. Novoselov, and N. Leconte,
  *Beyond the flat-band paradigm: open orbits as a route to correlated states
  in dispersive moire minibands*, Proc. Natl. Acad. Sci. USA, accepted (2026),
  [arXiv:2610.06110](https://arxiv.org/abs/2610.06110) (graphene/hBN
  minibands and their orbits in a magnetic field).

## License

GRABNES is free software, distributed under the terms of the GNU General
Public License, version 3 or (at your option) any later version. The full
text is in [`LICENSE`](LICENSE). The program comes without any warranty.

Some files were written by other authors and keep their own notices (GNU GPL,
GNU LGPL, two-clause BSD); they are listed in
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md). The reasoning behind
the choice of license is recorded in
[`docs/development/licensing-audit.md`](docs/development/licensing-audit.md).

## Contact and contributions

Please report reproducible problems through the GitHub issue tracker of the
repository, with the commit hash, compiler version, input file, and log.

Before changing the code, read [`AGENTS.md`](AGENTS.md): the rules a change
must follow and the validation a pull request must show. It is written for
coding agents and human contributors alike; the full conventions are in
[`docs/development/coding-conventions.md`](docs/development/coding-conventions.md)
and [`docs/development/output-style.md`](docs/development/output-style.md).
