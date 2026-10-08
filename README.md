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

> **Status.** This repository is a release candidate that has not yet been
> announced. Four example calculations are supported and reproducible; the
> remaining functionality is research code of varying maturity (see
> [What is supported](#what-is-supported)).

## Capabilities

Demonstrated by the supported examples and covered by regression tests:

- Tight-binding Hamiltonians of pristine graphene and of commensurate
  twisted-bilayer graphene, with nearest-neighbor or longer-range intralayer
  hopping and a distance-dependent two-center interlayer hopping.
- Band structures along a k-path and densities of states by exact
  diagonalization.
- Density of states of large cells by Lanczos recursion with a random-phase
  state (the "Kubo" solver), reproducible with a fixed seed.

Implemented in the source, but not part of the supported set (use only after
your own validation):

- Kubo time evolution, diffusion, and conductivity.
- Hexagonal boron nitride and graphene/hBN systems, including an effective
  moire-potential model of graphene on hBN, and multilayer, encapsulated, and
  bulk geometries.
- Plane-wave-reduced Hamiltonians (TAPW), spin-orbit terms, Berry curvature
  and Chern numbers, spectral functions, sparse diagonalization, magnetic
  fields.
- Semiclassical analysis of band structures (Fermi contours, open orbits,
  Onsager quantization, cyclotron masses, Berry-phase and orbital-moment
  corrections, magnetic breakdown).

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

With GNU Fortran the OpenMP directives are not compiled (the compiler rejects
some of them), so the GNU build is serial; the Intel build is threaded.

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

## Testing and reproducibility

```sh
./tests/regression/smoke_test.sh      # build + example 01
./tests/regression/run_examples.sh    # build + the four examples + further checks
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
| Supported | Runs, is documented, and is reproduced by the test suite | The four examples above |
| Independently validated | Additionally checked against an analytic result or a second implementation | Hamiltonian assembly and bands of graphene, hBN monolayer, and twisted bilayers; Kubo DOS of graphene |
| Research / experimental | Present in the source without sufficient validation for general use | Everything else, including transport, TAPW, spin-orbit, Berry curvature, semiclassics, graphene/hBN and multilayer models |

## Known limitations

- One MPI process only; a run with more is refused. Parallelism is by OpenMP,
  with the Intel build.
- The Kubo solver is stochastic: results vary from run to run unless
  `setSeed .true.` and `seedValue` are given.
- The effective moire-potential model of graphene on hBN assembles a hopping
  table that is not exactly Hermitian; the solver reports the size of the
  asymmetry. See the validation document.
- Many input keywords exist only for research workflows and are undocumented.
  The example inputs and READMEs are the reference for the supported ones.

## Documentation

- [`examples/README.md`](examples/README.md): the supported calculations.
- [`lanczosKuboCode/README.md`](lanczosKuboCode/README.md): building and running.
- [`tests/regression/README.md`](tests/regression/README.md): the regression suite.
- [`docs/development/`](docs/development/): validation report, functionality
  status, release-readiness checklist, licensing audit, and development log;
  [`docs/development/README.md`](docs/development/README.md) is the index.
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
section will name them when they do. Please also cite the publications behind
the models you use; those of the supported examples are named in their
READMEs.

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
