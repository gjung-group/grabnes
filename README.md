# GRABNES

GRABNES (GRAphene and Boron Nitride Electronic Structure) is a Fortran code
for electronic-structure and quantum-transport calculations in layered
two-dimensional materials.

The public repository is being prepared for its first announced release. The
four examples below are validated laptop-sized calculations; advanced TAPW,
spin-orbit-coupling, Berry-curvature, and transport workflows should currently
be treated as research functionality that is not covered by any regression
test.

## Quick start

Requirements are a Fortran compiler (tested: GNU Fortran 12.2 and Intel
`ifort` 2021.6), an MPI compiler wrapper, BLAS/LAPACK, ARPACK, and GNU Make. Build the
solver with:

```sh
cd lanczosKuboCode
cp make.sys.example make.sys
# Adjust compiler and library settings in make.sys when necessary.
make
```

The executable is written to `lanczosKuboCode/bin/grabnes`. Run the smallest
validated example from the repository root with:

```sh
cd examples/01_graphene_bands
./run.sh
python3 plot.py
```

The test harness in [`grabnes_testrun`](grabnes_testrun/) builds the same
sources in a temporary directory and compares example results with their
reference data:

```sh
./grabnes_testrun/smoke_test.sh      # build + graphene bands
./grabnes_testrun/run_examples.sh    # build + all four examples + solver checks
```

Tested compilers, measured results, and known limitations are recorded in
[`docs/development/cluster-build-and-validation.md`](docs/development/cluster-build-and-validation.md).

## Validated examples

| Example | Calculation |
| --- | --- |
| [`01_graphene_bands`](examples/01_graphene_bands/) | Pristine graphene bands |
| [`02_graphene_dos`](examples/02_graphene_dos/) | Pristine graphene density of states |
| [`03_twisted_bilayer_bands`](examples/03_twisted_bilayer_bands/) | Small commensurate twisted-bilayer bands |
| [`04_twisted_bilayer_dos`](examples/04_twisted_bilayer_dos/) | Small commensurate twisted-bilayer density of states |

Each directory contains an input file, launcher, plotting script, explanation,
and reference data. GRABNES currently runs on a single MPI process; use OpenMP
threads (Intel build) for parallel execution. See the [examples guide](examples/README.md) for
runtime and executable-selection details.

## Repository layout

- `lanczosKuboCode/`: maintained Fortran solver and build system.
- `examples/`: small validated calculations intended for new users.
- `docs/`: user, theory, API, and development documentation.
- `tests/`: Python mathematical checks and validation plans.
- `tools/`: standalone analysis and debugging utilities.
- `grabnes_testrun/`: regression harness that builds `lanczosKuboCode/` out of
  tree and checks the examples against their reference data.
- `PyBinding/`, `LAMMPSNotebooks/`, and `usefulGeneralScripts/`: legacy and
  research workflows that are not yet part of the minimal supported interface.

## Documentation and tests

Developer documentation can be built as described in
[`docs/development/building-documentation.md`](docs/development/building-documentation.md).
The lightweight TAPW mathematical tests are run with:

```sh
python3 -m pip install -r tests/requirements.txt
python3 -m pytest tests/tapw
```

## Citation

Manuscript citation information will be added before the announced release.
Until then, please cite the repository URL and the exact Git commit used, and
contact the maintainers before relying on unpublished research functionality.
When papers associated with GRABNES are released, their citations and archival
software identifiers will be listed here.

## Development status

Known limitations and implementation notes live under
[`docs/development`](docs/development/). The detailed collaborator log is kept
in [`pre-announcement-changelog.md`](docs/development/pre-announcement-changelog.md)
until it is distilled into a release changelog.

Please report reproducible problems through the repository's GitHub issue
tracker and include the commit hash, compiler version, input file, and log.

## License

No software license has yet been added to this repository. A license must be
selected before the announced release so that reuse terms are explicit.
