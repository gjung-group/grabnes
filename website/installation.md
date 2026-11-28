# Installation

The code is maintained on GitHub, where you are free to clone or fork the
stable or development versions of the code (see the GitHub documentation on
how to do this):
[github.com/gjung-group/grabnes](https://github.com/gjung-group/grabnes).

```{note}
The requirements, build and test instructions on this page are quoted from the
[README](https://github.com/gjung-group/grabnes/blob/main/README.md) and
[`lanczosKuboCode/README.md`](https://github.com/gjung-group/grabnes/blob/main/lanczosKuboCode/README.md)
of the code repository. GRABNES has no numbered release yet; see
[](status.md) for what is supported and [](citing.md) for how to cite the
commit you used.
```

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
`tests/regression/config/intel.make.sys` an Intel one.

The files under `Sys/` are historical site configurations with absolute paths
and are kept only as starting points.

Both the GNU and the Intel builds are threaded with OpenMP. A build without
Intel MKL is compiled with `-DNO_MKL` (as in `make.sys.example`); the
shift-invert sparse solver, which needs MKL PARDISO, is then unavailable.

The build does not depend on the current directory and can be placed outside
the source tree:

```sh
make -C /path/to/grabnes/lanczosKuboCode \
     MAKE_SYS=/path/to/make.sys BUILD_DIR=/scratch/build BIN_DIR=/scratch/bin
```

Measured builds on a cluster with both compilers are described in
[](../docs/development/cluster-build-and-validation.md). To remove the build:

```sh
make clean
```

## Run

GRABNES accepts the input filename as its first argument:

```sh
bin/grabnes Gendata.in > job.out
```

Start it with a single MPI process. The MPI domain decomposition is disabled
in this version and the program stops when it is given more than one.

A calculation is controlled by one input file of `keyword value` lines
(`Gendata.in` in the examples).

See [](../docs/user-guide/input-file.md) for its format and [](../docs/user-guide/input-keys.md)
for every key.

## Quick start

```sh
cd examples/01_graphene_bands
./run.sh              # runs grabnes on Gendata.in, writes generate.bands
python3 plot.py       # writes graphene_bands.png
```

For a minimal calculation with reference output, see {ref}`example 01 <ex-01>`;
all supported calculations are listed in [](examples/index.md).

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
independently in Python.
[`tests/regression/README.md`](https://github.com/gjung-group/grabnes/blob/main/tests/regression/README.md)
describes every check; measured results for both compilers are in
[](../docs/development/cluster-build-and-validation.md).

## Running without installing: Google Colab

```{warning}
This notebook was set up for the Twistronics 2023 workshop. It reads shared
Google Drive folders (Intel compilers, inputs) prepared for that workshop and
writes its inputs with the key names used before October 2026, so it does not
reflect the current code and its resources may no longer be accessible. For the
current code, build it as described above and run the [](examples/index.md).
```

A build of the code was made available on Google Colab for the workshop,
where the code can be driven from a Python Jupyter notebook:

{button}`Twistronics 2023 workshop notebook (Google Colab) <https://colab.research.google.com/drive/12_ofrjJi4pr6U9j-ER1WocYJkqxEffs-?usp=sharing>`

The same notebook is shown, cell by cell, in the [](tutorials/index.md). For
production work we recommend compiling the code on your own computational
resources to leverage the OMP capabilities of the code.

## Related packages

Several other real-space tight-binding packages exist, each with its own
strengths, so one may be more suitable than another for the calculation of
specific observables. In alphabetical order:

- [Kwant](https://kwant-project.org/)
- [PyBinding](https://docs.pybinding.site/en/stable/)
- [PythTB](https://www.physics.rutgers.edu/pythtb/)
- [Quantum-Kite](https://quantum-kite.com/)
- [TBPLaS](http://www.tbplas.net/)
- [Twister](https://github.com/qtm-iisc/Twister)
