# GRABNES Fortran code

GRABNES (GRAphene and Boron Nitride Electronic Structure) is a Fortran code
for electronic-structure and quantum-transport calculations in layered
two-dimensional materials.

## Requirements

- GNU Fortran (tested: 12.2) or Intel Fortran (tested: `ifort` 2021.6)
- MPI with a Fortran compiler wrapper (`mpif90`, `mpiifort`, ...)
- BLAS and LAPACK
- ARPACK
- GNU Make

## Build

Create a local build configuration from the portable template and compile:

```sh
cp make.sys.example make.sys
make
```

The executable is written to `bin/grabnes`. The local `make.sys`, `build/`,
and `bin/` paths are ignored by Git; keep machine-specific paths in `make.sys`.
An Intel Fortran configuration that has been tested is
`../grabnes_testrun/config/intel.make.sys`. The files under `Sys/` are
historical site configurations with absolute paths and are kept only as
starting points.

The build does not depend on the current directory and can be placed outside
the source tree:

```sh
make -C /path/to/grabnes/lanczosKuboCode \
     MAKE_SYS=/path/to/make.sys BUILD_DIR=/scratch/build BIN_DIR=/scratch/bin
```

With GNU Fortran the OpenMP directives are not compiled, because the compiler
rejects several of them; see
[`../docs/development/cluster-build-and-validation.md`](../docs/development/cluster-build-and-validation.md).

## Run

GRABNES accepts the input filename as its first argument:

```sh
bin/grabnes Gendata.in > job.out
```

Start it with a single MPI process. The MPI domain decomposition is disabled
in this version and the program stops when it is given more than one.

For a minimal calculation with reference output, see
[`../examples/01_graphene_bands`](../examples/01_graphene_bands/). To build
and check all examples in one step, run `../grabnes_testrun/run_examples.sh`.

## Clean

```sh
make clean
```

The code includes band-structure, density-of-states, Kubo transport, TAPW,
Berry-curvature, Chern-number, and semiclassical-orbit functionality. The
public examples document supported workflows as they are validated.
