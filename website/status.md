# Status and limitations

This page is quoted from the
[README](https://github.com/gjung-group/grabnes/blob/main/README.md) of the code
repository.

```{important}
GRABNES is under active development. There is no numbered release yet, so
please cite the commit you used (see [](citing.md)). The examples in
[](examples/index.md) are tested on every change and reproduce their reference
data; other functionality is research code of varying maturity (see
[What is supported](#what-is-supported)).
```

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

[](../docs/development/functionality-status.md) states for each capability what
evidence exists.

(what-is-supported)=
## What is supported

| Level | Meaning | Applies to |
| --- | --- | --- |
| Supported | Runs, is documented, and is reproduced by the test suite | The nine examples above (06 to 09 by `run_long_examples.sh`) |
| Independently validated | Additionally checked against an analytic result or a second implementation | Hamiltonian assembly and bands of graphene, hBN monolayer, and twisted bilayers; Kubo DOS of graphene; TAPW bands of examples 05, 08 and 09 against exact diagonalization |
| Research / experimental | Present in the source without sufficient validation for general use | Everything else, including transport, spin-orbit, Berry curvature, atomistic graphene/hBN, and the models of the examples used with other inputs |

"The nine examples above" are examples 01 to 09 of the README; they are
described, with examples 10 and 11, in [](examples/index.md).

## Known limitations

- One MPI process only; a run with more is refused. Parallelism is by OpenMP.
- The Kubo solver is stochastic: results vary from run to run unless
  `Kubo.SetSeed .true.` and `Kubo.SeedValue` are given.
- The bond term of the effective moire-potential model of graphene on hBN
  (`Moire.OffDiag`) gives a hopping table that is not exactly Hermitian; the
  solver reports the size of the asymmetry. `Moire.OffDiagMidpoint .true.`
  evaluates the term at the bond midpoint and makes it Hermitian. See the
  [validation document](../docs/development/cluster-build-and-validation.md).
- Many input keys exist only for research workflows. Nearly all have a one-line
  description in [](../docs/user-guide/input-keys.md), but only the keys of the
  examples are tested; the example inputs and READMEs are the reference for
  the supported ones.

## Reporting problems

Please report reproducible problems through the
[GitHub issue tracker](https://github.com/gjung-group/grabnes/issues) of the
repository, with the commit hash, compiler version, input file, and log.
