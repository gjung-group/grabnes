# GRABNES examples

These examples are ordered from a fast installation check to more complete
moire calculations.

| Example | Purpose | Status |
| --- | --- | --- |
| [`01_graphene_bands`](01_graphene_bands/) | Pristine graphene band structure | Ready |
| [`02_graphene_dos`](02_graphene_dos/) | Pristine graphene density of states | Ready |
| [`03_twisted_bilayer_bands`](03_twisted_bilayer_bands/) | 76-atom commensurate twisted-bilayer bands | Ready |
| [`04_twisted_bilayer_dos`](04_twisted_bilayer_dos/) | 76-atom commensurate twisted-bilayer DOS | Ready |

## Quick start

First compile GRABNES so that the executable is available at
`lanczosKuboCode/bin/grabnes` (see the repository README). Then run:

```sh
cd examples/01_graphene_bands
./run.sh
python3 plot.py
```

Set `GRABNES_BIN` if the executable is located elsewhere:

```sh
GRABNES_BIN=/path/to/grabnes ./run.sh
```

Each example contains its input, a short explanation, a plotting script, and
reference data from a known-good calculation. Examples 01 and 02 run in a few
seconds; examples 03 and 04 use a small, 76-atom commensurate twisted bilayer
so they remain practical laptop checks rather than production calculations.

All four examples use exact diagonalization and must be run with a single MPI
process. To check every example against its reference data in one step,
without writing into this directory, use:

```sh
./tests/regression/run_examples.sh
```
