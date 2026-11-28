(examples)=
# Examples

These examples are ordered from a fast installation check to more complete
moire calculations. They are the supported, reproducible calculations of the
[GRABNES code repository](https://github.com/gjung-group/grabnes); every input
file is also downloadable below.

| Example | Purpose | Status |
| --- | --- | --- |
| {ref}`01 <ex-01>` | Pristine graphene band structure | Ready |
| {ref}`02 <ex-02>` | Pristine graphene density of states | Ready |
| {ref}`03 <ex-03>` | 76-atom commensurate twisted-bilayer bands | Ready |
| {ref}`04 <ex-04>` | 76-atom commensurate twisted-bilayer DOS | Ready |
| {ref}`05 <ex-05>` | Aligned graphene on hBN (effective moire model): bands of one valley | Ready, 20 s |
| {ref}`06 <ex-06>` | Same system in a magnetic field: density of states against field, 0 to 26 T | Ready, 4 min |
| {ref}`07 <ex-07>` | Graphene on hBN at four twist angles: density of states | Ready, 5 min |
| {ref}`08 <ex-08>` | Twisted bilayer graphene at 1.08 degrees: flat bands | Ready, 2 min |
| {ref}`09 <ex-09>` | Alternating-twist trilayer at 1.54 degrees, structure read from files | Ready, 3 min |
| {ref}`10 <ex-10>` | Twisted bilayer at 6 degrees: spectral function unfolded on the graphene lattice, per layer | Ready, 10 s |
| {ref}`11 <ex-11>` | Alternating-twist trilayer at 6 degrees: the same, per layer | Ready, 30 s |

## Quick start

First compile GRABNES so that the executable is available at
`lanczosKuboCode/bin/grabnes` (see [](../installation.md)). Then run:

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
Examples 05 to 09 are moire systems at realistic sizes. Examples 10 and 11 use small moire
cells, because their spectral functions need a direct diagonalization at every k-point.

Examples 01 to 04, 10 and 11 use exact diagonalization. Examples 05, 08 and 09 project the
Hamiltonian on plane waves around one graphene valley, which makes cells of
several thousand atoms a matter of minutes; each of them states how the
reduction was checked against exact diagonalization. Examples 06 and 07 use
the recursion method on cells of 0.4 and 3.5 million atoms; example 07 needs
1.5 GB of memory. The times in the table are for a build with gfortran on
four cores.

All examples must be run with a single MPI process. To check every example
against its reference data in one step, without writing into the example
directory, use:

```sh
./tests/regression/run_examples.sh
```

```{seealso}
The example folders in the
[code repository](https://github.com/gjung-group/grabnes/tree/main/examples)
also contain reference data; the input files and scripts are mirrored under
the download links of each example page.
```

:::{toc}
:context: children
:::
