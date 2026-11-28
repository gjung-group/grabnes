# Building moiré systems

## The input file

Every GRABNES run is driven by a single input file, conventionally named
`Gendata.in`. The file is *not* read sequentially by the code, i.e. the order
in which the input flags are provided does not matter. The usual format is
the name of the flag followed by its value (a numerical value, a boolean
`.false.` or `.true.`, a string, or a code block).

```{important}
The parser requires negative values to be given as `0.0-1.0` instead of
`-1.0`.
```

The two flags that define what is built are `TypeOfSystem`, which selects the
structure-generating method, and the cell-size flags `CellSize` and
`Supercell`.

| Method | Use case |
|---|---|
| `Graphene` | Build graphene with the built-in methods. Useful for quick rigid structures. |
| `TwistedBilayerBasedOnMoireCell` | Build a commensurate twisted bilayer from the four Hermann indices. |
| `ReadXYZ` | Read atomic positions and sublattice assignment from external files. Required for relaxed structures, e.g. from LAMMPS. |

`CellSize` is the number of repetitions of the (graphene) unit cell, useful
when building a larger moiré unit cell, while `Supercell` repeats the cell
defined by `CellSize` that many times in the $x$ and $y$ directions. The
latter makes it easy to repeat moiré cells that are larger than the initial
graphene cell, which is what you need when sampling a DOS with Lanczos
recursion.

## Commensurate cells and Hermann indices

Twisted bilayer graphene has garnered a lot of attention in recent years due
to the interesting properties it displays. Twisting the two graphene layers at
a specific angle makes the hybridization of the two Dirac cones lead to the
formation of a flat band, which can give rise to superconductivity and Mott
insulator signatures.

To simulate this type of moiré system one needs to build a commensurate cell
for each of the twist angles of interest. Unlike the untwisted systems, these
commensurate cells can be quite large. As a first approximation, the smaller
the twist angle, the larger the moiré unit cell. Only a discrete set of twist
angles exists, but the larger the unit cell, the more angles become
accessible.

A versatile way to describe a moiré system uses four indices, as introduced in
[this reference](https://iopscience.iop.org/article/10.1088/0953-8984/24/31/314210).
The commensurate moiré superlattices are fully defined by four integers $i$,
$j$, $i^\prime$ and $j^\prime$, where $i$ and $j$ refer to the bottom
graphene layer and $i^\prime$ and $j^\prime$ to the other layer. The two
lattice vectors of the commensurate supercell $\textbf{r}_1$, $\textbf{r}_2$
relate to the lattice vectors of the bottom layer $\textbf{a}_1$,
$\textbf{a}_2$ and of the top layer $\textbf{a}^\prime_1$,
$\textbf{a}^\prime_2$ through

$$
\begin{pmatrix}
\textbf{r}_{1} \\
\textbf{r}_{2}
\end{pmatrix}
= \textbf{M}^{\prime} \cdot
\begin{pmatrix}
\textbf{a}_{1} \\
\textbf{a}_{2}
\end{pmatrix}
= \textbf{M} \cdot
\begin{pmatrix}
\textbf{a}^{\prime}_{1} \\
\textbf{a}^{\prime}_{2}
\end{pmatrix}
$$

with the transformation matrices

$$
\textbf{M} =
\begin{pmatrix}
      i      &     j        \\
     -j      &     i + j
\end{pmatrix},
\quad
\textbf{M}^{\prime} =
\begin{pmatrix}
     i^\prime      & j^\prime       \\
    -j^\prime     & i^\prime + j^\prime
\end{pmatrix}.
$$

Defining $a_{G_1}$ and $a_{G_2}$ as the lattice constants of the two
constituent layers, the lattice mismatch ratio $\alpha$ follows from the four
indices as

$$
\alpha = \frac{a_{G_1}}{a_{G_2}} =\sqrt{\frac{(i^{\prime 2}+j^{\prime 2}+i^{\prime}j^{\prime})}{i^2+i^2+ij}}
$$

and the rotation angle $\theta$ as

$$
\cos(\theta) = \frac{1}{2 \alpha g}\left[ 2 i^\prime i + 2 j^\prime j + i^\prime j + j^\prime i \right]
$$

where $g=i^2+j^2+ij$. In this formulation $a_{G_1} = a_{G_2}$ is set to avoid
internal strain.

```{warning}
The denominator of $\alpha$ is reproduced above exactly as it appears in the
workshop notebook. Taken together with the definition of $g$ below it is
inconsistent for $i \neq j$, so it should probably read $i^2 + j^2 + ij$.
This has not been verified against the implementation, so the original
expression is kept unmodified — please confirm with the code maintainers
before changing it here.
```

A table with the relevant indices for common angles is given in
[this reference](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.115410).
The indices are supplied to the code through the `MoireCellParameters` flag,
e.g. `MoireCellParameters 31 30 30 31` for the first magic angle.

## Lattice relaxation and the LAMMPS interface

If you are familiar with the field of twisted bilayer graphene, you are
probably aware that relaxation effects play an important role in the accurate
assessment of electronic structure observables. Hence most calculations
require an energy minimization step, which can be routinely achieved with most
molecular dynamics packages.

GRABNES provides the means to interface with LAMMPS: LAMMPS reads in the
commensurate cell produced by this package, and this package reads in the
relaxed positions from the LAMMPS minimization output. Concretely, the
structure-generating step writes output files that can be used as input for a
LAMMPS calculation or directly as input for GRABNES. To go the relaxation
route:

1. Build the commensurate cell with the built-in moiré methods.
2. Convert the output to the format expected by LAMMPS and minimize the
   energy using registry-dependent force fields (EXX-RPA-informed, with
   REBO2 or GAP20 intralayer terms).
3. Feed the relaxed coordinates back into GRABNES with
   `TypeOfSystem ReadXYZ`, providing in the working directory the
   `generate.xyz` file with the atomic positions and the
   `sublatticesSorted.dat` file stating which sublattice each atom belongs
   to.

```{tip}
To find out what a flag does, grep the source code for it:

```bash
grep -r "SuperCell ." Src/
```

This is the quickest way to discover all occurrences of a flag and the
routines that consume it.
```
