# Commensurate cells from four indices (`Structure.MoireCellParameters`)

The cells of `TwistedBilayerBasedOnMoireCell`, `Graphene_Over_BN`, `BilayerGraphene_Over_BN`,
`MoireEncapsulatedBilayerBasedOnMoireCell`, `TrilayerBasedOnMoireCell` and of `Graphene` / `BoronNitride` with
`Structure.BasedOnMoireCellParameters .true.` are built from Hermann's four integers
(K. Hermann, J. Phys.: Condens. Matter 24, 314210 (2012)):

```
Structure.MoireCellParameters m1 m2 m3 m4
```

## Definition

Each layer has the lattice vectors a1 = a (1, 0) and a2 = a (1/2, sqrt(3)/2) of its own lattice constant a,
turned with the layer. The supercell vector is

    A1 = m1 a1 + m2 a2   (first lattice)   =   m3 a1' + m4 a2'   (second lattice),

and A2 is A1 turned by +60 degrees. In complex notation, with tau = exp(i pi/3), the lattice point m a1 + n a2 is
a (m + n tau): the pair is the two Eisenstein integers alpha = m1 + m2 tau and beta = m3 + m4 tau. Then

- the second lattice is the first one scaled by |alpha/beta| and turned by theta = arg(alpha/beta);
- the cell holds N(alpha) = m1^2 + m1 m2 + m2^2 unit cells of the first lattice and N(beta) of the second
  (2 N atoms per honeycomb layer).

**Any** four integers with (m1, m2) and (m3, m4) both nonzero define a valid cell: signs and zeros are allowed,
there is no restriction to the form (a, b, b, a), and the lattice ratio follows from the indices (the second lattice
constant is set to aG |alpha/beta|; `Structure.LatticeParameterBN` is overwritten). A zero pair is refused.
Negative integers are written `0-1` in `Gendata.in` (a bare `-1` aborts the input parser).

Equivalent index sets: (u alpha, u beta) for a unit u (the six powers of tau: (1,0), (0,1), (-1,1), (-1,0), (0,-1),
(1,-1) multiply as Eisenstein integers) is the same cell with A1 turned by a multiple of 60 degrees. Exchanging the
two pairs gives theta -> -theta, and so does the mirror (m1 + m2, -m2, m3 + m4, -m4).

## What the solver builds and prints

The pair of the larger norm is built first (layer 1, upper); for equal norms (twisted bilayer graphene) the pairs
are exchanged and (m3, m4) is layer 1. Layer 1 is at angle 0, the other layer at the signed angle that the indices
define. `job.out` has, besides the former lines (`cell: Angle:` is the unsigned angle, as before):

```
cell: Hermann cell: the (m3, m4) lattice is the (m1, m2) lattice turned by -38.213211 deg (signed, counterclockwise)
cell: Hermann cell: layer 1 (upper) is built from (m3, m4), layer 2 from (m1, m2), turned by 38.213211 deg ...
cell: Hermann cell: NOT primitive, it holds 3 primitive cells (built as given). Primitive cell with the same
      twist: 3 -1 2 1
```

## Primitive and non-primitive cells

The points common to both lattices form the lattice (alpha / gcd(alpha, beta)) Z[tau], where gcd is the greatest
common divisor in the Eisenstein integers (a Euclidean ring). Hence

- the cell is **primitive** iff gcd(alpha, beta) is a unit;
- otherwise it holds N(gcd) primitive cells, and (alpha/g, beta/g) is the primitive cell of the same geometry.

Removing a common integer factor of the four indices is not enough: (4, 1, 1, 4) has none and holds 3 primitive
cells; (5, 6, -1, 10) has none and holds 7 (a sqrt(7) x sqrt(7) supercell). For twisted bilayer graphene written as
(a, b, b, a) with coprime a > b, the cell is primitive iff a - b is not a multiple of 3; otherwise its primitive cell
with the **same** twist is

    ( (2a+b)/3, (b-a)/3, (a+2b)/3, (a-b)/3 ),

with 3 times fewer atoms. (The same expression with the two pairs exchanged, ((a+2b)/3, (a-b)/3, (2a+b)/3, (b-a)/3),
has the same |theta| and atom count but the opposite twist: it is the mirror structure.)

| angle (deg) | (a, b, b, a) | atoms | primitive, same twist | atoms |
|---|---|---:|---|---:|
| 21.7868 | 2 1 1 2 | 28 | 2 1 1 2 | 28 |
| 38.2132 | 4 1 1 4 | 84 | 3 -1 2 1 | 28 |
| 27.7958 | 5 2 2 5 | 156 | 4 -1 3 1 | 52 |
| 46.8264 | 7 1 1 7 | 228 | 5 -2 3 2 | 76 |
| 17.8966 | 7 4 4 7 | 372 | 6 -1 5 1 | 124 |

**The solver never reduces a cell.** A non-primitive cell is built exactly as requested (the job.out line above
says so and gives the primitive indices): such supercells are needed for states that break the translation symmetry
of the primitive cell (charge order, magnetism, Hartree-Fock). Products (g alpha, g beta) for any Eisenstein
integer g give the C6-symmetric supercells of the primitive cell: g = 2 the 2 x 2, g = 1 + tau the sqrt(3) x sqrt(3)
R30, g = 2 + tau the sqrt(7) x sqrt(7) cell. Rectangular or other non-hexagonal supercells cannot be written with four
indices.

## Completeness for twisted bilayer graphene

For two identical lattices every commensurate angle is reached, and with its primitive cell, once signed indices
are allowed: up to 12000 atoms there are 828 angles in (0, 60) deg, and the signed-index search gives the same 828
angles and the same primitive atom counts as the formulas of Mele (Phys. Rev. B 81, 161405 (2010)) and Lopes dos
Santos et al. (Phys. Rev. B 86, 155449 (2012)), two-way. With non-negative indices only, the former GenMoire search
reached 544 of them within 12000 atoms, 409 with their primitive cell (and 112 with its default filter
N(alpha - beta) <= 45). `tools/structure/hermann.py --selftest` checks all of this.

## Finding a cell: `usefulGeneralScripts/twistedBilayers/GenMoire.f90`

Two different questions:

- **Exact geometry given by indices** (`GenMoire cell m1 m2 m3 m4 [reduce]`): says whether the cell is primitive
  and gives the primitive cell of the same twist and lattice ratio (used instead of the given one with `reduce`).
- **Best commensurate approximation** of a target angle and lattice ratio
  (`GenMoire search phideg aBN degtol lattol [nCellMax nMoire1 nMoire2 [signed]]`, or the values at the top of the
  file): every primitive cell with at most nCellMax cells of the first lattice inside the tolerances, signed indices,
  each geometry once (beta in the sector m3 > 0, m4 >= 0 removes the units; non-primitive cells are skipped, since
  their primitive cell is met too). The selection criterion `angle`, `delta` or `both` is as before; equal errors go
  to the smaller cell, and the smallest cell within the tolerances is printed as well. By default only |angle| is
  compared, as before; with `signed` phideg is the signed rotation of the second layer. The program prints the
  `Gendata.in` line (negative values as `0-n`).

`tools/structure/hermann.py` does the same in Python (`describe`, `search_hermann`, `search_tbg_hermann`) with exact
integer arithmetic.

## Conventions that did not change

- Rotation centre: both layers are built from the cell origin, a hexagon centre of both honeycomb lattices (AA
  stacking of the hexagon centres at the origin, no relative translation unless `Structure.TwistedBLAddShift`). Since
  a turn by 60 degrees about a hexagon centre maps a honeycomb layer onto itself, theta and theta + 60 deg give the
  same structure, and theta and 60 deg - theta are mirror images. Mele's SE-even / SE-odd distinction (rotation about
  an atom) is therefore not reproduced by the stacking used here: that is a separate question from the completeness of
  the indices and is not changed by them.
- `Structure.TwistedBLAddShift` with `Structure.BernalShift` / `BridgeShift` shifts the second layer by
  (aG/3)/(m1 aG) of A1 + A2 (and `BilayerShiftFactor` along A1 + A2): this is the Bernal or bridge shift only for
  cells (n, 0, ...). For other indices the added shift is something else; this is unchanged.
- Terms that take the angle of the cell (the moire effective model with `Structure.BasedOnMoireCellParameters`,
  `Spectral.RotateReferenceSystem`, the analytic twisted-bilayer terms) use the unsigned angle, i.e. the sense of the
  (a, b, b, a) cells with a > b. When the built cell has the opposite sense, the run prints a warning with the mirror
  indices; the structure itself is always built as the indices say.
