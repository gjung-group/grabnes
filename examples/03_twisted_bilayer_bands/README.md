# 03: Twisted-bilayer graphene bands

This example constructs the `(m,n) = (3,2)` commensurate twisted-bilayer
graphene cell (76 atoms, approximately 13.17 degrees) and evaluates its bands
along `K - Gamma - M - K'` with the Koshino interlayer hopping model.

## What the reference data represent

This example is a regression test of one specific Hamiltonian; it is not
independent evidence that this Hamiltonian is the best description of twisted
bilayer graphene.

- Intralayer: nearest-neighbor hopping only, -2.9888 eV (the default
  `TB.Hopping` for a lattice parameter of 2.46 A). `Gendata.in` requests this
  explicitly with `Neigh.CutAtNN3` and by setting the second- and
  third-neighbor F2G2 hoppings to zero.
- Interlayer: two-center form with `vpppi0` = 3.5 eV, `vppsigma0` = 0.48 eV,
  decay length 0.184 a, interlayer distance 3.34 A, for pairs up to 9.69 A
  apart in the plane.

GRABNES uses the F2G2 intralayer model by default (remove those three lines
and set `TB.NeighLevels 5`). It moves the Dirac point to about -0.33 eV and
the band edges by more than 2 eV, and has no reference data yet. The active
terms, an independent reconstruction accurate to 1e-4 eV, and the assessment
of the F2G2 variant are documented in
[`docs/development/cluster-build-and-validation.md`](../../docs/development/cluster-build-and-validation.md).

The relatively large twist angle is intentional: it keeps this tutorial
calculation small enough for a laptop. Near-magic-angle cells contain thousands
of atoms and are production calculations, not installation checks.

## Run

```sh
./run.sh
python3 plot.py
```

The calculation writes `generate.bands` and `job.out`; the plotting step writes
`twisted_bilayer_bands.png`. If no calculated file exists, the plotter uses
`reference/bands.dat`.
