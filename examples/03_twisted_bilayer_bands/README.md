# 03: Twisted-bilayer graphene bands

This example constructs the `(m,n) = (3,2)` commensurate twisted-bilayer
graphene cell (76 atoms, approximately 13.17 degrees) and evaluates its bands
along `K - Gamma - M - K'` with the Koshino interlayer hopping model.

## What the reference data represent

This example is a regression test of one specific Hamiltonian; it is not
independent evidence that this Hamiltonian is the best description of twisted
bilayer graphene.

- Intralayer: nearest-neighbor hopping only (`TB.NeighLevels 1`), -2.9888 eV,
  the default `TB.Hopping` for a lattice parameter of 2.46 A.
- Interlayer: two-center form of Moon and Koshino, Phys. Rev. B 85, 195458
  (2012), with `vppsigma0` = 0.48 eV, decay length 0.184 a, interlayer distance
  3.34 A, for pairs up to 9.69 A apart in the plane. `vpppi0` is 3.5 eV rather
  than their 2.7 eV. This is a deliberate calibration aimed at a more realistic
  Dirac velocity; it acts on the velocity when the intralayer hopping also
  follows the two-center form (`KoshinoIntralayer .true.`), and has a
  negligible effect in this example, whose velocity is set by `TB.Hopping`.

`TB.NeighLevels` is the number of intralayer neighbor shells. With
`TB.NeighLevels 5` GRABNES uses its F2G2 intralayer model, which moves the
Dirac point to about -0.33 eV and the band edges by more than 2 eV; it is
checked against an independent implementation but has no reference data here.

The reference data were regenerated in October 2026 after the neighbor search
was corrected (changes up to 1e-4 eV). The solver reproduces an independently
built Hamiltonian for this example to 4e-16 eV per matrix element; details are
in
[`docs/development/cluster-build-and-validation.md`](../../docs/development/cluster-build-and-validation.md).

## Run

```sh
./run.sh
python3 plot.py
```

The calculation writes `generate.bands` and `job.out`; the plotting step writes
`twisted_bilayer_bands.png`. If no calculated file exists, the plotter uses
`reference/bands.dat`.
