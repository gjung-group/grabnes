# 03: Twisted-bilayer graphene bands

This example constructs the `(m,n) = (3,2)` commensurate twisted-bilayer
graphene cell (76 atoms, approximately 13.17 degrees) and evaluates its bands
along `K - Gamma - M - K'` with the Koshino interlayer hopping model.

The intralayer model is deliberately the simplest one, nearest-neighbor
hopping only. `Gendata.in` requests it explicitly with `Neigh.CutAtNN3` and by
setting the second- and third-neighbor F2G2 hoppings to zero; this is the
model the reference data were computed with. GRABNES uses the F2G2 intralayer
model by default (remove those three lines and set `TB.NeighLevels 5`), which
shifts the bands by up to about 2 eV and has no reference data yet.

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
