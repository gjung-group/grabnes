# Restructuring of `HamHopping` and `HamOnSite`

State and plan of the work on the branch `refactor/ham-structure`. The aim is
the convention of the original code: routines that fit on a few screens, each
with one task.

## How a change is verified

- The model survey (`tests/regression/model_survey.py`, about 430 cases) stores
  a fingerprint of the complete Hamiltonian of every case; a restructuring step
  must leave every fingerprint unchanged, with the checked gfortran build and
  with the Intel build.
- Coverage of the tests, measured with a `--coverage` build of gfortran on the
  survey and the physics checks (October 2026):

  | File | Executable lines run |
  | --- | --- |
  | `ham.F90` | 82 % (`HamHopping` 87 %, `HamOnSite` 89 %) |
  | `atoms.F90` | 74 % |
  | `neigh.F90` | 21 % (the legacy search routines are never run) |
  | `diag.F90` | 14 % (spectral functions, Chern numbers and `DiagH0TAPW_withBlockH` are not run at all) |

  The parts of `ham.F90` that no test reaches: `diagotBG`, the
  `FrankMagneticField` branch, `manyBubbles` / `bigBubble`, `TBG.DiagPRB`,
  `BNBN.Diag`. They are moved, where they are moved, without any change.
  **`diag.F90` needs tests before it can be restructured.**
- Four production inputs (graphene/hBN at zero field and in a field, a
  four-layer cell, a strained cell) are rerun and compared byte for byte.

## Step 1 (done): self-contained blocks become internal procedures

`tools/maintenance/outline_block.py` moves a block that is not inside a loop
into an internal procedure (`contains`) of the same routine. The statements are
not changed: an internal procedure uses the variables of its routine, so
nothing is passed, renamed or declared again. The tool refuses blocks with
`return`, `goto` or statement labels, blocks inside loops, and blocks inside an
OpenMP region that is opened outside them.

- `HamOnSite`: the body went from 2575 to about 490 lines; 15 procedures
  (`HamOnSiteSpeciesEnergies`, `HamOnSiteMoirePotential`,
  `HamOnSiteLayerPotential`, `HamOnSitePNP`, `HamOnSiteSquareFunction`,
  `HamOnSiteTrigonalCDW`, `HamOnSiteMoireCDW`, ...). What is left reads the
  switches and calls one procedure per term. Not moved: the `Bubbles` block
  (it contains a statement label).
- `HamHopping`: 11 procedures (`HamHopShellsFromRigid`,
  `HamHopF2G2sParameters`, `HamHopInterlayerParameters`, `HamHopSecondMoire`,
  `HamHopRandomStrain`, `HamHopRealStrain`, `HamHopPeriodicStrain`,
  `HamHopRealisticBubbles`, `HamHopMagneticField`, `HamHopHaldane`,
  `HamHopWriteTables`).

`HamOnSiteMoirePotential` (1300 lines) and `HamHopRealStrain` (316 lines) are
still long; they can be split further in the same way, by the model they
serve.

## Step 2 (not started): the main loop of `HamHopping`

The loop over atoms and neighbors is one OpenMP loop of 4500 lines:

```
do i = 1, nAt                              ! one OpenMP parallel loop
   bond terms of the moire models that depend on the atom only  (650 lines)
   do j = 1, Nneigh(i)
      if (same layer) then
         do ilvl = 1, tbnn                 ! intralayer models, by shell  (2255 lines)
      else
         interlayer models, by TypeOfBL    (1529 lines)
```

Its branches cannot be moved into internal procedures: inside a parallel
region the private variables of the routine would be shared in the procedure.
They need procedures with argument lists, and the numbers show why this is a
design task and not a mechanical one:

- `HamHopping` declares 901 local variables; the main loop uses 646 of them.
- 287 are only read in the loop: the parameters of the models, set before it.
- 359 are assigned in the loop, 195 of which are in the `PRIVATE` clauses.
  The remaining ones are mostly model constants that are assigned inside the
  loop (`p1a0AA = 7.639_dp`, ...): every thread writes the same value, which
  works but is a data race in the strict sense and is repeated for every bond.

Order of work proposed:

1. **Parameters out of the loop.** Move the assignments of constants that do
   not depend on the bond to the parameter procedures that run before the loop
   (`HamHopInterlayerParameters` and one procedure per intralayer model).
   This removes the shared assignments from the parallel loop and is verified
   by the survey, one model at a time.
2. **One derived type per model family** for the parameters (intralayer F2G2,
   Koshino, Kaxiras, Srivani, the graphene/hBN distance laws, ...), filled by
   those procedures. The loop then reads `par%...` instead of 287 separate
   variables.
3. **One function per model** that returns the hopping of a bond from the bond
   geometry and the parameter type: `HamIntralayerF2G2(par, ilvl, species, d)`,
   `HamInterlayerKoshino(par, d, dz)`, ... These are module procedures with a
   short argument list and no access to the variables of `HamHopping`, so they
   can be tested alone.
4. The loop keeps the selection of the model and the assembly of `hopp`.

Each of these steps changes statements, unlike step 1, so each is done for one
model, checked with the survey and the production inputs, and committed
before the next.

## Not planned here

- `diag.F90`: tests first (TAPW with spin-orbit terms, spectral functions,
  Chern numbers have none).
- `AtomsPos` (1192 lines) and `NeighList` (504 lines) can follow step 1 as they
  are; the legacy search routines of `neigh.F90` are candidates for removal
  rather than restructuring.
