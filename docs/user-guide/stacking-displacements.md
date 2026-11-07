# Stacking displacements (moiré terms of relaxed structures)

Some moiré terms change with the local stacking between two layers: the graphene/hBN terms (`GBN.OffDiag`, the
on-site harmonic terms) and the twisted-bilayer terms (`TBG.OffDiag`, `TBG.Diag`). For each atom they need the
**stacking displacement** d: the in-plane offset between the atom and the matching atom of the adjacent layer.

## Convention

- d = 0 at AA stacking. For graphene/hBN that is carbon A over boron and carbon B over nitrogen; for two graphene
  layers it is A over A and B over B.
- The sign is "graphene minus hBN" for a graphene/hBN interface, and "lower layer minus upper layer" for two layers
  of the same material. This matches the analytic form d(r) = ((1+ε)R(θ) − 1) r used when no displacements are given
  (ε < 0 for graphene on hBN).
- The graphene layer of `Stack.EncapsulatedThreeLayers` has one displacement for each of its two hBN interfaces.

## Three ways to give d

1. **Analytic (default).** d(r) = ((1+ε)R(θ) − 1) r, from `Moire.LatticePercent` (ε) and `GBN.Angle` / `TBG.Angle`
   (θ). This assumes a rigid cell with an AA point at the origin of the coordinates.
2. **From the structure (recommended for relaxed cells).** `Structure.DisplacementsFromPositions .true.`, together
   with `GBN.UseDisplacementFile .true.` or `TBG.UseDisplacementFile .true.` and `Structure.ReadRigidXYZ .true.`.
   grabnes computes d itself from the rigid structure (`Structure.RigidPositions`, default `generateInit.xyz`) and
   the relaxed structure (the usual structure file, `generate.xyz`):
   - the partner of each atom is the atom of the adjacent layer with the matching sublattice, nearest in plane in
     the rigid structure;
   - d is the in-plane vector between the two in the relaxed structure.

   No ε, θ or origin is involved, so the result does not depend on how the cell is oriented or where its origin is.
   No separate file is needed.
3. **From a file.** `GBN.UseDisplacementFile` / `TBG.UseDisplacementFile` with `Structure.DisplacementsFile`
   (default `displacements.txt`): one line per atom (`dx dy dz`; six columns for the graphene layer of
   `Stack.EncapsulatedThreeLayers`), in the convention above. `Stack.GlobalTwist*` rotate the file values into the
   frame of the structure.

## Conditions of option 2

- `Stack.GlobalTwist`, `GlobalTwist2`, `GlobalPhiL1`, `L2a`, `L2b`, `L3` must be 0: the displacements are already
  in the frame of the structure.
- The partner must be among the atom's interlayer neighbours. Its in-plane distance is at most about 1.45 Å, so
  `Neigh.LayerDistFactor` must be at least about 1.1. The run stops with a message when a partner is missing.
- An atom between two layers of its own material (a middle layer of multilayer graphene) is refused, because it is
  not defined which interface its term belongs to.
- `SuperCell` must be 1 (as for `Structure.ReadRigidXYZ`).

The run prints the number of atoms and the largest |d| found.

## Note on displacement files made before November 2026

Before November 2026 the twisted-bilayer terms used the opposite sign for file values, while the graphene/hBN terms
used the convention above. Both now use the same convention, so twisted-bilayer runs with a displacement file made
before then change the sign of their moiré terms.
