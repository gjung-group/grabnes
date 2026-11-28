---
edit_url: null
---

(api-ham-hamhopping)=
# `HamHopping`

```{index} single: HamHopping (subroutine in ham)
```

```fortran
subroutine HamHopping()
```

- **Module:** {ref}`ham <api-ham>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/ham.F90, line 5373`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5373)

Assemble hopping integrals (intra/inter-layer) into hopp(:,:).

Computes hopping amplitudes using selected bilayer/single-layer models, including twist angles, strain/bubble effects, and optional magnetic phases. Supports OpenMP for performance.

```{note}
Works with parameters set in HamInit and structures from neigh/atoms.
```

**Side effects**

- Updates hopp(:,:) in place
- May access/transform geometry data
- Applies magnetic field phase factors if enabled

## Input keys read (183)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Strain.ShellsFromRigidPositions`
  - GBNtwoLayers: take the shell of every intralayer bond from its in-plane length in a rigid reference structure (Strain.ShellsRigidFile, same atom order); the strain law keeps the 3-D length.
* - `Moire.Strain`
  - Here: the moire-strain hopping law is applied in HamHopping.
* - `MagField.FrankMagneticField`
  - Peierls phases from a bond-direction gauge for a rectangular graphene cell (HamHopMagneticField), with MagField.Integer flux quanta.
* - `MagField.Integer`
  - Here: number of flux quanta of MagField.FrankMagneticField.
* - `Interlayer.ChangeLatticeParameterForSrivaniModel`
  - BLSrivani: compute rbar with aG = 2.43898 Angstrom, the lattice constant of the fit.
* - `Intralayer.RemoveF2G2Flag`
  - Switch Intralayer.F2G2Model off (use the 8-shell parameter set).
* - `Intralayer.F2G2Model`
  - Here: F2G2 intralayer hoppings (otherwise the 8-shell parameter set).
* - `Stack.ThreeLayerShort`
  - Three-layer stacks: short-range interlayer model with hoppings tB1A2 (atoms on top of each other) and tA1A2 / tA1B2 (first ring).
* - `TBG.OffDiag`
  - Here: the tBG bond term is added in HamHopping.
* - `TBG.OffDiagPRB`
  - TBG.OffDiag: alternative parameter set (see HamHopping).
* - `Moire.LatticePercent`
  - Here: lattice mismatch of the moire bond terms.
* - `Moire.PotCab`
  - Here: amplitude of the off-diagonal moire term on graphene bonds \[eV\].
* - `Moire.PotCabBN`
  - Amplitude of the off-diagonal moire term on hBN bonds \[eV\].
* - `Interlayer.BilayerOneParameter`
  - Bernal bilayer with one interlayer hopping, Intralayer.BilayertAB1 (default 0.361 eV).
* - `Intralayer.BilayerThreeParameters`
  - Bernal bilayer with three interlayer hoppings tAB1, tAB3, tAB4.
* - `Moire.BfactorCab`
  - Moire.DistanceDependentEffectiveModel: decay B of Cab \[1/Angstrom\].
* - `Structure.InterlayerDistance`
  - Here: reference distance z0 of the distance-dependent Cab \[Angstrom\].
* - `Intralayer.BilayertAB1`
  - Bernal-bilayer hopping t\_AB1 (first set: one-parameter model) \[eV\] of the bilayer F2G2 / few-parameter models; the default depends on the parameter set read here.
* - `Intralayer.BilayertAB3`
  - Bernal-bilayer hopping t\_AB3 \[eV\] of the bilayer F2G2 / few-parameter models; the default depends on the parameter set read here.
* - `Intralayer.BilayertAB4`
  - Bernal-bilayer hopping t\_AB4 \[eV\] of the bilayer F2G2 / few-parameter models; the default depends on the parameter set read here.
* - `Intralayer.BilayertAA1`
  - Bernal-bilayer hopping t\_AA1 \[eV\] of the bilayer F2G2 / few-parameter models; the default depends on the parameter set read here.
* - `Intralayer.BilayertAA3`
  - Bernal-bilayer hopping t\_AA3 \[eV\] of the bilayer F2G2 model (stored as tAA2).
* - `Intralayer.BilayertBA0`
  - Bernal-bilayer hopping t\_BA0 (dimer) \[eV\] of the bilayer F2G2 / few-parameter models; the default depends on the parameter set read here.
* - `Intralayer.BilayertBA2`
  - Bernal-bilayer hopping t\_BA2 \[eV\] of the bilayer F2G2 / few-parameter models; the default depends on the parameter set read here.
* - `Intralayer.BilayertBA5`
  - Bernal-bilayer hopping t\_BA5 \[eV\] of the bilayer F2G2 / few-parameter models; the default depends on the parameter set read here.
* - `Trilayer.Gamma3`
  - Trilayer.FanZhang: gamma3 \[eV\]; also the former source of gamma1 when Trilayer.Gamma1 is absent (warning).
* - `Trilayer.Gamma1`
  - Trilayer.FanZhang: vertical dimer hopping gamma1 \[eV\].
* - `Trilayer.Gamma4`
  - Trilayer.FanZhang: gamma4 \[eV\].
* - `Trilayer.Gamma2`
  - Trilayer.FanZhang: gamma2 \[eV\].
* - `Moire.PotPhiab`
  - Here: phase of the off-diagonal moire bond term \[rad\].
* - `Moire.PotPhiabG`
  - Phase of the off-diagonal moire term on graphene bonds \[degrees\] (0 with TBG.OffDiag).
* - `Moire.PotPhiabBN`
  - Phase of the off-diagonal moire term on hBN bonds \[degrees\].
* - `Interlayer.TypeOfBL`
  - Here: selects the interlayer hopping branch of HamHopping.
* - `Interlayer.TypeOfSL`
  - Intralayer model name for single-layer systems (same names as Interlayer.TypeOfBL, e.g. HTC).
* - `Interlayer.AddExponentialDecayForDihedral`
  - Multiply the interlayer hopping of that branch by exp(-(d - Structure.InterlayerDistance) / Interlayer.BLdelta), d = 3-D bond length.
* - `Interlayer.MinDelta`
  - Read but not used at present.
* - `Interlayer.TwistedBLtbt`
  - Interlayer AA amplitude t\_bt of the analytic twisted-bilayer interlayer model \[eV\].
* - `Intralayer.KoshinoIntralayer`
  - Moon-Koshino intralayer hopping (Vpppi0 also used in-plane; its default becomes 2.7 eV).
* - `Interlayer.Vpppi0`
  - Two-centre V\_pp\_pi^0 \[eV\]: pi part of the interlayer hopping in every two-centre model, and the intralayer first-neighbour hopping with KoshinoIntralayer. Default 2.7 with KoshinoIntralayer, 3.5 otherwise.
* - `Interlayer.Vppsigma0`
  - Two-centre V\_pp\_sigma^0 of the interlayer hopping \[eV\].
* - `Interlayer.BLdelta`
  - Two-centre decay length delta of the interlayer hopping \[Angstrom\] (default 0.184 aG).
* - `Interlayer.CorrugatedInterlayerTwoCenter`
  - Two-centre interlayer hopping with local p\_z axes (corrugated sheets) instead of the global z.
* - `Structure.InterlayerDistance`
  - Here: interlayer distance of the Mayou model's q\_sigma \[Angstrom\].
* - `Intralayer.UseOldGrapheneF2G2`
  - Older graphene F2G2 parameter set (nonzero t4).
* - `Interlayer.DeactivateInterlayer`
  - Drop the interlayer bonds of the model branches that test it.
* - `Interlayer.DeactivateInterlayer12`
  - Drop the interlayer bonds between layers 1 and 2.
* - `Interlayer.DeactivateInterlayer23`
  - Drop the interlayer bonds between layers 2 and 3.
* - `Interlayer.DeactivateInterlayer34`
  - Drop the interlayer bonds between layers 3 and 4.
* - `Interlayer.DeactivateInterlayert2GBN1to2`
  - Stack.T2GBN: drop the interlayer bonds between layers 1 and 2.
* - `Interlayer.DeactivateInterlayert2GBN2to3`
  - Stack.T2GBN: drop the interlayer bonds between layers 2 and 3.
* - `Interlayer.DeactivateInterlayert2BG2to3`
  - Stack.T2BG: drop the interlayer bonds between layers 2 and 3.
* - `Interlayer.DeactivateInterlayert3BG1to2`
  - Stack.T3BG: drop the interlayer bonds between layers 1 and 2.
* - `Interlayer.DeactivateInterlayert3BG2to3`
  - Stack.T3BG: drop the interlayer bonds between layers 2 and 3.
* - `Interlayer.DeactivateInterlayert3BG3to4`
  - Stack.T3BG: drop the interlayer bonds between layers 3 and 4.
* - `Interlayer.DeactivateInterlayert3BG4to5`
  - Stack.T3BG: drop the interlayer bonds between layers 4 and 5.
* - `Interlayer.DeactivateInterlayert3BG5to6`
  - Stack.T3BG: drop the interlayer bonds between layers 5 and 6.
* - `Interlayer.DeactivateInterlayerBG`
  - Drop the interlayer bonds between layers 2 and 3 (bilayer-graphene branch).
* - `Interlayer.DeactivateInterlayerTwisted`
  - Drop the bonds between adjacent layers in the twisted-stack branch.
* - `Strain.PeriodicStrainPeriod`
  - Number of periods of the periodic strain / strained moire in the moire cell.
* - `Strain.StrainedMoireMaxDisplacement`
  - Strain.StrainedMoire: amplitude u\_max of the sinusoidal displacement \[Angstrom\].
* - `Strain.StrainedMoire`
  - Interlayer distances computed with a sinusoidal in-plane displacement u\_max sin(2 pi n x / L\_moire).
* - `Structure.SuperCell`
  - Here: the moire length is \|a1\| / SuperCell.
* - `Structure.TwistedBLAddShift`
  - Here: read for the interlayer branch of shifted bilayers.
* - `Intralayer.ForceBilayerF2G2Intralayer`
  - Here: the intralayer hoppings take the Bernal-bilayer F2G2 set below.
* - `Intralayer.Bilayert2KA`
  - F2G2 intralayer hopping t2 on sublattice A (second shell) \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.Bilayert2KB`
  - F2G2 intralayer hopping t2 on sublattice B (second shell) \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.Bilayert3K`
  - F2G2 intralayer hopping t3 (third shell) \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.Bilayert4K`
  - F2G2 intralayer hopping t4 (fourth shell) \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.Bilayert5KA`
  - F2G2 intralayer hopping t5 on sublattice A (fifth shell) \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.Bilayert5KB`
  - F2G2 intralayer hopping t5 on sublattice B (fifth shell) \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.BilayertK6A`
  - F2G2 intralayer hopping t6 (sixth shell) \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.BilayertK7A`
  - F2G2 intralayer hopping t7 (seventh shell) \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.BilayertK8A`
  - F2G2 intralayer hopping t8 (eighth shell) \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Stack.T3BG`
  - Here: selects the T3BG intralayer parameter set.
* - `Neigh.LayerDistFactor`
  - Here: read by the hopping set-up.
* - `Interlayer.KaxirasCutoff`
  - BLKaxiras / BLSrivani: cut-off of the interlayer hopping (compared as a squared length).
* - `Moire.OffDiag`
  - Here: the off-diagonal moire bond term is added in HamHopping.
* - `GBN.OffDiag`
  - Here: the G/hBN off-diagonal bond term is added in HamHopping.
* - `Structure.AddDisplacements`
  - Here: read for the bond terms.
* - `Moire.OffDiagMidpoint`
  - Here: the bond term is evaluated at the bond midpoint.
* - `Interlayer.AddPressureDependence`
  - BLKaxiras / BLSrivani: pressure (interlayer-distance) dependent parameters of Carr et al. (default .true.).
* - `Interlayer.SwitchV3Sign`
  - BLKaxiras / BLSrivani: change the sign of the V3 term.
* - `Interlayer.OldParameterSet`
  - BLKaxiras / BLSrivani: older parameter set of the fitting functions.
* - `Interlayer.DeactivateV6`
  - BLKaxiras / BLSrivani: drop the V6 term.
* - `Interlayer.DeactivateV3`
  - BLKaxiras / BLSrivani: drop the V3 term.
* - `Intralayer.MayouIntralayer`
  - Intralayer hoppings of the Mayou model.
* - `Interlayer.LegacyAngleSearch`
  - BLKaxiras / BLSrivani: bond angles from the index-based search with a 1.43 Angstrom cut-off used before November 2026 (non-Hermitian on relaxed cells).
* - `Stack.EncapsulatedFourLayersF2G2`
  - Stack.EncapsulatedFourLayers: Bernal F2G2 coupling between the two graphene layers (2-3).
* - `Interlayer.OnlyV0`
  - BLKaxiras / BLSrivani: keep only the V0 term (V3, V6 set to zero).
* - `TBG.Angle`
  - Here: twist angle of the analytic tBG bond terms \[degrees\].
* - `Structure.RotationAngle`
  - Rotation angle \[degrees\] of the interlayer model branches that use it.
* - `Interlayer.FindThetasGeometrically`
  - BLKaxiras / BLSrivani: bond angles theta12, theta21 from the geometry of the bond.
* - `Interlayer.RenormalizeHoppingFactorAAp`
  - Factor on AA' interlayer hoppings; overwritten per bond by the distance law in GBNtwoLayers.
* - `Interlayer.RenormalizeHoppingFactorBBp`
  - Factor on BB' interlayer hoppings; overwritten per bond by the distance law in GBNtwoLayers.
* - `Interlayer.RenormalizeHoppingFactorABp`
  - Factor on AB' interlayer hoppings; overwritten per bond by the distance law in GBNtwoLayers.
* - `Interlayer.RenormalizeHoppingFactorBAp`
  - Factor on BA' interlayer hoppings; overwritten per bond by the distance law in GBNtwoLayers.
* - `Interlayer.NewFittingFunctions`
  - BLKaxiras / BLSrivani: sublattice-resolved fitting functions a exp(-((r-b)/c)^4) cos(d r - h).
* - `Interlayer.RenormalizeCoupling`
  - Multiply the interlayer hoppings by Interlayer.CouplingFactor.
* - `Interlayer.CouplingFactor`
  - Interlayer.RenormalizeCoupling: factor on every interlayer hopping. 0.895 puts the magic angle of relaxed tBG at 1.08 deg; a rigid lattice needs about 1.008.
* - `Interlayer.CouplingFactor2`
  - Interlayer.DifferentCouplings: factor on the 2-3 interface (CouplingFactor on 1-2 and 3-4).
* - `Interlayer.DifferentCouplings`
  - Use Interlayer.CouplingFactor2 for the 2-3 interface of a four-layer stack.
* - `Interlayer.KoshinoSR`
  - Koshino two-centre interlayer hopping scaled per bond by exp((\|dz\| - 3.25)/1.34).
* - `Interlayer.Onlyvppsigma`
  - Two-centre interlayer hopping without the pi term (V\_pp\_pi = 0).
* - `Structure.BernalReadXYZ`
  - ReadXYZ Bernal bilayer: use the Bernal interlayer branch (dimer sites get the F2G2 AB set).
* - `Intralayer.DeactivateIntrasublattice`
  - Drop the hoppings between sites of the same sublattice.
* - `Intralayer.DeactivateIntersublattice`
  - Drop the hoppings between sites of different sublattices.
* - `Intralayer.DeactivateIntraSublatticeForC`
  - Drop same-sublattice hoppings between carbon atoms.
* - `Intralayer.T2Value`
  - Read but not used at present.
* - `Hopping.tA1B2`
  - Multilayer interlayer hopping t\_A1B2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA1A2`
  - Multilayer interlayer hopping t\_A1A2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tB1B2`
  - Multilayer interlayer hopping t\_B1B2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA1A3`
  - Multilayer interlayer hopping t\_A1A3 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tB1A3`
  - Multilayer interlayer hopping t\_B1A3 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA1B2_2`
  - Multilayer interlayer hopping t\_A1B2\_2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA1A2_2`
  - Multilayer interlayer hopping t\_A1A2\_2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tB1B2_2`
  - Multilayer interlayer hopping t\_B1B2\_2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA1A3_2`
  - Multilayer interlayer hopping t\_A1A3\_2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tB1A3_2`
  - Multilayer interlayer hopping t\_B1A3\_2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA1B3`
  - Multilayer interlayer hopping t\_A1B3 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tB1A2`
  - Multilayer interlayer hopping t\_B1A2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA1A1_1`
  - Multilayer interlayer hopping t\_A1A1\_1 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tB1B1_1`
  - Multilayer interlayer hopping t\_B1B1\_1 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA2A2_1`
  - Multilayer interlayer hopping t\_A2A2\_1 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA1B3_1`
  - Multilayer interlayer hopping t\_A1B3\_1 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tB1A2_1`
  - Multilayer interlayer hopping t\_B1A2\_1 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA1A1_2`
  - Multilayer interlayer hopping t\_A1A1\_2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tB1B1_2`
  - Multilayer interlayer hopping t\_B1B1\_2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA2A2_2`
  - Multilayer interlayer hopping t\_A2A2\_2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tA1B3_2`
  - Multilayer interlayer hopping t\_A1B3\_2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Hopping.tB1A2_2`
  - Multilayer interlayer hopping t\_B1A2\_2 \[eV\] (sublattice and layer index; \_1/\_2: further neighbour shells) of the stack models in HamHopping.
* - `Interlayer.SetnLayersToZero`
  - Set the number of layers to 0 (small-bulk neighbour search on a non-bulk system).
* - `Moire.PotCabG`
  - Amplitude of the off-diagonal moire term on graphene bonds in the stack models \[eV\].
* - `Interlayer.LegacyBernalTAB1`
  - Bernal interlayer branch: the coupling used before November 2026 (gamma1 = F2G2 tAB1 = 0.1391 eV and nothing else) instead of the full Bernal F2G2 interlayer set.
* - `Strain.RealStrain`
  - Intralayer hoppings t exp(-beta (d/d\_ref - 1)), d = 3-D bond length (Pereira et al., PRB 80 045401); d\_ref per shell from Strain.RealStrainReferenceLatticeConstant(BN). Off: constant hoppings.
* - `Strain.OnlyFirstNeighborRealStrain`
  - Strain.RealStrain on the first-neighbour bonds only.
* - `Strain.RealStrainBeta`
  - Strain.RealStrain: decay exponent beta (default 3.37).
* - `Strain.PeriodicStrain`
  - Sinusoidal periodic strain of the hoppings (Strain.PeriodicStrainU0, PeriodicStrainPeriod).
* - `Strain.RealisticBubbles`
  - Strain bubbles with in-plane strain and an on-site shift (Strain.Bubble\* keys).
* - `Moire.Fourier`
  - Here: bond part of the Moire.Fourier term, added before the Peierls phase.
* - `Output.WriteDataFiles`
  - Here: write the hopping and on-site tables (generate.e, \<prefix>.s.mag, ...).
```

**`HamHopHaldane`** (internal procedure of `HamHopping`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L7043))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Haldane.PaperOrientation`
  - Haldane: multiply the model's own second-neighbour hoppings by exp(+-i phi) (otherwise add T2 exp(+-i phi)).
* - `Haldane.T2Complex`
  - Haldane: imaginary second-neighbour amplitude (see HamHopHaldane).
* - `Haldane.T2`
  - Haldane second-neighbour amplitude T2 \[eV\] (units of g0 before 2026-10-14); not used with Haldane.PaperOrientation.
* - `Haldane.OppositePhase`
  - Haldane: opposite sign of the phase phi.
* - `Haldane.BothLayers`
  - Deprecated (use Haldane.LayerControl / HaldaneLayers): the Haldane term on every layer.
* - `Haldane.Flux`
  - Haldane.SpecifyFlux: flux giving the phase 2 pi Flux / Phi\_0 (2 pi Flux / (pi/3) with Haldane.SetFluxQ).
* - `Haldane.SetFluxQ`
  - Haldane.SpecifyFlux: the phase is 2 pi Haldane.Flux / (pi/3).
* - `Haldane.Phase`
  - Haldane.SpecifyPhase: Haldane phase phi \[rad\].
```

**`HamHopRealisticBubbles`** (internal procedure of `HamHopping`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L7229))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Strain.BigBubble`
  - Strain.RealisticBubbles: one big bubble.
* - `Strain.ManyBubbles`
  - Strain.RealisticBubbles: many bubbles at random positions (Strain.BubblePercentage).
* - `Strain.BubbleInPlaneStrain`
  - Strain.RealisticBubbles: include the in-plane strain of each bubble.
* - `Strain.BubblePercentage`
  - Strain.ManyBubbles: percentage of the atoms that are bubble centres.
* - `Strain.BubbleGaussian`
  - Strain.RealisticBubbles: add a Gaussian potential at each bubble (Strain.BubbleGaussPot\*).
* - `Strain.BubbleSigmaR`
  - Strain.RealisticBubbles: width of the bubble profile \[Angstrom\].
* - `Strain.BubbleC`
  - Strain.RealisticBubbles: amplitude C of the in-plane strain.
* - `Strain.BubbleShift`
  - Strain.RealisticBubbles: on-site shift inside a bubble \[eV\] (units of g0 before 2026-11-07).
* - `Strain.BubbleRadius`
  - Strain.RealisticBubbles: radius of a bubble \[Angstrom\].
```

**`HamHopPeriodicStrain`** (internal procedure of `HamHopping`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L7444))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Strain.PeriodicStrainU0`
  - Strain.PeriodicStrain: amplitude u0 of the displacement \[Angstrom\].
```

**`HamHopRealStrain`** (internal procedure of `HamHopping`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L7500))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Strain.RealStrainReferenceLatticeConstant`
  - Strain.RealStrain: graphene lattice constant of the reference bond lengths \[Angstrom\]; the lattice of the hopping fit, 2.4389777651302801 for the F2G2 sets.
* - `Strain.RealStrainReferenceLatticeConstantBN`
  - Strain.RealStrain: hBN lattice constant of the reference bond lengths \[Angstrom\].
* - `Strain.RealStrainLegacyFourthShell`
  - Strain.RealStrain: fixed 3.75696 Angstrom reference for the fourth shell, as before November 2026.
```

**`HamHopShellsFromRigid`** (internal procedure of `HamHopping`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L7693))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Strain.ShellsRigidFile`
  - Strain.ShellsFromRigidPositions: rigid reference structure (xyz, same atom order).
```

**`HamHopInterlayerParameters`** (internal procedure of `HamHopping`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L7911))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Interlayer.UseBNGSrivani`
  - BLSrivani: the G/hBN parameter set.
* - `Interlayer.Oppositedxdy`
  - BLKaxiras / BLSrivani: opposite sign of the bond vector in the angle calculation.
* - `Interlayer.SublatticeDependent`
  - BLKaxiras / BLSrivani: sublattice-dependent parameters.
* - `Interlayer.UseThetaIJ`
  - BLKaxiras / BLSrivani: angle variant theta\_ij (see HamHopInterlayerParameters).
* - `Interlayer.UseTheta`
  - BLKaxiras / BLSrivani: angle variant theta (see HamHopInterlayerParameters).
* - `Interlayer.SublatticeIndependent`
  - BLKaxiras / BLSrivani: sublattice-independent parameters.
* - `Interlayer.UseOnlyVAB`
  - BLKaxiras / BLSrivani: use the V\_AB parameters for every sublattice pair.
```

**`HamHopF2G2sParameters`** (internal procedure of `HamHopping`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L8253))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Intralayer.SingleLayert2KSL`
  - F2G2 intralayer hopping t2, second shell \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.SingleLayert3K`
  - F2G2 intralayer hopping t3, third shell \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.SingleLayert4K`
  - F2G2 intralayer hopping t4, fourth shell \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.SingleLayert5KSL`
  - F2G2 intralayer hopping t5, fifth shell \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.hBNt2B`
  - F2G2 intralayer hopping t2 of hBN on boron \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.hBNt2N`
  - F2G2 intralayer hopping t2 of hBN on nitrogen \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.hBNt3`
  - F2G2 intralayer hopping t3 of hBN \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.hBNt4`
  - F2G2 intralayer hopping t4 of hBN \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.hBNt5B`
  - F2G2 intralayer hopping t5 of hBN on boron \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.hBNt5N`
  - F2G2 intralayer hopping t5 of hBN on nitrogen \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.hBNt6`
  - F2G2 intralayer hopping t6 of hBN \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.hBNt7`
  - F2G2 intralayer hopping t7 of hBN \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
* - `Intralayer.hBNt8`
  - F2G2 intralayer hopping t8 of hBN \[eV\]; the default depends on the parameter set read here (graphene, hBN, G/hBN fits).
```

**`HamHopBondTermsTBG`** (internal procedure of `HamHopping`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L12444))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Moire.PotCabG`
  - Here (tBG bond terms): amplitude of the off-diagonal term on graphene bonds \[eV\].
```

## See also

{ref}`HamInit <api-ham-haminit>`

## Calls

{ref}`AtomsSetCart <api-atoms-atomssetcart>`, {ref}`AtomsSetFrac <api-atoms-atomssetfrac>`, [ComputePzNormals](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4295) ({ref}`ham <api-ham>`), [cross](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L13566) ({ref}`ham <api-ham>`), [distanceDependentC](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4062) ({ref}`ham <api-ham>`), [fitFunSrivani](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4041) ({ref}`ham <api-ham>`), {ref}`GaussPotDefinedPositions <api-gauss-gausspotdefinedpositions>`, [GetUnit](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L606) ({ref}`file <api-file>`), {ref}`HaldaneEnabledForLayer <api-ham-haldaneenabledforlayer>`, [HamCheckFinite](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4432) ({ref}`ham <api-ham>`), [HamCheckHermiticity](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4491) ({ref}`ham <api-ham>`), [HamHaldaneChirality](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L13063) ({ref}`ham <api-ham>`), [HamHopMoireFourier](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5009) ({ref}`ham <api-ham>`), [HamHopValleyZeeman](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5061) ({ref}`ham <api-ham>`), {ref}`HamInputEnergyEV <api-ham-haminputenergyev>`, [HamInputHBN](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L13105) ({ref}`ham <api-ham>`), [HamMinImageXY](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L176) ({ref}`ham <api-ham>`), {ref}`HamPeierlsPeriodicGauge <api-ham-hampeierlsperiodicgauge>`, [HamPrintNormalStats](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4357) ({ref}`ham <api-ham>`), [HamRigidBondXY](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L153) ({ref}`ham <api-ham>`), [HamSOCBondTables](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L12929) ({ref}`ham <api-ham>`), [HamStrainReferenceDistance](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L201) ({ref}`ham <api-ham>`), [HamUpdateNormals](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4578) ({ref}`ham <api-ham>`), {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), [interlayerBLAA](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5328) ({ref}`ham <api-ham>`), [interlayerBLAB](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5261) ({ref}`ham <api-ham>`), [interlayerBLBA](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5294) ({ref}`ham <api-ham>`), [interlayerTwoCenter](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4094) ({ref}`ham <api-ham>`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, [MoireBondDeltaMidpoint](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4822) ({ref}`ham <api-ham>`), {ref}`num2str <api-mio-interface-num2str>`, [offdiago](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4752) ({ref}`ham <api-ham>`), [offdiago2](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4783) ({ref}`ham <api-ham>`), [offdiagoGBN](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5103) ({ref}`ham <api-ham>`), [offdiagotBG](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5133) ({ref}`ham <api-ham>`), {ref}`MIO_StringComp <api-string-interface-stringcomp>` (= `StringComp`), {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_PrintDebug <api-sys-sysprintdebug>` (= `SysPrintDebug`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`), {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

{ref}`CalcSelect <api-calc-calcselect>`
