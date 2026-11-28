---
edit_url: null
---

(api-ham-hamonsite)=
# `HamOnSite`

```{index} single: HamOnSite (subroutine in ham)
```

```fortran
subroutine HamOnSite()
```

- **Module:** {ref}`ham <api-ham>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/ham.F90, line 756`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L756)

Compute and add on-site energy terms to H0.

Controls the on-site contributions to the Hamiltonian. Values are added to H0 (never reinitialized) to preserve composition of effects (e.g., strain + gates). Handles multiple layer configurations and optional random/moire strain.

```{note}
Always update H0 additively; do not reinitialize.
```

**Side effects**

- Updates H0(:) in place
- May switch between fractional and cartesian coordinates temporarily
- Applies optional Gaussian/interface potentials

## Input keys read (201)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Stack.GBNtwoLayersF2G2s`
  - Stack.GBNtwoLayers with F2G2 intralayer parameters for both graphene and hBN.
* - `Stack.GBNtwoLayers`
  - Graphene on one hBN layer (ReadXYZ cells): intralayer shells from the in-plane bond length, distance-dependent C-B and C-N interlayer hoppings.
* - `Stack.BNBNtwoLayers`
  - Two hBN layers: refused (no interlayer hopping is implemented).
* - `Stack.T2GBN`
  - Stack model flag: selects the twisted bilayer graphene on hBN (three layers) branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.T2BG`
  - Stack model flag; read here but not used in HamOnSite.
* - `Stack.T3BG`
  - Stack model flag: selects the six-layer T3BG branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.EncapsulatedThreeLayers`
  - Stack model flag: selects the hBN-encapsulated graphene (three layers) branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.EncapsulatedFourLayers`
  - Stack model flag: selects the hBN-encapsulated bilayer (four layers) branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.EncapsulatedFiveLayers`
  - Stack model flag: selects the five-layer encapsulated stack branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.EncapsulatedSixLayers`
  - Stack model flag: selects the six-layer encapsulated stack branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.EncapsulatedSevenLayers`
  - Stack model flag: selects the seven-layer encapsulated stack branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Potential.UseLayerSpecificOnsiteEnergyTerms`
  - Use the per-layer carbon on-site energies TB.OnSiteC1\_LB, C2\_LB, C1\_LT, C2\_LT.
* - `Stack.T3GwithBN`
  - Stack model flag: selects the twisted trilayer graphene with hBN branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.BNt2GBN`
  - Stack model flag: selects the hBN / twisted bilayer graphene / hBN branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Moire.DontUseInplaneMoire`
  - Skip the harmonic in-plane moire on-site term for GBNtwoLayers and the encapsulated stacks.
* - `Stack.TwoLayers`
  - Stack model flag: selects the two-layer graphene branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.ThreeLayers`
  - Stack model flag: selects the three-layer graphene branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.FourLayersSandwiched`
  - Stack model flag: selects the four-layer sandwiched graphene branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.HelicalTwistedMBM`
  - Stack model flag: selects the helical-twisted four-layer (mono-bi-mono) branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.FiveLayersSandwiched`
  - Stack model flag: selects the five-layer graphene branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.SixLayersSandwiched`
  - Stack model flag: selects the six-layer graphene branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.SevenLayersSandwiched`
  - Stack model flag: selects the seven-layer graphene branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.EightLayersSandwiched`
  - Stack model flag: selects the eight-layer graphene branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.TenLayersSandwiched`
  - Stack model flag: selects the ten-layer graphene branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Stack.TwentyLayersSandwiched`
  - Stack model flag: selects the twenty-layer graphene branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Intralayer.ForceBilayerF2G2Intralayer`
  - Multilayer stacks: Bernal-bilayer F2G2 intralayer parameters instead of the single-layer ones.
* - `Stack.MiddleTwist`
  - Single-layer F2G2 parameters in every layer of the multilayer stacks; default .true. unless Intralayer.ForceBilayerF2G2Intralayer. A threeLayers stack then gets the averaged on-site energy.
* - `Stack.FourLayers`
  - Stack model flag: selects the four-layer graphene branches of the on-site and hopping terms (HamOnSite, HamHopping).
* - `Structure.ReadRigidXYZ`
  - Here: the rigid positions enter the on-site terms of the displacement-based models.
* - `Moire.Potential`
  - Moire on-site potential (HamOnSiteMoirePotential): the G/hBN effective model needs one of MoireH0AndHZ, MoireOnlyH0, MoireOnlyHZ, MoireNoH0AndHZ. No effect for ReadXYZ + GBNtwoLayers.
* - `Moire.PreFactor`
  - Moire.Potential: factor on the amplitudes of the moire potential.
* - `TBG.Diag`
  - Stacking-dependent on-site term of twisted bilayer graphene (from the local stacking displacements).
* - `TBG.DiagPRB`
  - TBG.Diag: use the amplitudes Interlayer.CAA, CBB, ... instead of the default set.
* - `BNBN.Diag`
  - Stacking-dependent on-site term of an hBN bilayer.
* - `Moire.LatticePercent`
  - Lattice mismatch epsilon of the moire potential, d(r) = ((1+eps) R(theta) - 1) r; commensurate G/hBN in a CellSize cell wants -1/CellSize.
* - `Moire.BilayerElectricField`
  - Layer potential (displacement field) of Moire.BilayerElectricShift; the profile follows the stack flag (e.g. Stack.FiveLayersSandwiched). Without this flag the shift does nothing.
* - `Moire.BilayerElectricShift`
  - Moire.BilayerElectricField: potential difference between the outer layers \[eV\].
* - `Moire.BilayerElectricFieldInvert`
  - Moire.BilayerElectricField: invert the field for atoms beyond half the cell in x.
* - `Structure.CellSize`
  - Here: half the cell length in x, CellSize\*SuperCell\*aG/2, for Moire.BilayerElectricFieldInvert.
* - `Trilayer.FanZhang`
  - On-site terms of the Fan-Zhang trilayer model (HamOnSiteTrilayerFanZhang).
* - `Potential.AddSublatticeMassterm`
  - Sublattice mass: Species 1 gets -m/2, Species 2 +m/2, m = Potential.SublatticeMassterm.
* - `Potential.OnlyBottomLayerMassTerm`
  - Potential.AddSublatticeMassterm on the bottom layer only.
* - `Potential.SublatticeMassterm`
  - Potential.AddSublatticeMassterm: mass m \[eV\].
* - `Potential.AddOnsiteEnergyShift`
  - Add Potential.OnsiteEnergyShift to every on-site energy.
* - `Potential.OnsiteEnergyShift`
  - Potential.AddOnsiteEnergyShift: shift \[eV\].
* - `Strain.Bubbles`
  - Disorder: shift the on-site energy by Disorder.OnsiteShift on randomly chosen cells of a CheckerDivider x CheckerDivider checkerboard (fraction Disorder.CheckerDensity).
* - `Disorder.OnsiteShift`
  - Strain.Bubbles: on-site shift of the chosen cells \[eV\] (units of g0 before 2026-10-14).
* - `Disorder.CheckerDivider`
  - Strain.Bubbles: the cell is divided into D x D checkerboard cells.
* - `Disorder.CheckerDensity`
  - Strain.Bubbles: fraction of checkerboard cells that are shifted.
* - `Potential.PNP`
  - Step p-n-p potential along x: -PNPAmp outside, +PNPAmp inside the central region (HamOnSitePNP).
* - `Potential.PNPKink`
  - Smooth p-n-p junction along x: Potential.PNPAmp x tanh steps of width Potential.PNPDelta at L/4 and 3L/4 (HamOnSitePNPKink).
* - `Potential.SinusModulation`
  - Legacy: add sin(2 pi N x / L) with amplitude 1 g0 (no amplitude key); do not use.
* - `Potential.SinusNumberOfPeriod`
  - Potential.SinusModulation: number N of periods in the cell.
* - `Potential.SinusModulationUsingPeriod`
  - Add Potential.SinusFactor x sin(2 pi x / Potential.SinusModulationPeriod); x is Cartesian, so it is lattice-periodic only for suitable periods.
* - `Potential.SinusModulationPeriod`
  - Period of the sine modulation \[Angstrom\] (absolute length, not a fraction of the cell).
* - `Potential.SinusModulationAddMassTerm`
  - Replace the sine modulation by its staggered form SinusMassTermFactor x sin(...) x (-1)^Species.
* - `Potential.SinusFactor`
  - Amplitude of the sine modulation \[eV\] (units of g0 before 2026-10-14).
* - `Potential.SinusMassTermFactor`
  - Amplitude of the staggered sine modulation \[eV\] (units of g0 before 2026-10-14).
* - `Potential.SinusModulationUsingPeriodYDirection`
  - As Potential.SinusModulationUsingPeriod, along y.
* - `Potential.CosinusModulationUsingPeriod`
  - Add Potential.CosinusFactor x cos(2 pi x / Potential.CosinusModulationPeriod).
* - `Potential.CosinusModulationPeriod`
  - Period of the cosine modulation \[Angstrom\].
* - `Potential.CosinusModulationAddMassTerm`
  - Replace the cosine modulation by CosinusMassTermFactor x cos(...) x (-1)^Species.
* - `Potential.CosinusFactor`
  - Amplitude of the cosine modulation \[eV\] (units of g0 before 2026-10-14).
* - `Potential.CosinusMassTermFactor`
  - Amplitude of the staggered cosine modulation \[eV\] (units of g0 before 2026-10-14).
* - `Potential.CosinusModulationUsingPeriodYDirection`
  - As Potential.CosinusModulationUsingPeriod, along y.
* - `Potential.SquareFunction`
  - Staggered potential -+(-1)^Species A following the sign of sin(2 pi N x/L) sin(2 pi N y/L) (HamOnSiteSquareFunction).
* - `Potential.SquareFunction2`
  - Staggered square-wave potential on honeycomb-sized regions (HamOnSiteSquareFunction2).
* - `Potential.Zterm1D`
  - Staggered potential that changes sign with sin(2 pi N x / L) (HamOnSiteZterm1D).
* - `Potential.Zterm1DKink`
  - Staggered potential with smooth tanh kinks along x (HamOnSiteZterm1DKink); written for two periods.
* - `Potential.SquareChecker2219`
  - Staggered checkerboard potential with fixed periods 22 x 19 in the cell (HamOnSiteSquareChecker).
* - `Disorder.SublatticeDisorder`
  - Random on-site disorder on a fraction Disorder.SublattPct of the Species-2 atoms.
* - `Disorder.SublattAmp`
  - Disorder.SublatticeDisorder: width w of the random energies \[eV\] (units of g0 before 2026-11-07).
* - `Disorder.SublattPct`
  - Disorder.SublatticeDisorder: fraction of the atoms affected.
* - `Disorder.Anderson`
  - Anderson disorder: random on-site energies of width Disorder.AndersonAmp on every atom.
* - `Disorder.AndersonAmp`
  - Disorder.Anderson: width w of the uniform distribution \[eV\] (units of g0 before 2026-11-07).
* - `Disorder.DeltaDisorder`
  - Delta-like impurity disorder (HamOnSiteDeltaDisorder; Disorder.DeltaAmp, DeltaSkewFactor).
* - `Disorder.GaussDisorder`
  - Gaussian impurities of strength Disorder.GaussPotStrength (Disorder.GaussPot\* keys).
* - `Disorder.GaussPotStrength`
  - Disorder.GaussDisorder: amplitude w of the Gaussian impurities \[eV\] (units of g0 before 2026-11-07).
* - `Potential.FourLayerOnsiteShifts`
  - Add the on-site shifts Stack.FourLayerShift1..4 to the atoms of layers 1..4 (four-layer stacks only), e.g. for a screened displacement-field profile.
* - `Stack.FourLayerShift1`
  - Potential.FourLayerOnsiteShifts: on-site shift of layer 1 \[eV\].
* - `Stack.FourLayerShift2`
  - Potential.FourLayerOnsiteShifts: on-site shift of layer 2 \[eV\].
* - `Stack.FourLayerShift3`
  - Potential.FourLayerOnsiteShifts: on-site shift of layer 3 \[eV\].
* - `Stack.FourLayerShift4`
  - Potential.FourLayerOnsiteShifts: on-site shift of layer 4 \[eV\].
* - `Potential.HelicalTwistedMBM_CDW`
  - C3-symmetric trigonal CDW, (A/3) \[cos 2 pi S f1 + cos 2 pi S f2 + cos 2 pi S (f1+f2)\] in fractional coordinates, S = SuperCell; Potential.CDWAmplitude, CDWUseMassTerm.
* - `Moire.CDW`
  - On-site CDW \[A + A\_m (-1)^Species\] (1/nQ) sum\_q cos(2 pi (h f1 + k f2)/D + phase) on every atom, q from the block moireCDW.Qvectors in units of the cell's reciprocal lattice (moireCDW.\* keys).
* - `Moire.Fourier`
  - General moire term given by Fourier components: on-site part here, bond part in HamHopping.
```

**`HamOnSiteMoireCDW`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L1358))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `moireCDW.Amplitude`
  - Moire.CDW: scalar amplitude A \[eV\].
* - `moireCDW.MassAmplitude`
  - Moire.CDW: staggered amplitude A\_m, multiplied by (-1)^Species \[eV\].
* - `moireCDW.Phase`
  - Moire.CDW: phase of every cosine \[rad\].
* - `moireCDW.Denominator`
  - Moire.CDW: denominator D of the wave vectors q = (h, k)/D; the cell must be periodic for them.
* - `Structure.SuperCellX`
  - Here: checks that the moireCDW wave vectors are periodic in the supercell.
```

**`HamOnSiteTrigonalCDW`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L1446))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Potential.CDWAmplitude`
  - Potential.HelicalTwistedMBM\_CDW: amplitude A \[eV\].
* - `Potential.CDWPeriod`
  - Read but not used: the trigonal CDW period is set by Structure.SuperCell.
* - `Potential.CDWUseMassTerm`
  - Potential.HelicalTwistedMBM\_CDW: multiply the CDW by (-1)^Species.
```

**`HamOnSiteDeltaDisorder`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L1499))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Disorder.DeltaAmp`
  - Disorder.DeltaDisorder: amplitude w \[eV\] (units of g0 before 2026-10-14).
* - `Disorder.DeltaSkewFactor`
  - Disorder.DeltaDisorder: skew factor of the distribution.
```

**`HamOnSiteSquareChecker`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L1546))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Potential.AmplitudeOfSquare`
  - Amplitude A of the staggered square / checker potentials \[eV\] (units of g0 before 2026-10-14).
```

**`HamOnSiteZterm1DKink`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L1590))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Potential.Zterm1DAmp`
  - Potential.Zterm1DKink: amplitude of the tanh profile \[eV\] (units of g0 before 2026-10-14).
* - `Potential.Zterm1DDelta`
  - Potential.Zterm1DKink: width of the tanh steps \[Angstrom\].
```

**`HamOnSiteZterm1D`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L1631))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Potential.AmplitudeOfSquare3`
  - Potential.Zterm1D: amplitude \[eV\] (units of g0 before 2026-10-14).
```

**`HamOnSiteSquareFunction2`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L1664))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Potential.NumberOfWidthHoneycomb`
  - Potential.SquareFunction2: size of the regions in graphene cells.
* - `Potential.AmplitudeOfSquare2`
  - Potential.SquareFunction2: amplitude \[eV\] (units of g0 before 2026-10-14).
* - `Potential.PhaseOfSquareX`
  - Potential.SquareFunction2: phase of the x factor \[rad\].
* - `Potential.PhaseOfSquareY`
  - Potential.SquareFunction2: phase of the y factor \[rad\].
* - `Potential.TwoDimension`
  - Potential.SquareFunction2: modulation along x and y (otherwise x only).
* - `Potential.AddZTerm`
  - Potential.SquareFunction(2): add the staggered z term.
* - `Potential.ArmChairShape`
  - Potential.SquareFunction2: regions of armchair shape.
```

**`HamOnSiteSquareFunction`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L1760))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Potential.TwoDimensional`
  - Potential.SquareFunction: modulation along x and y (otherwise x only).
```

**`HamOnSitePNPKink`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L1832))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Potential.PNPAmp`
  - Amplitude of the p-n-p potential \[eV\] (units of g0 before 2026-10-14).
```

**`HamOnSiteTrilayerFanZhang`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L1939))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Trilayer.U1`
  - Trilayer.FanZhang: on-site potential of layer 1 \[eV\].
* - `Trilayer.U2`
  - Trilayer.FanZhang: on-site potential of layer 2 \[eV\].
* - `Trilayer.U3`
  - Trilayer.FanZhang: on-site potential of layer 3 \[eV\].
* - `Trilayer.Delta`
  - Trilayer.FanZhang: Delta term of the model \[eV\].
```

**`HamOnSiteMoirePotential`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L2082))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Moire.Sachs`
  - Moire potential of Sachs et al., PRB 84, 195414 (2011), amplitudes Moire.PotA, PotB, PotC.
* - `Moire.PotA`
  - Moire.Sachs: amplitude A \[eV\].
* - `Moire.PotB`
  - Moire.Sachs: amplitude B \[eV\].
* - `Moire.PotC`
  - Moire.Sachs: constant term C \[eV\].
* - `Moire.Jeil`
  - G/hBN moire potential of Jung et al., PRB 89, 205414 (2014): amplitudes Moire.PotC0, PotCz, PotCab and phases Moire.PotPhi0, PotPhiz, PotPhiab (the effective model).
* - `Moire.PreFactor`
  - Here: factor on Moire.PotC0, PotCz, PotCab.
* - `Moire.PotC0`
  - Moire.Jeil: amplitude C0 of the scalar moire term \[eV\].
* - `Moire.PotCz`
  - Moire.Jeil: amplitude Cz of the mass (sublattice-odd) term \[eV\].
* - `Moire.PotCab`
  - Moire.Jeil: amplitude Cab of the off-diagonal (bond) term \[eV\].
* - `Moire.PotPhi0`
  - Moire.Jeil: phase of the C0 term \[rad\].
* - `Moire.PotPhiz`
  - Moire.Jeil: phase of the Cz term \[rad\].
* - `Moire.PotPhiab`
  - Moire.Jeil: phase of the Cab term \[rad\].
* - `Moire.DistanceDependentEffectiveModel`
  - Moire.Jeil amplitudes scaled with the local interlayer distance z: C(z) = C exp(-B (z - z0)).
* - `Moire.BfactorC0`
  - Moire.DistanceDependentEffectiveModel: decay B of C0 \[1/Angstrom\].
* - `Moire.BfactorCz`
  - Moire.DistanceDependentEffectiveModel: decay B of Cz \[1/Angstrom\].
* - `Structure.InterlayerDistance`
  - Here: reference distance z0 of Moire.DistanceDependentEffectiveModel \[Angstrom\].
* - `Potential.SublatticeBasis`
  - Moire potential from the AA/BB sublattice form (diagoH0FromAABB) instead of the C0/phase form.
* - `Structure.AddDisplacements`
  - Add the stacking displacements read from the file to the analytic ones (needs a displacement file).
* - `Moire.Trilayer`
  - Two moire potentials of a trilayer (needs Moire.Sachs; Moire.PotA1..B2, Moire.Trin1/m1/n2/m2).
* - `Moire.RemoveTopMoireInL2`
  - Leave out the moire term from the top interface in layer 2 (see HamOnSiteMoirePotential).
* - `BNBN.UseDisplacementFile`
  - BNBN.Diag with per-atom stacking displacements from Structure.DisplacementsFile.
* - `GBN.UseHarmonicApprox`
  - Stack.GBNtwoLayers: harmonic (first-shell) G/hBN on-site moire term; TBG.Diag and BNBN.Diag are then switched off.
* - `GBN.Angle`
  - Twist angle of the analytic G/hBN stacking displacements \[degrees\].
* - `Stack.GlobalTwist`
  - Rotate the stacking displacements read from the file by this angle before use \[degrees\].
* - `Stack.GlobalTwist2`
  - As Stack.GlobalTwist, for the second set of displacements \[degrees\].
* - `Stack.GlobalPhiL1`
  - Rotation of the displacements of layer 1 \[degrees\] (layer-resolved stacks).
* - `Stack.GlobalPhiL2`
  - Rotation of the displacements of layer 2 \[degrees\].
* - `Stack.GlobalPhiL2a`
  - Rotation of the bottom-interface displacements of layer 2 \[degrees\] (encapsulated layers).
* - `Stack.GlobalPhiL2b`
  - Rotation of the top-interface displacements of layer 2 \[degrees\].
* - `Stack.GlobalPhiL3`
  - Rotation of the displacements of layer 3 \[degrees\].
* - `Structure.RotationAngle`
  - Rotation angle \[degrees\] used by the hopping terms (HamHopping); read here but not used.
* - `TBG.SwitchDxDy`
  - TBG.Diag: exchange the x and y components of the stacking displacement.
* - `Interlayer.CAA`
  - Stacking-dependent on-site term (TBG.Diag, BNBN.Diag, or the G/hBN default set): amplitude C\_AA on sublattice A of layer 1 \[eV\]; the defaults depend on the model.
* - `Interlayer.CBB`
  - Stacking-dependent on-site term (TBG.Diag, BNBN.Diag, or the G/hBN default set): amplitude C\_BB on sublattice B of layer 1 \[eV\]; the defaults depend on the model.
* - `Interlayer.CAA0`
  - Stacking-dependent on-site term (TBG.Diag, BNBN.Diag, or the G/hBN default set): constant part C\_AA0 on sublattice A of layer 1 \[eV\]; the defaults depend on the model.
* - `Interlayer.CBB0`
  - Stacking-dependent on-site term (TBG.Diag, BNBN.Diag, or the G/hBN default set): constant part C\_BB0 on sublattice B of layer 1 \[eV\]; the defaults depend on the model.
* - `Interlayer.CApAp`
  - Stacking-dependent on-site term (TBG.Diag, BNBN.Diag, or the G/hBN default set): amplitude C\_A'A' on sublattice A of layer 2 \[eV\]; the defaults depend on the model.
* - `Interlayer.CBpBp`
  - Stacking-dependent on-site term (TBG.Diag, BNBN.Diag, or the G/hBN default set): amplitude C\_B'B' on sublattice B of layer 2 \[eV\]; the defaults depend on the model.
* - `Interlayer.CApAp0`
  - Stacking-dependent on-site term (TBG.Diag, BNBN.Diag, or the G/hBN default set): constant part on sublattice A of layer 2 \[eV\]; the defaults depend on the model.
* - `Interlayer.CBpBp0`
  - Stacking-dependent on-site term (TBG.Diag, BNBN.Diag, or the G/hBN default set): constant part on sublattice B of layer 2 \[eV\]; the defaults depend on the model.
* - `Interlayer.PhiAA`
  - Stacking-dependent on-site term: phase of the AA harmonic \[degrees\]; the defaults depend on the model.
* - `Interlayer.PhiBB`
  - Stacking-dependent on-site term: phase of the BB harmonic \[degrees\]; the defaults depend on the model.
* - `Interlayer.PhiApAp`
  - Stacking-dependent on-site term: phase of the A'A' harmonic \[degrees\]; the defaults depend on the model.
* - `Interlayer.PhiBpBp`
  - Stacking-dependent on-site term: phase of the B'B' harmonic \[degrees\]; the defaults depend on the model.
* - `TBG.Angle`
  - Twist angle of the analytic tBG stacking displacements \[degrees\] (on-site and bond terms).
* - `TBG.LegacyAnalyticOnSite`
  - Analytic tBG on-site term as before November 2026 (GBN.Angle, d in Angstrom); warning printed.
* - `Output.WriteDataFiles`
  - Here: also writes the moire on-site terms to Hjj1, Hjj2.
* - `Moire.PotA1`
  - Moire.Trilayer: amplitude A of the first moire \[eV\] (default Moire.PotA).
* - `Moire.PotB1`
  - Moire.Trilayer: amplitude B of the first moire \[eV\] (default Moire.PotB).
* - `Moire.PotA2`
  - Moire.Trilayer: amplitude A of the second moire \[eV\] (default Moire.PotA).
* - `Moire.PotB2`
  - Moire.Trilayer: amplitude B of the second moire \[eV\] (default Moire.PotB).
* - `Moire.Trin1`
  - Moire.Trilayer: n1 of the moire periods n1 aG/m1 and n1 n2 aG/m2.
* - `Moire.Trim1`
  - Moire.Trilayer: m1 of the first period n1 aG/m1 (non-zero).
* - `Moire.Trin2`
  - Moire.Trilayer: n2 of the second period n1 n2 aG/m2.
* - `Moire.Trim2`
  - Moire.Trilayer: m2 of the second period n1 n2 aG/m2 (non-zero).
* - `Structure.SingleLayerXYZ`
  - ReadXYZ: the structure is a single layer (selects the single-layer branch of the moire potential).
* - `Trilayer.FanZhang`
  - Here: selects the trilayer branch of the moire potential.
* - `Moire.OnlyH0`
  - Moire potential: apply only the scalar term H0 (one of OnlyH0, OnlyHZ, NoH0AndHZ, H0AndHZ is required).
* - `Moire.OnlyHZ`
  - Moire potential: apply only the mass term HZ.
* - `Moire.NoH0AndHZ`
  - Moire potential: apply neither H0 nor HZ (bond terms only).
* - `Moire.H0AndHZ`
  - Moire potential: apply both the scalar H0 and the mass HZ terms.
* - `Moire.BilayerTopAngle`
  - Encapsulated bilayer: twist angle of the top hBN relative to the graphene \[degrees\].
* - `Moire.BilayerBottomAngle`
  - Encapsulated bilayer: twist angle of the bottom hBN relative to the graphene \[degrees\].
* - `Moire.BLDeactivateUpperLayer`
  - Encapsulated bilayer: no moire term on layer 2.
* - `Moire.TDBLDeactivateUpperLayers`
  - Multilayer: no moire term on any layer above layer 1.
* - `Moire.Kekule`
  - Kekule-type moire on-site term (HamOnSiteMoireKekule; Moire.KekuleAngle, KekuleEpsFactor).
* - `Moire.AddSecondMoire`
  - Second moire on the same graphene (effective model): same amplitudes and phases, slot 2 = Moire.Twisted2, TwistAngle2, LatticePercent x LatticePercentFactor, SecondMoireMassFactor.
* - `Moire.LatticePercentFactor`
  - Moire.AddSecondMoire: lattice mismatch of the second moire = Moire.LatticePercent x this factor.
* - `Moire.SwitchHzjj`
  - Read but not used at present.
* - `Moire.TwistAngle`
  - Moire.Twisted: twist angle theta of the first moire, d(r) = ((1+eps) R(theta) - 1) r \[degrees\].
* - `Moire.TwistAngle2`
  - Moire.Twisted2: twist angle of the second moire (Moire.AddSecondMoire) \[degrees\].
* - `Moire.Twisted`
  - The first moire includes the twist Moire.TwistAngle (otherwise lattice mismatch only).
* - `Moire.Twisted2`
  - The second moire includes the twist Moire.TwistAngle2.
* - `Intralayer.DeactivateASubLattice`
  - Moire potential: skip the moire term on the Species-2 atoms.
* - `Intralayer.DeactivateBSubLattice`
  - Moire potential: skip the moire term on the Species-1 atoms.
* - `Moire.LayerShift1`
  - Moire potential: added to the y component of the stacking displacement of the first (bottom) moire \[Angstrom\].
* - `Moire.LayerShift2`
  - Moire potential: shift of the second moire pattern \[Angstrom\].
* - `Moire.SecondMoireRotateFirst`
  - Moire.AddSecondMoire: build the displacement with the rotation applied first (default).
* - `Moire.SecondMoireMassFactor`
  - Moire.AddSecondMoire: factor (sign) of the mass term of the second moire; \< 0 flips BN to NB.
* - `Moire.FirstMoireMassFactor`
  - Factor (sign) of the mass term of the first moire; \< 0 flips BN to NB.
* - `Potential.MinZ`
  - Read but not used at present.
* - `Moire.SymmetricPot`
  - Symmetric moire potential with amplitudes Moire.SymmA, Moire.SymmB.
* - `Moire.SymmA`
  - Moire.SymmetricPot: amplitude A \[eV\].
* - `Moire.SymmB`
  - Moire.SymmetricPot: amplitude B \[eV\].
```

**`HamOnSiteMoireKekule`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L3648))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Moire.KekuleAngle`
  - Moire.Kekule: angle of the Kekule pattern \[degrees\].
* - `Moire.KekuleEpsFactor`
  - Moire.Kekule: the lattice mismatch enters as Moire.LatticePercent / this factor.
```

**`HamOnSiteSpeciesEnergies`** (internal procedure of `HamOnSite`, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L3703))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Interlayer.TypeOfBL`
  - Interlayer hopping model: None, Koshino, Jeil, Mayou, BLKaxiras, BLSrivani, ... (see HamHopping); also selects matching on-site energies (e.g. BLKaxiras with Interlayer.UseBNGKaxiras).
* - `Interlayer.UseBNGKaxiras`
  - With TypeOfBL BLKaxiras: the G/hBN interlayer model (Eq. 4) with the hBN 8-shell intralayer set and its B / N on-site energies; needs a stack flag (e.g. Stack.TwoLayers).
```

## See also

{ref}`HamInit <api-ham-haminit>`

## Calls

{ref}`AtomsSetCart <api-atoms-atomssetcart>`, {ref}`AtomsSetFrac <api-atoms-atomssetfrac>`, [diago](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4696) ({ref}`ham <api-ham>`), [diagoGBN](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4724) ({ref}`ham <api-ham>`), [diagoH0FromAABB](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4621) ({ref}`ham <api-ham>`), [diagoHzFromAABB](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4660) ({ref}`ham <api-ham>`), [diagotBG](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5164) ({ref}`ham <api-ham>`), [distanceDependentC](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4062) ({ref}`ham <api-ham>`), [fitFunSrivani](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4041) ({ref}`ham <api-ham>`), {ref}`GaussPot <api-gauss-gausspot>`, {ref}`HamInputEnergyEV <api-ham-haminputenergyev>`, [HamMeasureLayerTwist](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4124) ({ref}`ham <api-ham>`), [HamOnSiteMoireFourier](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4978) ({ref}`ham <api-ham>`), [harmonicApprox](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5193) ({ref}`ham <api-ham>`), {ref}`MIO_InputBlock <api-input-interface-inputblock>` (= `InputBlock`), {ref}`MIO_InputFindBlock <api-input-interface-inputfindblock>` (= `InputFindBlock`), {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_InputSearchLabel <api-input-inputsearchlabel>` (= `InputSearchLabel`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`num2str <api-mio-interface-num2str>`, [randomInRange](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4023) ({ref}`ham <api-ham>`), {ref}`RandSeedFromInput <api-random-randseedfrominput>`, {ref}`MIO_StringComp <api-string-interface-stringcomp>` (= `StringComp`), {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_PrintDebug <api-sys-sysprintdebug>` (= `SysPrintDebug`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`), {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

{ref}`CalcSelect <api-calc-calcselect>`, {ref}`HamInit <api-ham-haminit>`
