---
edit_url: null
---

(api-calc)=
# `calc`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/calc.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`atoms <api-atoms>`, {ref}`cell <api-cell>`, {ref}`constants <api-constants>`, {ref}`diag <api-diag>`, {ref}`ham <api-ham>`, {ref}`kubo <api-kubo>`, {ref}`kuboarrays <api-kuboarrays>`, {ref}`kubosubs <api-kubosubs>`, {ref}`magf <api-magf>`, {ref}`math <api-math>`, {ref}`mio <api-mio>`, {ref}`moireBLShift <api-moireblshift>`, {ref}`name <api-name>`, {ref}`neigh <api-neigh>`, {ref}`scf <api-scf>`, {ref}`tbpar <api-tbpar>`
- **Used by:** {ref}`grabnes <api-grabnes>`
- **Public interface:** 1 procedure
- **Input keys:** 93 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`CalcSelect <api-calc-calcselect>`
  - subroutine
  - —
```

## Input keys read by private procedures

Procedures that are not part of the public interface but read input keys.

**`CalcKubo`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Kubo.RecursionNumber`
  - Number of Lanczos steps (length of the continued fraction); sets the energy resolution. No reorthogonalisation: the useful length is bounded by the number of distinct levels the start vector sees.
* - `Kubo.NumberofTimeSteps`
  - Number of time steps of the Kubo time evolution (KuboTEvol, only without Calculate.OnlyDOS).
* - `Kubo.TimeStep`
  - Time step of the Kubo time evolution, in units of hbar/g0 (g0 = TB.Hopping).
* - `Kubo.NumberofPolynomials`
  - Number of Chebyshev polynomials in the expansion of U(Kubo.TimeStep) (KuboCn; \<prefix>.cn.dat).
* - `Kubo.NumberofEnergyPoints`
  - Number of points of the energy grid from Kubo.EnergyMin to Kubo.EnergyMax at which the DOS is written.
* - `Kubo.Epsilon`
  - Recursion DOS: broadening eta of the continued fraction, the Green function is taken at E + i eta \[eV\].
* - `Kubo.WriteNumberofSteps`
  - The time evolution writes its output (\<prefix>.Z2\_\<step>.dat) every this many time steps.
* - `Kubo.EnergyMin`
  - Lower end of the energy grid of the recursion DOS \[eV\].
* - `Kubo.EnergyMax`
  - Upper end of the energy grid of the recursion DOS \[eV\].
* - `Kubo.PDOS`
  - Local DOS of selected atoms, layers or species instead of the DOS of the whole cell. The selection is the FIRST of PDOSLayer, PDOSLayerAndSpecies, PDOSDist, PDOSMoireSC1, PDOSMoireSuperMoire, ...2, ...3, PDOSAtomList, PDOSByNumber, PDOSPNP that is set; PDOSMoireSuperMoire(2,3) default to .true. Output in folders PDOS\_\*.
* - `Kubo.PDOSxmin`
  - Kubo.PDOSDist: lower x bound of the selected atoms \[Angstrom\].
* - `Kubo.PDOSxmax`
  - Kubo.PDOSDist: upper x bound of the selected atoms \[Angstrom\].
* - `Kubo.PDOSymin`
  - Kubo.PDOSDist: lower y bound of the selected atoms \[Angstrom\].
* - `Kubo.PDOSymax`
  - Kubo.PDOSDist: upper y bound of the selected atoms \[Angstrom\].
* - `Kubo.PDOSDist`
  - Local DOS of every atom with PDOSxmin \< x \< PDOSxmax, PDOSymin \< y \< PDOSymax and index >= PDOSMinValue, one folder PDOS\_atom\<index> each.
* - `Kubo.PDOSLayer`
  - DOS of each layer 1..Kubo.NumberOfLayers, from a random-phase vector on the atoms of that layer (folders PDOS\_layer\<NN>).
* - `Kubo.PDOSLayerAndSpecies`
  - DOS of each layer and Species, from a random-phase vector on those atoms (folders PDOS\_layer\<NN>\_sub\<NN>; combinations without atoms are skipped).
* - `Kubo.PDOSMoireSC1`
  - Local DOS of every atom of the first SuperCell copy (both fractional coordinates below 1/SuperCell) with index >= PDOSMinValue, one folder PDOS\_atom\<index> each.
* - `Kubo.PDOSMoireSuperMoire`
  - Local DOS of the Species-1 atoms whose distance along both lattice vectors (fractional coordinate times CellSize\*SuperCell\*aG) is below Kubo.PDOSMoireSuperMoireLength. Default .true.
* - `Kubo.PDOSMoireSuperMoire2`
  - As Kubo.PDOSMoireSuperMoire, for distances from 1 to 2 times Kubo.PDOSMoireSuperMoireLength.
* - `Kubo.PDOSMoireSuperMoire3`
  - As Kubo.PDOSMoireSuperMoire, for distances from 2 to 3 times Kubo.PDOSMoireSuperMoireLength.
* - `Kubo.PDOSMoireSuperMoireLength`
  - Side of the regions selected by Kubo.PDOSMoireSuperMoire(2,3) \[Angstrom\].
* - `Kubo.PDOSMinValue`
  - Atoms with a smaller index are skipped by Kubo.PDOSDist, PDOSMoireSC1 and PDOSMoireSuperMoire(2,3).
* - `Structure.CellSize`
  - Here: sets the cell length CellSize\*SuperCell\*aG used by the Kubo.PDOSMoireSuperMoire selections.
* - `Kubo.PDOSAtomList`
  - Local DOS of the six atoms listed in Kubo.PDOSList.
* - `Kubo.PDOSByNumber`
  - Local DOS of the single atom Kubo.PDOSNumber.
* - `Kubo.PDOSIgnoreLayer1and4`
  - With Kubo.PDOSMoireSC1: skip the atoms of layers 1 and 4.
* - `Kubo.NumberOfLayers`
  - Number of layers of Kubo.PDOSLayer / PDOSLayerAndSpecies (layers 1..N; a warning is printed when the structure has more).
* - `Kubo.NumberOfBNAtoms1`
  - Kubo.PDOSLayer with Stack.EncapsulatedFourLayers: assumed atoms of layer 1 per input cell (times SuperCell^2). The start vector is normalised afterwards, so the value does not change the result.
* - `Kubo.NumberOfBNAtoms2`
  - As Kubo.NumberOfBNAtoms1, for layer 4.
* - `Kubo.NumberOfCAtoms`
  - As Kubo.NumberOfBNAtoms1, for layers 2 and 3 (Stack.EncapsulatedFourLayers) or layers 1-6 (Stack.T3BG).
* - `Stack.EncapsulatedFourLayers`
  - Stack model flag (ham.F90). Here: Kubo.PDOSLayer takes the layer sizes from Kubo.NumberOf\*Atoms.
* - `Stack.EncapsulatedSixLayers`
  - Stack model flag (ham.F90); read here but not used in this routine.
* - `Stack.T3BG`
  - Stack model flag (ham.F90). Here: Kubo.PDOSLayer takes the layer sizes from Kubo.NumberOfCAtoms.
* - `Kubo.PDOSPNP`
  - Local DOS of two atoms of a p-n-p potential (Potential.PNPKink): the first with x just above L/4 and the first with x just above L/4 + 2\*Potential.PNPDelta, L = CellSize\*SuperCell\*aG.
* - `Kubo.PDOSList`
  - Kubo.PDOSAtomList: exactly six atom indices (1-based, the order of generate.pos).
* - `Kubo.PDOSNumber`
  - Kubo.PDOSByNumber: index of the atom (1-based, the order of generate.pos).
* - `Potential.PNPDelta`
  - Width of the tanh steps of the p-n-p potential \[Angstrom\]. Here: offset of the second atom of Kubo.PDOSPNP (L/4 + 2\*PNPDelta).
* - `Calculate.Polynomials`
  - Compute the Chebyshev coefficients of U(Kubo.TimeStep) (\<prefix>.cn.dat) also with Calculate.OnlyDOS.
* - `Calculate.OnlyDOS`
  - Stop after the recursion DOS. Without it the Kubo time evolution (KuboTEvol) follows, which costs far more.
```

**`CalcTunn`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L595))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Tunn.BigKVecX`
  - First component of the k of Tunn.Calc, in units of the reciprocal lattice vectors (as in Bands.Path).
* - `Tunn.BigKVecY`
  - Second component of the k of Tunn.Calc, in units of the reciprocal lattice vectors.
* - `Tunn.BigKVecZ`
  - Third component of the k of Tunn.Calc, in units of the reciprocal lattice vectors.
* - `Tunn.NumberOfMoires`
  - Tunn.Calc: the written sums are divided by this number.
```

**`CalcDiag`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L678))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Calculate.DOS`
  - DOS by exact diagonalisation on a k-grid (DiagDOS; DiagPDOS with Calculate.PDOS).
* - `Calculate.PDOS`
  - With Calculate.DOS: projected DOS (DiagPDOS) instead of the total DOS.
* - `Calculate.Bands`
  - Band structure along Bands.Path (DiagBands, writes generate.bands). Bands.GrapheneUnitCell, Bands.AroundGrapheneK and Diag.RashbaCalculation select other band routines.
* - `Calculate.Chern`
  - Refused: this dense finite-difference route never computed a Chern number. Use Diag.ChernCalculation.
* - `Diag.RashbaCalculation`
  - With Calculate.Bands: DiagBandsRashba instead of DiagBands (refused at present, the routine needs updating).
* - `Diag.MoireAngle`
  - TAPW: orientation of the graphene lattice in the frame of the cell \[degrees\] (the canonical lattice, a1 along x, is rotated by it). 0 for a cell aligned with graphene; GenMoire prints it as 'Angle of commensurate cell'. A wrong value gives a smooth but empty low-energy spectrum.
* - `Diag.GGridRotationAngle`
  - TAPW: rotate the plane-wave (G) list by this angle after it is built \[degrees\]; 0 = no rotation.
* - `Diag.N_G`
  - TAPW: number of distinct G lengths (shells of the cell's reciprocal lattice, around the valley) in the basis; M = (number of G) x (layer, sublattice labels). N\_G 8 gives 55 G on a hexagonal cell; the default 180 is very large.
* - `Diag.ChernCalculation`
  - Chern numbers of the TAPW bands from the Kubo formula on a Diag.ChernGridX x ChernGridY k-grid.
* - `Diag.ChernGridX`
  - Diag.ChernCalculation: number of k-points along the first reciprocal lattice vector.
* - `Diag.ChernGridY`
  - Diag.ChernCalculation: number of k-points along the second reciprocal lattice vector.
* - `Diag.Chern_EF`
  - Diag.ChernCalculation: Fermi energy \[eV\].
* - `Diag.TriangularTruncation`
  - TAPW: G list from a triangular truncation (generate\_triangular\_G\_list) instead of shells of equal G length.
* - `Diag.CheckTAPWUnitary`
  - TAPW diagnostic: G list covering the whole first Brillouin zone (sized from Diag.PhysicalTwistAngle), and a check that the transformation is unitary (run with Diag.TAPWDebug).
* - `Diag.PhysicalTwistAngle`
  - Twist angle \[degrees\] used only to size the G list of Diag.CheckTAPWUnitary.
* - `Diag.UseKprimeValley`
  - TAPW basis around K' = (1/3,2/3) of the graphene reciprocal lattice instead of K = (2/3,1/3).
* - `Diag.TAPWDebug`
  - Extra TAPW diagnostics in the output (Hermiticity, G list, ...); the bands are unchanged.
* - `Diag.socDebug`
  - Verbose output of the spinful (SOC) block TAPW path.
* - `Diag.forceBlockTAPW`
  - Test switch: use the spinful block TAPW path even without Rashba SOC.
* - `Diag.UseRigidPositions`
  - TAPW: build the plane waves on the rigid positions of generateInit.xyz (fixed name, same atom order); the Hamiltonian keeps the relaxed positions. For relaxed cells Diag.TAPWLowdin is the usual choice.
* - `Diag.TAPWLowdin`
  - TAPW: orthonormalise the basis, X \<- X (X^H X)^(-1/2). Plane waves on displaced (relaxed) atoms are not orthonormal, and the plain projection then shifts the levels; prints the range of eig(X^H X).
* - `Diag.TAPWBothValleys`
  - TAPW: G shells around both K and K' in one basis (M doubles), so that intervalley coupling enters. Not with Diag.TriangularTruncation, CheckTAPWUnitary or the sparse TAPW path.
* - `Diag.TAPWValleyDecouple`
  - Diagnostic for Diag.TAPWBothValleys: set the K-K' block to zero.
* - `Diag.Calculate3DTAPWBands`
  - Bands on a Diag.3DBandsGridX x 3DBandsGridY k-grid instead of Bands.Path (DiagBands); rows kx ky kz E\_1 ... E\_M in generate.bands.
* - `Diag.GammaCentred3DGrid`
  - Diag.Calculate3DTAPWBands: grid through Gamma, so that high-symmetry points on it are sampled. .false. = the Monkhorst-Pack grid used before, offset by half a step.
* - `Diag.3DBandsGridX`
  - Diag.Calculate3DTAPWBands: number of k-points along the first reciprocal lattice vector.
* - `Diag.3DBandsGridY`
  - Diag.Calculate3DTAPWBands: number of k-points along the second reciprocal lattice vector.
* - `Diag.GWeights`
  - TAPW: write the plane-wave (G) composition of the eigenstates to \<prefix>.GWeights: per k and band, the energy and the Diag.GWeightsTop largest weights.
* - `Diag.GWeightsBandMin`
  - Diag.GWeights: first band written (1-based).
* - `Diag.GWeightsBandMax`
  - Diag.GWeights: last band written (1-based).
* - `Diag.GWeightsTop`
  - Diag.GWeights: number of largest plane-wave weights written per state.
* - `Diag.GWeightsHam`
  - With Diag.GWeights: also write the projected TAPW Hamiltonian H(k) \[eV\], before diagonalisation, to \<prefix>.TAPWHam (M x M per k: large).
* - `Diag.LayerWeights`
  - TAPW: write the weight of each eigenstate on each TAPW label (layer and Species), summed over the atoms, to \<prefix>.LayerWeights.
* - `Diag.LayerWeightsBandMin`
  - Diag.LayerWeights: first band written (1-based; 0 with LayerWeightsBandMax 0 = all bands).
* - `Diag.LayerWeightsBandMax`
  - Diag.LayerWeights: last band written (1-based; 0 with LayerWeightsBandMin 0 = all bands).
* - `Calculate.3DBands`
  - Bands on the k-grid Diag.KGrid by exact diagonalisation (Diag3DBands).
* - `Calculate.Spectral`
  - Unfolded spectral function along Bands.Path (Spectral.\* keys).
* - `Calculate.SpectralEnergyCut`
  - Unfolded spectral function on a k-grid at fixed energies (constant-energy cuts).
* - `Spectral.DirectBandVal`
  - Use the \_v2 routines for Calculate.Spectral / SpectralEnergyCut, which also write the band energy next to the unfolded weights.
* - `Calculate.SpectralEnergyCutNickDale`
  - Constant-energy cuts of the unfolded spectral function, NickDale variant of Calculate.SpectralEnergyCut.
* - `Bands.GrapheneUnitCell`
  - Bands with Bloch phases from the graphene primitive cell (DiagBandsG); refused unless the simulation cell is that cell.
* - `Bands.AroundGrapheneK`
  - Bands.Path points are taken relative to graphene K = (2/3,1/3) of the graphene reciprocal lattice.
* - `Hubbard.EnableSCF`
  - Self-consistent Hubbard term in two-spin runs. The module is unfinished: a two-spin diagonalisation is refused unless this is .false. (the default is .true.).
```
