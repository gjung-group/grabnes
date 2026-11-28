---
edit_url: null
---

(api-diag)=
# `diag`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/diag.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`atoms <api-atoms>`, {ref}`cell <api-cell>`, {ref}`constants <api-constants>`, {ref}`ham <api-ham>`, {ref}`interface <api-interface>`, `iso_fortran_env`, {ref}`magf <api-magf>`, {ref}`math <api-math>`, {ref}`mio <api-mio>`, {ref}`name <api-name>`, {ref}`neigh <api-neigh>`, `omp_lib`, {ref}`scf <api-scf>`, {ref}`tbpar <api-tbpar>`
- **Used by:** {ref}`calc <api-calc>`
- **Public interface:** 16 procedures, 32 variables
- **Input keys:** 63 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`Diag3DBands <api-diag-diag3dbands>`
  - subroutine
  - —
* - {ref}`DiagBands <api-diag-diagbands>`
  - subroutine
  - —
* - {ref}`DiagBandsAroundK <api-diag-diagbandsaroundk>`
  - subroutine
  - —
* - {ref}`DiagBandsG <api-diag-diagbandsg>`
  - subroutine
  - —
* - {ref}`DiagBandsRashba <api-diag-diagbandsrashba>`
  - subroutine
  - —
* - {ref}`DiagChern <api-diag-diagchern>`
  - subroutine
  - —
* - {ref}`DiagDOS <api-diag-diagdos>`
  - subroutine
  - —
* - {ref}`DiagInit <api-diag-diaginit>`
  - subroutine
  - —
* - {ref}`DiagPDOS <api-diag-diagpdos>`
  - subroutine
  - —
* - {ref}`DiagSpectralFunction <api-diag-diagspectralfunction>`
  - subroutine
  - —
* - {ref}`DiagSpectralFunctionKGrid <api-diag-diagspectralfunctionkgrid>`
  - subroutine
  - —
* - {ref}`DiagSpectralFunctionKGridInequivalent <api-diag-diagspectralfunctionkgridinequivalent>`
  - subroutine
  - —
* - {ref}`DiagSpectralFunctionKGridInequivalent_v2 <api-diag-diagspectralfunctionkgridinequivalent_v2>`
  - subroutine
  - —
* - {ref}`DiagSpectralFunctionKGridInequivalentEnergyCut <api-diag-diagspectralfunctionkgridinequivalentenergycut>`
  - subroutine
  - —
* - {ref}`DiagSpectralFunctionKGridInequivalentEnergyCut_v2 <api-diag-diagspectralfunctionkgridinequivalentenergycut_v2>`
  - subroutine
  - —
* - {ref}`DiagSpectralFunctionKGridInequivalentEnergyCutNickDale <api-diag-diagspectralfunctionkgridinequivalentenergycutnickdale>`
  - subroutine
  - —
```

## Module variables

```{list-table}
:header-rows: 1

* - Name
  - Declaration
  - Initial value
  - Description
* - `calculate3DTAPWBands`
  - `logical, save`
  - `.false.`
  - —
* - `calculateChern`
  - `logical, save`
  - `.false.`
  - —
* - `checkTAPWUnitary`
  - `logical, save`
  - `.false.`
  - —
* - `fermi_energy`
  - `real(dp), save`
  - `0.0_dp`
  - —
* - `forceBlockTAPW`
  - `logical, save`
  - `.false.`
  - —
* - `gammaCentred3D`
  - `logical, save`
  - `.true.`
  - —
* - `gGridRotationAngle`
  - `real(dp), save`
  - `0.0_dp`
  - —
* - `gw_b1`
  - `integer, save`
  - `0`
  - —
* - `gw_b2`
  - `integer, save`
  - `0`
  - —
* - `gw_ntop`
  - `integer, save`
  - `8`
  - —
* - `gWeightsHam`
  - `logical, save`
  - `.false.`
  - Diag.GWeightsHam (opt-in, needs Diag.GWeights): the projected TAPW Hamiltonian before diagonalisation -> \<prefix>.TAPWHam, for offline first-order (two-level) analysis of avoided crossings.
* - `gWeightsTAPW`
  - `logical, save`
  - `.false.`
  - Plane-wave composition of TAPW eigenstates (opt-in diagnostic, Diag.GWeights) -> \<prefix>.GWeights: per k and band in gw\_b1..gw\_b2, the gw\_ntop largest weights sum\_label \|c\_(G,label)\|^2.
* - `layerWeightsTAPW`
  - `logical, save`
  - `.false.`
  - Layer / sublattice composition of TAPW eigenstates (opt-in, Diag.LayerWeights) -> \<prefix>.LayerWeights: per k and band in lw\_b1..lw\_b2 (0 0 = all), the atom-space weight sum\_{i in label} \|psi\_i\|^2, psi = X c.
* - `lw_b1`
  - `integer, save`
  - `0`
  - —
* - `lw_b2`
  - `integer, save`
  - `0`
  - —
* - `M_tapw`
  - `integer, save`
  - `0`
  - —
* - `moireAngle`
  - `real(dp), save`
  - `33.004491598883078_dp`
  - —
* - `nk_3D_x`
  - `integer, save`
  - `10`
  - —
* - `nk_3D_y`
  - `integer, save`
  - `10`
  - —
* - `nk_chern_x`
  - `integer, save`
  - `10`
  - —
* - `nk_chern_y`
  - `integer, save`
  - `10`
  - —
* - `physicalTwistAngle`
  - `real(dp), save`
  - `1.08_dp`
  - —
* - `skipGRotation`
  - `logical, save`
  - `.false.`
  - —
* - `socDebug`
  - `logical, save`
  - `.false.`
  - —
* - `tapwBothValleys`
  - `logical, save`
  - `.false.`
  - Both valleys in one TAPW basis (opt-in): the G-list is the union of the shells around K and K', so X^H H X also carries the intervalley block.  Default .false. = the single-valley basis.
* - `tapwDebug`
  - `logical, save`
  - `.false.`
  - —
* - `tapwLowdin`
  - `logical, save`
  - `.false.`
  - —
* - `tapwNG`
  - `integer, save`
  - `180`
  - —
* - `tapwValleyDecouple`
  - `logical, save`
  - `.false.`
  - —
* - `useKprimeValley`
  - `logical, save`
  - `.false.`
  - —
* - `useRigidPositions`
  - `logical, save`
  - `.false.`
  - —
* - `useTriangularTruncation`
  - `logical, save`
  - `.false.`
  - —
```

## Input keys read by private procedures

Procedures that are not part of the public interface but read input keys.

**`DiagSpectralWindow`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2088))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Spectral.Emin`
  - Lower end of the energy grid of the spectral function \[eV\] (units of g0 before November 2026).
* - `Spectral.Emax`
  - Upper end of the energy grid of the spectral function \[eV\] (units of g0 before November 2026).
```

**`DiagSpectralEpsilon`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2103))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Kubo.Epsilon`
  - Spectral function: Gaussian width of Spectral.UseGaussianBroadening / GaussianConvolution \[eV\] (units of g0 before November 2026).
```

**`DiagHamSparse`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L5180))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Diag.SparseUseShift`
  - Diag.SparseSolver: shift-invert mode (MKL PARDISO), the Bands.SparseNeig levels nearest Diag.SparseSetShift; output \<prefix>.SparseEig.
* - `Diag.SparseSetShift`
  - Diag.SparseUseShift: energy around which the levels are computed \[eV\].
* - `Diag.SparseSaveRitz`
  - Diag.SparseSolver: also compute the eigenvectors (Ritz vectors) in ARPACK.
* - `Diag.SparseCount`
  - Diag.SparseUseShift: count exactly the states below the shift (so level indices are absolute) and check completeness ('complete T/F' per k in \<prefix>.SparseEig).
* - `Diag.SparseSetTol`
  - Diag.SparseSolver: ARPACK tolerance (default 1e-9 in shift-invert mode, 0.1 otherwise).
```

**`DiagH0TAPW`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L5902))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `TAPW.UseTriangularDistanceOrder`
  - Diag.TriangularTruncation: order the triangular G list by distance (default) instead of grid order.
```

**`DiagHamChern`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L10766))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Kubo.Epsilon`
  - Here (DiagHamChern, refused): regularisation of the energy denominators.
```

**`DiagSpectralWeightNishi`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11043))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Spectral.numberOfTS`
  - Spectral.Nishi: number of sub-cell translations (set to 6 by the routine).
```

**`CalculateChernTAPW`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L12066))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Diag.ChernEnergyWindow`
  - Diag.ChernCalculation: energy window around Diag.Chern\_EF for the band selection \[eV\].
* - `Diag.ChernAllOccupiedBands`
  - Diag.ChernCalculation: sum over all bands below the Fermi energy.
* - `Diag.ChernNBands`
  - Diag.ChernCalculation: number of bands included.
* - `Diag.ChernUseOptionA`
  - Diag.ChernCalculation: Berry curvature from the velocities in the TB basis instead of the TAPW basis.
* - `Diag.ChernBandsBelowFermi`
  - Diag.ChernCalculation: bands below the Fermi energy included.
* - `Diag.ChernBandsAboveFermi`
  - Diag.ChernCalculation: bands above the Fermi energy included.
* - `Diag.ChernTargetBands`
  - Diag.ChernCalculation: number of bands closest to the Fermi energy whose Chern numbers are reported.
```

**`CalculateBerryAtKpointFromStored`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L12802))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Diag.ChernEnergyThreshold`
  - Diag.ChernCalculation: pairs of levels closer than this \[units of g0\] are left out of the Kubo sum.
* - `Diag.ChernEtaEV`
  - Diag.ChernCalculation: read but not used at present.
* - `Diag.ChernEps`
  - Diag.ChernCalculation: eps in the denominator (E\_n - E\_m)^2 + eps of the Kubo sum \[units of g0^2\].
```
