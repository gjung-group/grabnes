---
edit_url: null
---

(api-ham-haminit)=
# `HamInit`

```{index} single: HamInit (subroutine in ham)
```

```fortran
subroutine HamInit()
```

- **Module:** {ref}`ham <api-ham>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/ham.F90, line 242`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L242)

Initialize Hamiltonian parameters and neighbor lists.

Reads input parameters (strain, Zeeman terms, models), prepares neighbor lists and internal arrays, and sets module state. On-site terms are computed in HamOnSite; hopping is assembled elsewhere.

**Side effects**

- Initializes module-level variables
- Builds neighbor lists
- Configures strain and magnetic field parameters

## Input keys read (46)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Zeeman.Term`
  - Spin Zeeman shift +-Zeeman.Factor x Zeeman.Bmag on the spin sector Zeeman.Spin; with Zeeman.SpinPolarized both sectors. Follows SOC.LayerControl unless Zeeman.FollowSOCLayers .false.
* - `Zeeman.PseudoTerm`
  - Pseudo-Zeeman term (Zeeman.PseudoFactor): not implemented, no effect.
* - `Zeeman.Spin`
  - Zeeman.Term in a one-spin run: spin sector computed, 1 (up) or -1 / 2 (down).
* - `Zeeman.Factor`
  - Zeeman.Term: Zeeman energy per tesla \[eV/T\]; the shift is Factor x Zeeman.Bmag.
* - `Zeeman.PseudoFactor`
  - Zeeman.PseudoTerm: factor \[eV/T\] (the term is not implemented).
* - `Zeeman.SpinPolarized`
  - Two spin sectors (nspin = 2), needed by the spinful terms (Zeeman, SOC). Refused unless Hubbard.EnableSCF .false.
* - `SOC.Intrinsic`
  - Spin-independent sublattice mass +-SOC.LambdaI (not Kane-Mele, see SOC.KaneMele).
* - `SOC.LambdaI`
  - SOC.Intrinsic: amplitude \[eV\].
* - `SOC.Ising`
  - Ising (valley-Zeeman) spin-orbit term +-SOC.LambdaIsing.
* - `SOC.LambdaIsing`
  - SOC.Ising: amplitude \[eV\].
* - `SOC.Rashba`
  - Rashba spin-orbit term on the first-neighbour bonds (spinful block TAPW path).
* - `SOC.LambdaR`
  - SOC.Rashba: amplitude \[eV\]; the K-point splitting is 3 LambdaR (2/3 of the Gmitra lambda\_R).
* - `SOC.PIA`
  - Pseudospin-inversion-asymmetry SOC term (Gmitra), spinful TAPW only.
* - `SOC.LambdaPIA`
  - SOC.PIA: amplitude for both sublattices \[eV\].
* - `SOC.LambdaPIAA`
  - SOC.PIA: amplitude on sublattice A \[eV\] (default SOC.LambdaPIA).
* - `SOC.LambdaPIAB`
  - SOC.PIA: amplitude on sublattice B \[eV\] (default SOC.LambdaPIA).
* - `SOC.RashbaConsistentPhase`
  - .true. (default since 2026-10-16): spinful bond terms carry the same Bloch phase as the hopping. .false. restores the former Rashba phase, exact only at K and Gamma (warning printed).
* - `SOC.RashbaAngle`
  - SOC.Rashba: rotation of the Rashba spin texture (radial component) \[degrees\].
* - `SOC.KaneMele`
  - Kane-Mele intrinsic spin-orbit term on second-neighbour bonds, SOC.LambdaKM.
* - `SOC.LambdaKM`
  - SOC.KaneMele: amplitude lambda\_KM \[eV\] (gap 2 lambda\_KM in a single layer).
* - `SOC.Interlayer`
  - Interlayer spin-orbit term lambda\_4 between same-sublattice sites of adjacent untwisted graphene layers, i (lambda/3) c\_ij sgn(dz) (d x s)\_z.
* - `SOC.LambdaInterlayer`
  - SOC.Interlayer: amplitude \[eV\].
* - `Zeeman.FollowSOCLayers`
  - With SOC.LayerControl: the Zeeman term acts on the same layers as SOC (.false.: every layer).
* - `SOC.LayerControl`
  - Apply the SOC terms only to the layers listed in SOCLayers (e.g. 'SOCLayers 1 5'); otherwise every layer.
* - `Haldane.LayerControl`
  - Apply the Haldane term only to the layers listed in HaldaneLayers; otherwise see Haldane.NNN.
* - `Haldane.NNN`
  - Haldane second-neighbour term (Haldane.\* keys) on every layer.
* - `Neigh.fastNN`
  - Neighbour search fastNN: refused (it reads past the end of its bins).
* - `Neigh.fastNNnotsquare`
  - Binned neighbour search, the default and best tested one; .false. selects NeighList (see Neigh.CompareAll).
* - `Neigh.fastNNnotsquareNotRectangle`
  - Neighbour search fastNNnotsquareNotRectangle: refused (it reads past the end of its bins).
* - `Neigh.CutAtNN3`
  - Deprecated: limits the search to three intralayer shells; use TB.NeighLevels 3.
* - `Neigh.IntralayerRadius`
  - In-plane radius of the intralayer neighbour search \[Angstrom\], replacing the shell table's: for strained or corrugated cells (strained cells need 5.1134 to reproduce the search before October 2026).
* - `Intralayer.F2G2Model`
  - F2G2 intralayer model (5 shells); .false.: the 8-shell parameter set. A warning lists shells not searched.
* - `Neigh.LayerNeighbors`
  - Search interlayer neighbours (any value but 0); 0: no interlayer pairs.
* - `Stack.Bulk`
  - The cell is periodic along z too (graphite-like stack): interlayer pairs across the z boundary are kept. DiagHamSparse samples k\_z = 0 only.
* - `Stack.BulkSmall`
  - Stack.Bulk: neighbour search for small cells (fastNNnotsquareBulkSmall).
* - `Stack.NonBulkSmall`
  - Slab (non-periodic z) with the neighbour search for small cells.
* - `Structure.InterlayerDistance`
  - Here: read by the neighbour search set-up.
* - `Neigh.LayerDistFactor`
  - In-plane radius of the interlayer neighbour search: 1.1 aG/sqrt(3) times this factor.
* - `Output.ReadDataFiles`
  - Here: the on-site energies are read from generate.e (units of g0) instead of computed.
* - `Moire.Strain`
  - Modify the hoppings by a model moire strain (Moire.StrainFactor, Moire.HeightAmp).
* - `Moire.StrainFactor`
  - Moire.Strain: beta of the hopping law t = exp(-beta (d/d\_CC - 1)).
* - `Moire.HeightAmp`
  - Moire.Strain: amplitude of the height modulation \[Angstrom\].
* - `Strain.RandomStrain`
  - Random Gaussian height bumps (Strain.GaussHeight\*): refused at present.
* - `Output.PrintBubble`
  - Strain bubble diagnostic around a hard-coded atom: refused at present.
* - `Strain.BubbleSigmaR`
  - Output.PrintBubble (refused): Gaussian width of the height bump \[Angstrom\].
* - `Strain.BubbleHeight`
  - Output.PrintBubble (refused): height of the Gaussian bump \[Angstrom\].
```

## See also

{ref}`HamOnSite <api-ham-hamonsite>`

## Calls

{ref}`AtomsSetCart <api-atoms-atomssetcart>`, {ref}`AtomsSetFrac <api-atoms-atomssetfrac>`, {ref}`fastNN <api-neigh-fastnn>`, {ref}`fastNNnotsquare <api-neigh-fastnnnotsquare>`, {ref}`fastNNnotsquareBulk <api-neigh-fastnnnotsquarebulk>`, {ref}`fastNNnotsquareBulkSmall <api-neigh-fastnnnotsquarebulksmall>`, {ref}`fastNNnotsquareNotRectangle <api-neigh-fastnnnotsquarenotrectangle>`, {ref}`GaussHeight <api-gauss-gaussheight>`, [HamDisplacementsFromPositions](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L13132) ({ref}`ham <api-ham>`), {ref}`HamOnSite <api-ham-hamonsite>`, {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_InputParameterVariable <api-input-interface-inputparametervariable>` (= `InputParameterVariable`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`NeighList <api-neigh-neighlist>`, {ref}`NeighShellCutoff2 <api-neigh-neighshellcutoff2>`, {ref}`NeighShellWarn2 <api-neigh-neighshellwarn2>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_PrintDebug <api-sys-sysprintdebug>` (= `SysPrintDebug`), {ref}`MIO_Section <api-sys-syssection>` (= `SysSection`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`), {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

{ref}`CalcSelect <api-calc-calcselect>`
