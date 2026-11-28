---
edit_url: null
---

(api-diag-diagspectralfunctionkgridinequivalent)=
# `DiagSpectralFunctionKGridInequivalent`

```{index} single: DiagSpectralFunctionKGridInequivalent (subroutine in diag)
```

```fortran
subroutine DiagSpectralFunctionKGridInequivalent()
```

- **Module:** {ref}`diag <api-diag>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/diag.F90, line 2572`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2572)

## Input keys read (14)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Spectral.WeiKu`
  - Unfolded spectral weight per layer and sublattice on a graphene reference lattice (Wei Ku method), files \<prefix>.spectral\<layer>\<A\|B>.
* - `Spectral.UseGaussianBroadening`
  - Spectral.WeiKu: spread each level over the energy grid with a Gaussian of width Kubo.Epsilon.
* - `Spectral.WeiKuOld`
  - Spectral.WeiKu: use the older weight routine (A and B weights only).
* - `Spectral.Nishi`
  - Unfolding weights of the Nishi routine (sub-cell translations).
* - `Spectral.Lee`
  - Unfolding weights of the Lee routine.
* - `Spectral.FoldByOne`
  - Fold the graphene-lattice k into the cell's zone by subtracting G10 + G01.
* - `Spectral.FoldByZero`
  - Do not fold the graphene-lattice k (use it as it is).
* - `Spectral.UseCoordinates`
  - Spectral.Path points in Cartesian coordinates \[1/Angstrom\] instead of reciprocal-lattice units.
* - `Structure.MoireCellParameters`
  - Here: sets the angle of Spectral.RotateReferenceSystem.
* - `Spectral.RotateReferenceSystem`
  - Rotate the graphene reference lattice by the angle of Structure.MoireCellParameters.
* - `Spectral.AlignReferenceSystem`
  - Rotate the graphene reference lattice by Spectral.AlignmentAngle (a layer unfolds under its own orientation).
* - `Spectral.AlignmentAngle`
  - Spectral.AlignReferenceSystem: orientation of the reference lattice \[degrees\].
* - `Spectral.topBottomRatio`
  - Passed to the unfolding weight routines, which do not use it at present.
* - `Spectral.energyGridResolution`
  - A level adds its weight to every energy point closer than this \[eV\]: points further apart than twice it lose levels (a warning is printed).
```

## Calls

[convolve](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L12004) ({ref}`diag <api-diag>`), [DiagSpectralCheckGrid](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2117) ({ref}`diag <api-diag>`), [DiagSpectralEpsilon](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2103) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalentLee](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11693) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalentMoreOrbitals](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11542) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalentNishi](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11866) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalentOld](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11286) ({ref}`diag <api-diag>`), [DiagSpectralWindow](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2088) ({ref}`diag <api-diag>`), {ref}`MIO_InputBlock <api-input-interface-inputblock>` (= `InputBlock`), {ref}`MIO_InputFindBlock <api-input-interface-inputfindblock>` (= `InputFindBlock`), {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_InputSearchLabel <api-input-inputsearchlabel>` (= `InputSearchLabel`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_PrintDebug <api-sys-sysprintdebug>` (= `SysPrintDebug`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`)

## Called by

[CalcDiag](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L678) ({ref}`calc <api-calc>`)
