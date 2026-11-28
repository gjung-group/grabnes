---
edit_url: null
---

(api-sys-syssubsection)=
# `SysSubsection`

```{index} single: SysSubsection (subroutine in sys)
```

```fortran
subroutine SysSubsection(title, module)
```

- **Module:** {ref}`sys <api-sys>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/sys.F90, line 391`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/sys.F90#L391)

Print a sub-section title inside a section: "module: --- title ---".

## Parameters

`title` : `character(len=*), intent(in)`
: —

`module` : `character(len=*), intent(in)`
: —

## Calls

{ref}`SysPrint <api-sys-interface-sysprint>`

## Called by

[CalculateChernTAPW](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L12066) ({ref}`diag <api-diag>`), {ref}`Diag3DBands <api-diag-diag3dbands>`, {ref}`DiagBands <api-diag-diagbands>`, {ref}`DiagBandsAroundK <api-diag-diagbandsaroundk>`, {ref}`DiagBandsG <api-diag-diagbandsg>`, {ref}`DiagBandsRashba <api-diag-diagbandsrashba>`, {ref}`DiagChern <api-diag-diagchern>`, {ref}`DiagDOS <api-diag-diagdos>`, [DiagH0TAPW](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L5902) ({ref}`diag <api-diag>`), {ref}`DiagPDOS <api-diag-diagpdos>`, {ref}`DiagSpectralFunction <api-diag-diagspectralfunction>`, {ref}`DiagSpectralFunctionKGrid <api-diag-diagspectralfunctionkgrid>`, {ref}`DiagSpectralFunctionKGridInequivalent <api-diag-diagspectralfunctionkgridinequivalent>`, {ref}`DiagSpectralFunctionKGridInequivalent_v2 <api-diag-diagspectralfunctionkgridinequivalent_v2>`, {ref}`DiagSpectralFunctionKGridInequivalentEnergyCut <api-diag-diagspectralfunctionkgridinequivalentenergycut>`, {ref}`DiagSpectralFunctionKGridInequivalentEnergyCut_v2 <api-diag-diagspectralfunctionkgridinequivalentenergycut_v2>`, {ref}`DiagSpectralFunctionKGridInequivalentEnergyCutNickDale <api-diag-diagspectralfunctionkgridinequivalentenergycutnickdale>`, {ref}`GaussHeight <api-gauss-gaussheight>`, {ref}`GaussPot <api-gauss-gausspot>`, {ref}`GaussPotDefinedPositions <api-gauss-gausspotdefinedpositions>`, {ref}`HamHopping <api-ham-hamhopping>`, {ref}`HamInit <api-ham-haminit>`, [HamMeasureLayerTwist](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4124) ({ref}`ham <api-ham>`), {ref}`HamOnSite <api-ham-hamonsite>`, {ref}`InterfacePot1 <api-interface-interfacepot1>`, {ref}`KuboTEvol <api-kubo-kubotevol>`
