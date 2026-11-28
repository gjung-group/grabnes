---
edit_url: null
---

(api-diag-diagspectralfunctionkgridinequivalentenergycutnickdale)=
# `DiagSpectralFunctionKGridInequivalentEnergyCutNickDale`

```{index} single: DiagSpectralFunctionKGridInequivalentEnergyCutNickDale (subroutine in diag)
```

```fortran
subroutine DiagSpectralFunctionKGridInequivalentEnergyCutNickDale()
```

- **Module:** {ref}`diag <api-diag>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/diag.F90, line 4698`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L4698)

## Input keys read (3)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Spectral.KGridCutX`
  - NickDale cut: length of the kx line from graphene K \[1/Angstrom\].
* - `Spectral.KGridCutY`
  - NickDale cut: the line is at ky = K\_y + KGridCutY \[1/Angstrom\].
* - `Spectral.KGridLowerGridHalf`
  - NickDale cut: the kx line runs from K\_x - KGridCutX to K\_x instead of K\_x to K\_x + KGridCutX.
```

## Calls

[convolve](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L12004) ({ref}`diag <api-diag>`), [DiagSpectralCheckGrid](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2117) ({ref}`diag <api-diag>`), [DiagSpectralEpsilon](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2103) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalent](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11412) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalentNishi](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11866) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalentOld](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11286) ({ref}`diag <api-diag>`), [DiagSpectralWindow](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2088) ({ref}`diag <api-diag>`), {ref}`MIO_InputBlock <api-input-interface-inputblock>` (= `InputBlock`), {ref}`MIO_InputFindBlock <api-input-interface-inputfindblock>` (= `InputFindBlock`), {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_InputSearchLabel <api-input-inputsearchlabel>` (= `InputSearchLabel`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_PrintDebug <api-sys-sysprintdebug>` (= `SysPrintDebug`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`), {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

[CalcDiag](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L678) ({ref}`calc <api-calc>`)
