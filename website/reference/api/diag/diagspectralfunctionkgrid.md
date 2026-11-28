---
edit_url: null
---

(api-diag-diagspectralfunctionkgrid)=
# `DiagSpectralFunctionKGrid`

```{index} single: DiagSpectralFunctionKGrid (subroutine in diag)
```

```fortran
subroutine DiagSpectralFunctionKGrid()
```

- **Module:** {ref}`diag <api-diag>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/diag.F90, line 2333`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2333)

## Input keys read (1)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Structure.CellSize`
  - Here: the graphene-lattice path is divided by CellSize to give the k of the cell.
```

## Calls

[convolve](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L12004) ({ref}`diag <api-diag>`), [DiagSpectralEpsilon](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2103) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKu](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11161) ({ref}`diag <api-diag>`), [DiagSpectralWindow](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2088) ({ref}`diag <api-diag>`), {ref}`MIO_InputBlock <api-input-interface-inputblock>` (= `InputBlock`), {ref}`MIO_InputFindBlock <api-input-interface-inputfindblock>` (= `InputFindBlock`), {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`)
