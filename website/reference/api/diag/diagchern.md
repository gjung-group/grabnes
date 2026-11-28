---
edit_url: null
---

(api-diag-diagchern)=
# `DiagChern`

```{index} single: DiagChern (subroutine in diag)
```

```fortran
subroutine DiagChern()
```

- **Module:** {ref}`diag <api-diag>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/diag.F90, line 1561`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L1561)

## Calls

[DiagCheckThreads](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L584) ({ref}`diag <api-diag>`), [DiagHamChern](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L10766) ({ref}`diag <api-diag>`), {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`)

## Called by

[CalcDiag](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L678) ({ref}`calc <api-calc>`)
