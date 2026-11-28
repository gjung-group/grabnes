---
edit_url: null
---

(api-diag-diagdos)=
# `DiagDOS`

```{index} single: DiagDOS (subroutine in diag)
```

```fortran
subroutine DiagDOS()
```

- **Module:** {ref}`diag <api-diag>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/diag.F90, line 178`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L178)

## Input keys read (5)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Diag.KGrid`
  - Number of k-points along each reciprocal lattice vector (three integers) of the k-grid of the DOS, projected DOS, 3-D bands and spectral cuts; Monkhorst-Pack points for the DOS.
* - `Kubo.Epsilon`
  - Diagonalisation DOS: width of the Gaussian broadening of every level \[eV\].
* - `Kubo.NumberofEnergyPoints`
  - Diagonalisation DOS: number of points of the energy grid from DOS.Emin to DOS.Emax.
* - `DOS.Emin`
  - Lower end of the energy grid of the diagonalisation DOS \[eV\].
* - `DOS.Emax`
  - Upper end of the energy grid of the diagonalisation DOS \[eV\].
```

## Calls

[DiagCheckThreads](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L584) ({ref}`diag <api-diag>`), [DiagHam](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L5128) ({ref}`diag <api-diag>`), [GetUnit](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L606) ({ref}`file <api-file>`), {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_PrintDebug <api-sys-sysprintdebug>` (= `SysPrintDebug`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`)

## Called by

[CalcDiag](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L678) ({ref}`calc <api-calc>`)
