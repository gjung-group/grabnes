---
edit_url: null
---

(api-magf-magfinit)=
# `MagfInit`

```{index} single: MagfInit (subroutine in magf)
```

```fortran
subroutine MagfInit()
```

- **Module:** {ref}`magf <api-magf>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/magf.F90, line 20`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/magf.F90#L20)

## Input keys read (8)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `MagField.Value`
  - Perpendicular magnetic field \[T\], rounded to an integer number of flux quanta through the cell (a warning is printed when it rounds to zero). Peierls phases on every hopping.
* - `MagField.Integer`
  - Perpendicular magnetic field as m flux quanta through the whole cell: B = m Phi\_0 / A. Use MagField.Integer 0 for no field.
* - `MagField.InitialB`
  - Field sweep from MagField.InitialB to MagField.FinalB \[T\], each rounded to whole flux quanta.
* - `MagField.FinalB`
  - Last field of the sweep started by MagField.InitialB \[T\].
* - `MagField.InitialInteger`
  - Field sweep in flux quanta through the cell, from InitialInteger to FinalInteger: one calculation per value, in folders m\<m>\_B\<B>. Both 0 = no field.
* - `MagField.FinalInteger`
  - Last value of the MagField.InitialInteger sweep (flux quanta through the cell).
* - `MagField.Step`
  - Step of a field sweep \[T\], rounded to whole flux quanta (otherwise MagField.IntegerStep).
* - `MagField.IntegerStep`
  - Step of a field sweep, in flux quanta through the cell.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_InputSearchLabel <api-input-inputsearchlabel>` (= `InputSearchLabel`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

{ref}`CalcSelect <api-calc-calcselect>`
