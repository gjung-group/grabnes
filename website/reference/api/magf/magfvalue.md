---
edit_url: null
---

(api-magf-magfvalue)=
# `MagfValue`

```{index} single: MagfValue (subroutine in magf)
```

```fortran
subroutine MagfValue(mB)
```

- **Module:** {ref}`magf <api-magf>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/magf.F90, line 100`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/magf.F90#L100)

## Parameters

`mB` : `integer, intent(in)`
: —

## Input keys read (2)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Zeeman.Bmag`
  - Field of the Zeeman term \[T\], separate from the orbital field; the Zeeman energy is ZeemanFactor x Zeeman.Bmag.
* - `Run.TypeOfSystem`
  - Here: used to name the output folders of a field sweep.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Section <api-sys-syssection>` (= `SysSection`)

## Called by

{ref}`CalcSelect <api-calc-calcselect>`
