---
edit_url: null
---

(api-moireblshift-moireblshiftvalue)=
# `moireBLShiftValue`

```{index} single: moireBLShiftValue (subroutine in moireBLShift)
```

```fortran
subroutine moireBLShiftValue(mS)
```

- **Module:** {ref}`moireBLShift <api-moireblshift>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/moireBLShift.F90, line 55`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/moireBLShift.F90#L55)

## Parameters

`mS` : `integer, intent(in)`
: —

## Input keys read (12)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Moire.BilayerBottomXShift`
  - Moire.BLShiftDefinedValues: x shift of the bottom layer (written into the output prefix).
* - `Moire.BilayerBottomYShift`
  - Moire.BLShiftDefinedValues: y shift of the bottom layer.
* - `Moire.BilayerTopXShift`
  - Moire.BLShiftDefinedValues: x shift of the top layer.
* - `Moire.BilayerTopYShift`
  - Moire.BLShiftDefinedValues: y shift of the top layer.
* - `Moire.BilayerBottomXMin`
  - Moire.BLShiftRange: first x shift of the bottom layer.
* - `Moire.BilayerBottomXMax`
  - Moire.BLShiftRange: last x shift of the bottom layer.
* - `Moire.BilayerBottomYMin`
  - Moire.BLShiftRange: first y shift of the bottom layer.
* - `Moire.BilayerBottomYMax`
  - Moire.BLShiftRange: last y shift of the bottom layer.
* - `Moire.BilayerTopXMin`
  - Moire.BLShiftRange: first x shift of the top layer.
* - `Moire.BilayerTopXMax`
  - Moire.BLShiftRange: last x shift of the top layer.
* - `Moire.BilayerTopYMin`
  - Moire.BLShiftRange: first y shift of the top layer.
* - `Moire.BilayerTopYMax`
  - Moire.BLShiftRange: last y shift of the top layer.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Section <api-sys-syssection>` (= `SysSection`)

## Called by

{ref}`CalcSelect <api-calc-calcselect>`
