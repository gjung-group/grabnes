---
edit_url: null
---

(api-moireblshift-moireblshiftinit)=
# `moireBLShiftInit`

```{index} single: moireBLShiftInit (subroutine in moireBLShift)
```

```fortran
subroutine moireBLShiftInit()
```

- **Module:** {ref}`moireBLShift <api-moireblshift>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/moireBLShift.F90, line 17`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/moireBLShift.F90#L17)

## Input keys read (3)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Moire.BLShiftDefinedValues`
  - moireBLShift.Calc: one calculation with the shifts Moire.Bilayer{Bottom,Top}{X,Y}Shift.
* - `Moire.BLShiftRange`
  - moireBLShift.Calc: sweep each shift from its Moire.Bilayer\*Min to \*Max in Moire.BLShiftRangeSteps steps.
* - `Moire.BLShiftRangeSteps`
  - Moire.BLShiftRange: number of steps of the shift sweep.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`

## Called by

{ref}`CalcSelect <api-calc-calcselect>`
