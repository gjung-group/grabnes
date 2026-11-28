---
edit_url: null
---

(api-magf-haldphasevalue)=
# `HaldPhaseValue`

```{index} single: HaldPhaseValue (subroutine in magf)
```

```fortran
subroutine HaldPhaseValue(mPhi)
```

- **Module:** {ref}`magf <api-magf>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/magf.F90, line 213`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/magf.F90#L213)

## Parameters

`mPhi` : `integer, intent(in)`
: —

## Calls

{ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Section <api-sys-syssection>` (= `SysSection`)

## Called by

{ref}`CalcSelect <api-calc-calcselect>`
