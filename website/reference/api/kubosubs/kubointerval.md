---
edit_url: null
---

(api-kubosubs-kubointerval)=
# `KuboInterval`

```{index} single: KuboInterval (subroutine in kubosubs)
```

```fortran
subroutine KuboInterval(nRecurs, a, b, ac, bc)
```

- **Module:** {ref}`kubosubs <api-kubosubs>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubosubs.F90, line 290`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L290)

## Parameters

`nRecurs` : `integer, intent(in)`
: —

`a` : `real(dp), dimension(nRecurs), intent(in)`
: —

`b` : `real(dp), dimension(nRecurs), intent(in)`
: —

`ac` : `real(dp), intent(out)`
: —

`bc` : `real(dp), intent(out)`
: —

## Calls

{ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)

## Called by

[CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`)
