---
edit_url: null
---

(api-kubo-kuboinitwfpdos)=
# `KuboInitWFPDOS`

```{index} single: KuboInitWFPDOS (subroutine in kubo)
```

```fortran
subroutine KuboInitWFPDOS(Psi, PDOSAtomNumber)
```

- **Module:** {ref}`kubo <api-kubo>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubo.F90, line 87`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubo.F90#L87)

## Parameters

`Psi` : `complex(dp), dimension(inode1:), intent(out)`
: —

`PDOSAtomNumber` : `integer`
: —

## Calls

{ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)

## Called by

[CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`)
