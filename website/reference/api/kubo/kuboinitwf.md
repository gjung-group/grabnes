---
edit_url: null
---

(api-kubo-kuboinitwf)=
# `KuboInitWF`

```{index} single: KuboInitWF (subroutine in kubo)
```

```fortran
subroutine KuboInitWF(Psi)
```

- **Module:** {ref}`kubo <api-kubo>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubo.F90, line 24`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubo.F90#L24)

## Parameters

`Psi` : `complex(dp), dimension(inode1:), intent(out)`
: —

## Input keys read (2)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Kubo.SetSeed`
  - Seed the random-phase start vector with Kubo.SeedValue; otherwise from the clock (runs not reproducible).
* - `Kubo.SeedValue`
  - Kubo.SetSeed: seed of the random-phase start vector; use different values to average over vectors.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)

## Called by

[CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`)
