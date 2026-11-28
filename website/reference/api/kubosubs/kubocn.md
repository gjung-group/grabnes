---
edit_url: null
---

(api-kubosubs-kubocn)=
# `KuboCn`

```{index} single: KuboCn (subroutine in kubosubs)
```

```fortran
subroutine KuboCn(dT, ac, bc, nPol, c)
```

- **Module:** {ref}`kubosubs <api-kubosubs>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubosubs.F90, line 352`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L352)

## Parameters

`dT` : `real(dp), intent(in)`
: —

`ac` : `real(dp), intent(in)`
: —

`bc` : `real(dp), intent(in)`
: —

`nPol` : `integer, intent(in)`
: —

`c` : `complex(dp), dimension(nPol), intent(out)`
: —

## Calls

[GetUnit](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L606) ({ref}`file <api-file>`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`MPIBarrier <api-mpi-mpibarrier>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), [TQL2](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L725) ({ref}`kubosubs <api-kubosubs>`)

## Called by

[CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`)
