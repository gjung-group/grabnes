---
edit_url: null
---

(api-kubosubs-kubofrac)=
# `KuboFrac`

```{index} single: KuboFrac (subroutine in kubosubs)
```

```fortran
subroutine KuboFrac(a, b, nRecurs, norm, eps, nEn, Emin, Emax, name, name2)
```

- **Module:** {ref}`kubosubs <api-kubosubs>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubosubs.F90, line 222`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L222)

## Parameters

`a` : `real(dp), dimension(nRecurs), intent(in)`
: —

`b` : `real(dp), dimension(nRecurs), intent(in)`
: —

`nRecurs` : `integer, intent(in)`
: —

`norm` : `real(dp), intent(in)`
: —

`eps` : `real(dp), intent(in)`
: —

`nEn` : `integer, intent(in)`
: —

`Emin` : `real(dp), intent(in)`
: —

`Emax` : `real(dp), intent(in)`
: —

`name` : `character(len=*), intent(in)`
: —

`name2` : `character(len=*), intent(in)`
: —

## Calls

[GetUnit](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L606) ({ref}`file <api-file>`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`MPIBarrier <api-mpi-mpibarrier>`

## Called by

{ref}`KuboDOS <api-kubo-kubodos>`, {ref}`KuboTEvol <api-kubo-kubotevol>`
