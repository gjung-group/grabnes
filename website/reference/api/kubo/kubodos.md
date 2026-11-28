---
edit_url: null
---

(api-kubo-kubodos)=
# `KuboDOS`

```{index} single: KuboDOS (subroutine in kubo)
```

```fortran
subroutine KuboDOS(Psi, Psin, Psinm1, a, b, H, H0, hopp, NList, nRecurs, eps, nEn, Emin, Emax)
```

- **Module:** {ref}`kubo <api-kubo>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubo.F90, line 276`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubo.F90#L276)

## Parameters

`Psi` : `complex(dp), dimension(inode1:), intent(out)`
: —

`Psin` : `complex(dp), dimension(inode1:), intent(out)`
: —

`Psinm1` : `complex(dp), dimension(inode1:), intent(out)`
: —

`a` : `real(dp), dimension(nRecurs), intent(out)`
: —

`b` : `real(dp), dimension(nRecurs), intent(out)`
: —

`H` : `complex(dp), dimension(inode1:), intent(out)`
: —

`H0` : `real(dp), dimension(inode1:), intent(in)`
: —

`hopp` : `complex(dp), dimension(:,inode1:), intent(in)`
: —

`NList` : `integer, dimension(1:,inode1:), intent(in)`
: —

`nRecurs` : `integer, intent(in)`
: —

`eps` : `real(dp), intent(in)`
: —

`nEn` : `integer, intent(in)`
: —

`Emin` : `real(dp)`
: —

`Emax` : `real(dp)`
: —

## Calls

{ref}`KuboFrac <api-kubosubs-kubofrac>`, {ref}`KuboRecursion <api-kubosubs-kuborecursion>`, {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)

## Called by

[CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`)
