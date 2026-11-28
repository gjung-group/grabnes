---
edit_url: null
---

(api-kubo-kubotevol)=
# `KuboTEvol`

```{index} single: KuboTEvol (subroutine in kubo)
```

```fortran
subroutine KuboTEvol(Psi, ZUPsi, Psin, Psinm1, XpnPsi, XpnPsim1, c, dT, nT, az, bz, H, H0, hopp, NList, NeighD, ac, bc, nEn, Emin, Emax, eps, nRecurs, nPol, nWr)
```

- **Module:** {ref}`kubo <api-kubo>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubo.F90, line 321`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubo.F90#L321)

## Parameters

`Psi` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`ZUPsi` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`Psin` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`Psinm1` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`XpnPsi` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`XpnPsim1` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`c` : `complex(dp), dimension(nPol), intent(in)`
: —

`dT` : `real(dp), intent(in)`
: —

`nT` : `integer, intent(in)`
: —

`az` : `real(dp), dimension(nRecurs), intent(inout)`
: —

`bz` : `real(dp), dimension(nRecurs), intent(inout)`
: —

`H` : `complex(dp), dimension(inode1:), intent(out)`
: —

`H0` : `real(dp), dimension(inode1:), intent(in)`
: —

`hopp` : `complex(dp), dimension(:,inode1:), intent(in)`
: —

`NList` : `integer, dimension(1:,inode1:), intent(in)`
: —

`NeighD` : `real(dp), dimension(1:,1:,inode1:), intent(in)`
: —

`ac` : `real(dp), intent(in)`
: —

`bc` : `real(dp), intent(in)`
: —

`nEn` : `integer, intent(in)`
: —

`Emin` : `real(dp), intent(in)`
: —

`Emax` : `real(dp), intent(in)`
: —

`eps` : `real(dp), intent(in)`
: —

`nRecurs` : `integer, intent(in)`
: —

`nPol` : `integer, intent(in)`
: —

`nWr` : `integer, intent(in)`
: —

## Calls

{ref}`KuboEvol <api-kubosubs-kuboevol>`, {ref}`KuboFrac <api-kubosubs-kubofrac>`, {ref}`KuboRecursion <api-kubosubs-kuborecursion>`, {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`)

## Called by

[CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`)
