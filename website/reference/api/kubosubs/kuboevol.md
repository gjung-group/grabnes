---
edit_url: null
---

(api-kubosubs-kuboevol)=
# `KuboEvol`

```{index} single: KuboEvol (subroutine in kubosubs)
```

```fortran
subroutine KuboEvol(Psi, ZUPsi, Psin, Psinm1, XpnPsi, XpnPsim1, c, H0, hopp, NList, NeighD, ac, bc, nPol, prnt, dir)
```

- **Module:** {ref}`kubosubs <api-kubosubs>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubosubs.F90, line 438`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L438)

One time step of Psi and of the commutator ZUPsi = \[U(t), X\] Psi\_0 (Roche-Mayou): Psi \<- U(dT) Psi, ZUPsi \<- U(dT) ZUPsi + \[U(dT), X\] Psi. The optional dir (1 = x, 2 = y) selects the direction; without it the direction comes from Kubo.TimeEvolutionInYDirection, as before.

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

`nPol` : `integer, intent(in)`
: —

`prnt` : `logical, intent(in)`
: —

`dir` : `integer, intent(in), optional`
: —

## Input keys read (1)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Kubo.TimeEvolutionInYDirection`
  - Kubo time evolution: position operator along y instead of x.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`KuboChebyshevEvolve <api-kubosubs-kubochebyshevevolve>`, {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)

## Called by

{ref}`KuboTEvol <api-kubo-kubotevol>`
