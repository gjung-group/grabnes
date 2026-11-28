---
edit_url: null
---

(api-neigh-fastnnnotsquarenotrectangle)=
# `fastNNnotsquareNotRectangle`

```{index} single: fastNNnotsquareNotRectangle (subroutine in neigh)
```

```fortran
subroutine fastNNnotsquareNotRectangle(natoms, x, y, z, aCC, cutoff2, cutoff2bis, A1, A2, maxnn)
```

- **Module:** {ref}`neigh <api-neigh>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/neigh.F90, line 1983`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L1983)

## Parameters

`natoms` : `integer, intent(in)`
: —

`x` : `real(dp), dimension(natoms), intent(in)`
: —

`y` : `real(dp), dimension(natoms), intent(in)`
: —

`z` : `real(dp), dimension(natoms), intent(in)`
: —

`aCC` : `real(dp)`
: —

`cutoff2` : `real(dp)`
: —

`cutoff2bis` : `real(dp)`
: —

`A1` : `real(dp), dimension(3)`
: —

`A2` : `real(dp), dimension(3)`
: —

`maxnn` : `integer, intent(in)`
: —

## Calls

{ref}`AtomsSetCart <api-atoms-atomssetcart>`, {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MPIBarrier <api-mpi-mpibarrier>`, [NeighCheckCount](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L42) ({ref}`neigh <api-neigh>`), [NeighCountsPrint](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L2708) ({ref}`neigh <api-neigh>`), {ref}`MIO_StringComp <api-string-interface-stringcomp>` (= `StringComp`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

{ref}`HamInit <api-ham-haminit>`
