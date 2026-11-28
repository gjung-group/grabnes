---
edit_url: null
---

(api-neigh-fastnnnotsquaresmall)=
# `fastNNnotsquareSmall`

```{index} single: fastNNnotsquareSmall (subroutine in neigh)
```

```fortran
subroutine fastNNnotsquareSmall(natoms, x, y, z, aCC, cutoff2, cutoff2bis, A1, A2, maxnn)
```

- **Module:** {ref}`neigh <api-neigh>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/neigh.F90, line 1477`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L1477)

## Parameters

`natoms` : `integer, intent(in)`
: —

`x` : `real(dp), dimension(natoms), intent(in)`
: —

`y` : `real(dp), dimension(natoms), intent(in)`
: —

`z` : `real(dp), dimension(natoms), intent(in)`
: —

`aCC` : `real(dp), intent(in)`
: —

`cutoff2` : `real(dp), intent(in)`
: —

`cutoff2bis` : `real(dp), intent(in)`
: —

`A1` : `real(dp), dimension(3), intent(in)`
: —

`A2` : `real(dp), dimension(3), intent(in)`
: —

`maxnn` : `integer, intent(in)`
: —

## Input keys read (1)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Neigh.Only000Cell`
  - fastNNnotsquareSmall: keep only the first-neighbour image in the home cell when assigning image cells.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, [NeighCheckCount](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L42) ({ref}`neigh <api-neigh>`), [NeighCountsPrint](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L2708) ({ref}`neigh <api-neigh>`), {ref}`MIO_StringComp <api-string-interface-stringcomp>` (= `StringComp`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)
