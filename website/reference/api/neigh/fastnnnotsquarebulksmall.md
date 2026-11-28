---
edit_url: null
---

(api-neigh-fastnnnotsquarebulksmall)=
# `fastNNnotsquareBulkSmall`

```{index} single: fastNNnotsquareBulkSmall (subroutine in neigh)
```

```fortran
subroutine fastNNnotsquareBulkSmall(natoms, x, y, z, cutoff2, cutoff2bis, A1, A2, A3)
```

- **Module:** {ref}`neigh <api-neigh>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/neigh.F90, line 1942`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L1942)

## Parameters

`natoms` : `integer, intent(in)`
: —

`x` : `real(dp), dimension(natoms), intent(in)`
: —

`y` : `real(dp), dimension(natoms), intent(in)`
: —

`z` : `real(dp), dimension(natoms), intent(in)`
: —

`cutoff2` : `real(dp), intent(in)`
: —

`cutoff2bis` : `real(dp), intent(in)`
: —

`A1` : `real(dp), dimension(3), intent(in)`
: —

`A2` : `real(dp), dimension(3), intent(in)`
: —

`A3` : `real(dp), intent(in)`
: —

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, [NeighFinishSearch](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L319) ({ref}`neigh <api-neigh>`), [NeighReadDataFiles](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L446) ({ref}`neigh <api-neigh>`), [NeighSearchLayered](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L117) ({ref}`neigh <api-neigh>`)

## Called by

{ref}`HamInit <api-ham-haminit>`
