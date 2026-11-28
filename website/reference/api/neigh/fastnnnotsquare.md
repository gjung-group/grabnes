---
edit_url: null
---

(api-neigh-fastnnnotsquare)=
# `fastNNnotsquare`

```{index} single: fastNNnotsquare (subroutine in neigh)
```

```fortran
subroutine fastNNnotsquare(natoms, x, y, z, cutoff2, cutoff2bis, A1, A2)
```

- **Module:** {ref}`neigh <api-neigh>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/neigh.F90, line 1433`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L1433)

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

## Input keys read (2)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Output.ReadDataFiles`
  - Here: the neighbour lists are read from the files v, dx, dy, dz instead of searched.
* - `Interlayer.AddSecondLayerInteractions`
  - Also pair atoms two layers apart (4.5 \< \|dz\| \< 7.5 Angstrom) in the neighbour search, for the models with next-layer hoppings.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, [NeighFinishSearch](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L319) ({ref}`neigh <api-neigh>`), [NeighReadDataFiles](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L446) ({ref}`neigh <api-neigh>`), [NeighSearchLayered](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L117) ({ref}`neigh <api-neigh>`)

## Called by

{ref}`HamInit <api-ham-haminit>`
