---
edit_url: null
---

(api-neigh-fastnn)=
# `fastNN`

```{index} single: fastNN (subroutine in neigh)
```

```fortran
subroutine fastNN(natoms, x, y, z, aCC, cutoff, A1, A2, maxnn)
```

- **Module:** {ref}`neigh <api-neigh>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/neigh.F90, line 2414`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L2414)

calculate nearest neighbors of a sample with a fast method that uses binning

## Parameters

`natoms` : `integer`
: —

`x` : `double precision, dimension(natoms)`
: —

`y` : `double precision, dimension(natoms)`
: —

`z` : `double precision, dimension(natoms)`
: —

`aCC` : `double precision`
: —

`cutoff` : `double precision`
: —

`A1` : `real(dp), dimension(3), intent(in)`
: —

`A2` : `real(dp), dimension(3), intent(in)`
: —

`maxnn` : `integer`
: —

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MPIBarrier <api-mpi-mpibarrier>`, [NeighCheckCount](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L42) ({ref}`neigh <api-neigh>`), [NeighCountsPrint](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L2708) ({ref}`neigh <api-neigh>`), {ref}`MIO_StringComp <api-string-interface-stringcomp>` (= `StringComp`)

## Called by

{ref}`HamInit <api-ham-haminit>`
