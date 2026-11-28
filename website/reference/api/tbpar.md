---
edit_url: null
---

(api-tbpar)=
# `tbpar`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/tbpar.f90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/tbpar.f90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`cell <api-cell>`, {ref}`mio <api-mio>`
- **Used by:** {ref}`calc <api-calc>`, {ref}`diag <api-diag>`, {ref}`ham <api-ham>`, {ref}`interface <api-interface>`, {ref}`neigh <api-neigh>`, {ref}`scf <api-scf>`
- **Public interface:** 1 procedure, 14 variables
- **Input keys:** 20 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`TBInit <api-tbpar-tbinit>`
  - subroutine
  - —
```

## Module variables

```{list-table}
:header-rows: 1

* - Name
  - Declaration
  - Initial value
  - Description
* - `e0_B`
  - `real(dp), save`
  - —
  - —
* - `e0_C`
  - `real(dp), save`
  - —
  - —
* - `e0_C1`
  - `real(dp), save`
  - —
  - —
* - `e0_C1_LB`
  - `real(dp), save`
  - —
  - —
* - `e0_C1_LT`
  - `real(dp), save`
  - —
  - —
* - `e0_C2`
  - `real(dp), save`
  - —
  - —
* - `e0_C2_LB`
  - `real(dp), save`
  - —
  - —
* - `e0_C2_LT`
  - `real(dp), save`
  - —
  - —
* - `e0_N`
  - `real(dp), save`
  - —
  - —
* - `g0`
  - `real(dp), save`
  - —
  - —
* - `gIntLay`
  - `real(dp), save`
  - —
  - —
* - `gn`
  - `real(dp), dimension(:,:,:), pointer`
  - —
  - —
* - `tbnn`
  - `integer, save`
  - —
  - —
* - `U`
  - `real(dp), dimension(4), save`
  - —
  - —
```
