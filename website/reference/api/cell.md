---
edit_url: null
---

(api-cell)=
# `cell`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/cell.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/cell.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`constants <api-constants>`, {ref}`math <api-math>`, {ref}`mio <api-mio>`, {ref}`name <api-name>`
- **Used by:** {ref}`atoms <api-atoms>`, {ref}`calc <api-calc>`, {ref}`diag <api-diag>`, {ref}`gauss <api-gauss>`, {ref}`grabnes <api-grabnes>`, {ref}`ham <api-ham>`, {ref}`hybrid <api-hybrid>`, {ref}`interface <api-interface>`, {ref}`magf <api-magf>`, {ref}`moireBLShift <api-moireblshift>`, {ref}`neigh <api-neigh>`, {ref}`scf <api-scf>`, {ref}`tbpar <api-tbpar>`
- **Public interface:** 1 procedure, 8 variables
- **Input keys:** 14 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`CellGet <api-cell-cellget>`
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
* - `aBN`
  - `real(dp), save`
  - —
  - —
* - `aG`
  - `real(dp), save`
  - —
  - —
* - `area`
  - `real(dp), save`
  - —
  - —
* - `rcell`
  - `real(dp), dimension(3,3), save`
  - —
  - —
* - `sCell`
  - `integer, save`
  - —
  - —
* - `sCell2`
  - `integer, save`
  - —
  - —
* - `ucell`
  - `real(dp), dimension(3,3), save`
  - —
  - —
* - `volume`
  - `real(dp), save`
  - —
  - —
```
