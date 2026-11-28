---
edit_url: null
---

(api-parallel)=
# `parallel`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/parallel.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/parallel.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`mio <api-mio>`
- **Used by:** {ref}`atoms <api-atoms>`, {ref}`gauss <api-gauss>`, {ref}`grabnes <api-grabnes>`, {ref}`ham <api-ham>`, {ref}`kubo <api-kubo>`, {ref}`neigh <api-neigh>`, {ref}`random <api-random>`
- **Public interface:** 1 procedure, 7 variables

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`ParallelDiv <api-parallel-paralleldiv>`
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
* - `nDiv`
  - `integer, save`
  - —
  - —
* - `nTh`
  - `integer, dimension(:), pointer`
  - —
  - —
* - `procID`
  - `integer, save`
  - —
  - —
* - `xDiv`
  - `integer, save`
  - —
  - —
* - `xSubdiv`
  - `integer, save`
  - —
  - —
* - `yDiv`
  - `integer, save`
  - —
  - —
* - `ySubdiv`
  - `integer, save`
  - —
  - —
```
