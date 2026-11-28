---
edit_url: null
---

(api-interface)=
# `interface`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/interface.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/interface.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`atoms <api-atoms>`, {ref}`cell <api-cell>`, {ref}`constants <api-constants>`, {ref}`math <api-math>`, {ref}`mio <api-mio>`, {ref}`neigh <api-neigh>`, {ref}`tbpar <api-tbpar>`
- **Used by:** {ref}`diag <api-diag>`, {ref}`ham <api-ham>`
- **Public interface:** 2 procedures, 7 variables
- **Input keys:** 13 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`InterfacePot1 <api-interface-interfacepot1>`
  - subroutine
  - —
* - {ref}`InterfacePot2 <api-interface-interfacepot2>`
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
* - `edgeH`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `edgeHopp`
  - `logical, save`
  - —
  - —
* - `edgeIndx`
  - `integer, dimension(:), pointer`
  - —
  - —
* - `NedgeCell`
  - `integer, dimension(:,:,:), pointer`
  - —
  - —
* - `nEdgeN`
  - `integer, dimension(:), pointer`
  - —
  - —
* - `NeI`
  - `integer, dimension(:,:), pointer`
  - —
  - —
* - `nQ`
  - `integer, save`
  - —
  - —
```
