---
edit_url: null
---

(api-scf)=
# `scf`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/scf.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/scf.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`atoms <api-atoms>`, {ref}`cell <api-cell>`, {ref}`constants <api-constants>`, {ref}`ham <api-ham>`, {ref}`mio <api-mio>`, {ref}`neigh <api-neigh>`, {ref}`tbpar <api-tbpar>`
- **Used by:** {ref}`calc <api-calc>`, {ref}`diag <api-diag>`
- **Public interface:** 2 procedures, 2 variables
- **Input keys:** 8 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`SCFGetCharge <api-scf-scfgetcharge>`
  - subroutine
  - —
* - {ref}`SCFInit <api-scf-scfinit>`
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
* - `charge`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `Zch`
  - `real(dp), parameter`
  - `1.0_dp`
  - —
```
