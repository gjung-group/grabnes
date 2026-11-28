---
edit_url: null
---

(api-magf)=
# `magf`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/magf.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/magf.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`cell <api-cell>`, {ref}`constants <api-constants>`, {ref}`mio <api-mio>`, {ref}`name <api-name>`
- **Used by:** {ref}`calc <api-calc>`, {ref}`diag <api-diag>`, {ref}`ham <api-ham>`
- **Public interface:** 4 procedures, 11 variables
- **Input keys:** 16 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`HaldPhaseInit <api-magf-haldphaseinit>`
  - subroutine
  - —
* - {ref}`HaldPhaseValue <api-magf-haldphasevalue>`
  - subroutine
  - —
* - {ref}`MagfInit <api-magf-magfinit>`
  - subroutine
  - —
* - {ref}`MagfValue <api-magf-magfvalue>`
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
* - `Bmag`
  - `real(dp), save`
  - `0.0_dp`
  - —
* - `BmagZeeman`
  - `real(dp), save`
  - `0.0_dp`
  - —
* - `HaldanePhase`
  - `real(dp), save`
  - —
  - —
* - `HaldPhase`
  - `logical, save`
  - `.false.`
  - —
* - `magfield`
  - `logical, save`
  - `.false.`
  - —
* - `mBf`
  - `integer, save`
  - `0`
  - —
* - `mBi`
  - `integer, save`
  - `0`
  - —
* - `mPhif`
  - `integer, save`
  - `0`
  - —
* - `mPhii`
  - `integer, save`
  - `0`
  - —
* - `mPhiStep`
  - `integer, save`
  - `0`
  - —
* - `mStep`
  - `integer, save`
  - `1`
  - —
```
