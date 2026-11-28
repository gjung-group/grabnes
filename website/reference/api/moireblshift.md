---
edit_url: null
---

(api-moireblshift)=
# `moireBLShift`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/moireBLShift.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/moireBLShift.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`cell <api-cell>`, {ref}`mio <api-mio>`, {ref}`name <api-name>`
- **Used by:** {ref}`calc <api-calc>`, {ref}`ham <api-ham>`
- **Public interface:** 2 procedures, 7 variables
- **Input keys:** 15 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`moireBLShiftInit <api-moireblshift-moireblshiftinit>`
  - subroutine
  - —
* - {ref}`moireBLShiftValue <api-moireblshift-moireblshiftvalue>`
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
* - `mSf`
  - `integer, save`
  - `0`
  - —
* - `mSi`
  - `integer, save`
  - `0`
  - —
* - `mSStep`
  - `integer, save`
  - `1`
  - —
* - `tauX1`
  - `real(dp), save`
  - —
  - —
* - `tauX2`
  - `real(dp), save`
  - —
  - —
* - `tauY1`
  - `real(dp), save`
  - —
  - —
* - `tauY2`
  - `real(dp), save`
  - —
  - —
```
