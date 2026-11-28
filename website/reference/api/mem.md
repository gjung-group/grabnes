---
edit_url: null
---

(api-mem)=
# `mem`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/MIO/mem_in.f90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mem_in.f90#L13)
- **Group:** {ref}`MIO library <api-group-mio>`
- **Uses:** {ref}`format <api-format>`, {ref}`precision <api-precision>`, {ref}`sys <api-sys>`
- **Used by:** {ref}`input <api-input>`, {ref}`mio <api-mio>`, {ref}`timer <api-timer>`
- **Public interface:** 3 procedures, 7 variables

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`MemCount <api-mem-memcount>`
  - subroutine
  - Counts the memory in use
* - {ref}`MemGetName <api-mem-memgetname>`
  - subroutine
  - Get names for the variable and the routine
* - {ref}`MemPrint <api-mem-memprint>`
  - subroutine
  - Print the memory usage
```

## Module variables

```{list-table}
:header-rows: 1

* - Name
  - Declaration
  - Initial value
  - Description
* - `Csz`
  - `integer, parameter`
  - `[CSZ]`
  - —
* - `Dsz`
  - `integer, parameter`
  - `[DSZ]`
  - —
* - `Isz`
  - `integer, parameter`
  - `[ISZ]`
  - —
* - `Lsz`
  - `integer, parameter`
  - `[LSZ]`
  - —
* - `Rsz`
  - `integer, parameter`
  - `[RSZ]`
  - —
* - `Ssz`
  - `integer, parameter`
  - `[SSZ]`
  - —
* - `Zsz`
  - `integer, parameter`
  - `[ZSZ]`
  - —
```
