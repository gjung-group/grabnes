---
edit_url: null
---

(api-random)=
# `random` (`random.F90`)

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/random.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/random.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`mio <api-mio>`, {ref}`name <api-name>`
- **Used by:** {ref}`gauss <api-gauss>`, {ref}`ham <api-ham>`
- **Public interface:** 2 procedures, 1 derived type
- **Input keys:** 2 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`RandSeedFromInput <api-random-randseedfrominput>`
  - subroutine
  - Seed the intrinsic random-number generator for the disorder models.
* - {ref}`RandTest <api-random-randtest>`
  - subroutine
  - —
```

## Derived types

(api-random-type-rand_t)=
### `rand_t`

```{list-table}
:header-rows: 1

* - Component
  - Declaration
  - Default
  - Description
* - `state`
  - `integer, dimension(ns)`
  - `default_seed`
  - —
```

Defined in [`lanczosKuboCode/Src/random.F90, line 21`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/random.F90#L21).
