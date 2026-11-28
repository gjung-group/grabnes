---
edit_url: null
---

(api-kuboarrays)=
# `kuboarrays`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/kuboarrays.f90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kuboarrays.f90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`mio <api-mio>`
- **Used by:** {ref}`calc <api-calc>`, {ref}`kubosubs <api-kubosubs>`
- **Public interface:** 12 variables

## Module variables

```{list-table}
:header-rows: 1

* - Name
  - Declaration
  - Initial value
  - Description
* - `a`
  - `real(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `b`
  - `real(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `c`
  - `complex(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `H`
  - `complex(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `Psi`
  - `complex(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `Psin`
  - `complex(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `Psinm1`
  - `complex(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `tempD`
  - `real(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `tempZ`
  - `complex(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `XpnPsi`
  - `complex(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `XpnPsim1`
  - `complex(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `ZUPsi`
  - `complex(dp), dimension(:), pointer`
  - `NULL()`
  - —
```
