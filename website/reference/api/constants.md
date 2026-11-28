---
edit_url: null
---

(api-constants)=
# `constants`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/constants.f90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/constants.f90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`mio <api-mio>`
- **Used by:** {ref}`atoms <api-atoms>`, {ref}`calc <api-calc>`, {ref}`cell <api-cell>`, {ref}`diag <api-diag>`, {ref}`ham <api-ham>`, {ref}`hybrid <api-hybrid>`, {ref}`interface <api-interface>`, {ref}`kubo <api-kubo>`, {ref}`kubosubs <api-kubosubs>`, {ref}`magf <api-magf>`, {ref}`neigh <api-neigh>`, {ref}`scf <api-scf>`
- **Public interface:** 10 variables

## Module variables

```{list-table}
:header-rows: 1

* - Name
  - Declaration
  - Initial value
  - Description
* - `cmplx_0`
  - `complex(dp), parameter`
  - `(0.0_dp, 0.0_dp)`
  - —
* - `cmplx_1`
  - `complex(dp), parameter`
  - `(1.0_dp, 0.0_dp)`
  - —
* - `cmplx_i`
  - `complex(dp), parameter`
  - `(0.0_dp, 1.0_dp)`
  - —
* - `fluxq`
  - `real(dp), parameter`
  - `hPlanck/qe`
  - —
* - `hbar`
  - `real(dp), parameter`
  - `hPlanck/twopi`
  - —
* - `hPlanck`
  - `real(dp), parameter`
  - `6.62606896d-34`
  - —
* - `pi`
  - `real(dp), parameter`
  - `3.14159265358979323846264338327950288419716939937510582097494459230781640_dp`
  - —
* - `qe`
  - `real(dp), parameter`
  - `1.602176565d-19`
  - —
* - `twopi`
  - `real(dp), parameter`
  - `2.0_dp*pi`
  - —
* - `uB`
  - `real(dp), parameter`
  - `5.7883818066d-5`
  - —
```
