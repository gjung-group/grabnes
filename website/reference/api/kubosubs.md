---
edit_url: null
---

(api-kubosubs)=
# `kubosubs`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/kubosubs.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`atoms <api-atoms>`, {ref}`constants <api-constants>`, {ref}`kuboarrays <api-kuboarrays>`, {ref}`math <api-math>`, {ref}`mio <api-mio>`, {ref}`name <api-name>`, {ref}`neigh <api-neigh>`
- **Used by:** {ref}`calc <api-calc>`, {ref}`kubo <api-kubo>`
- **Public interface:** 7 procedures
- **Input keys:** 1 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`KuboChebyshevEvolve <api-kubosubs-kubochebyshevevolve>`
  - subroutine
  - One Chebyshev time step of a single state: Z \<- U(dT) Z, with the coefficients c of KuboCn.
* - {ref}`KuboCn <api-kubosubs-kubocn>`
  - subroutine
  - —
* - {ref}`KuboEvol <api-kubosubs-kuboevol>`
  - subroutine
  - One time step of Psi and of the commutator ZUPsi = \[U(t), X\] Psi\_0 (Roche-Mayou): Psi \<- U(dT) Psi, ZUPsi \<- U(dT) ZUPsi + \[U(dT), X\] Psi.
* - {ref}`KuboFermi <api-kubosubs-kubofermi>`
  - subroutine
  - —
* - {ref}`KuboFrac <api-kubosubs-kubofrac>`
  - subroutine
  - —
* - {ref}`KuboInterval <api-kubosubs-kubointerval>`
  - subroutine
  - —
* - {ref}`KuboRecursion <api-kubosubs-kuborecursion>`
  - subroutine
  - Lanczos recursion from the normalised state Psi.
```
