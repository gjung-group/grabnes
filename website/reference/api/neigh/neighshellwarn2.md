---
edit_url: null
---

(api-neigh-neighshellwarn2)=
# `NeighShellWarn2`

```{index} single: NeighShellWarn2 (function in neigh)
```

```fortran
function NeighShellWarn2(nshell, aCC) result(warn2)
```

- **Module:** {ref}`neigh <api-neigh>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/neigh.F90, line 86`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L86)

Squared in-plane length beyond which an intralayer pair is closer to the search radius of NeighShellCutoff2 than a rigid lattice allows: the last quarter of the interval between the last requested shell and the radius. Pairs found there belong to a strained structure.

## Parameters

`nshell` : `integer, intent(in)`
: —

`aCC` : `real(dp), intent(in)`
: —

## Returns

`warn2` : `real(dp)`
: —

## Called by

{ref}`HamInit <api-ham-haminit>`
