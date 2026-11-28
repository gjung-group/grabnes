---
edit_url: null
---

(api-neigh-neighshellcutoff2)=
# `NeighShellCutoff2`

```{index} single: NeighShellCutoff2 (function in neigh)
```

```fortran
function NeighShellCutoff2(nshell, aCC) result(cut2)
```

- **Module:** {ref}`neigh <api-neigh>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/neigh.F90, line 63`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L63)

Squared in-plane search radius that contains the first nshell neighbour shells of a honeycomb lattice with bond length aCC.

The radius lies halfway between shell nshell and the next one, so it tolerates moderately strained or relaxed lattices.

## Parameters

`nshell` : `integer, intent(in)`
: —

`aCC` : `real(dp), intent(in)`
: —

## Returns

`cut2` : `real(dp)`
: —

## Calls

{ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`)

## Called by

{ref}`HamInit <api-ham-haminit>`
