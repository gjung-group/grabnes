---
edit_url: null
---

(api-scf-scfinit)=
# `SCFInit`

```{index} single: SCFInit (subroutine in scf)
```

```fortran
subroutine SCFInit(N)
```

- **Module:** {ref}`scf <api-scf>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/scf.F90, line 259`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/scf.F90#L259)

## Parameters

`N` : `integer, intent(in)`
: —

## Calls

{ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`)

## Called by

[CalcDiag](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L678) ({ref}`calc <api-calc>`), {ref}`SCFGetCharge <api-scf-scfgetcharge>`
