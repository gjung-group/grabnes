---
edit_url: null
---

(api-sys-sysdebugging)=
# `SysDebugging`

```{index} single: SysDebugging (function in sys)
```

```fortran
logical function SysDebugging(module)
```

- **Module:** {ref}`sys <api-sys>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/MIO/sys.F90, line 450`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/sys.F90#L450)

True if debug output is on for this module (or for every module). Use it to skip the computation of a diagnostic that is only printed.

## Parameters

`module` : `character(len=*), intent(in)`
: —

## Returns

`SysDebugging` : `logical`
: —

## Called by

[AtomsOrder](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/atoms.F90#L1890) ({ref}`atoms <api-atoms>`), {ref}`SysPrintDebug <api-sys-sysprintdebug>`
