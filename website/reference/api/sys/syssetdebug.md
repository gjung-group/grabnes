---
edit_url: null
---

(api-sys-syssetdebug)=
# `SysSetDebug`

```{index} single: SysSetDebug (subroutine in sys)
```

```fortran
subroutine SysSetDebug(on, module)
```

- **Module:** {ref}`sys <api-sys>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/sys.F90, line 411`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/sys.F90#L411)

Switch debug output on or off, for every module (no module given) or for one module.

## Parameters

`on` : `logical, intent(in)`
: —

`module` : `character(len=*), intent(in), optional`
: —

## Called by

{ref}`MIO_Initialize <api-mio-mio_initialize>`
