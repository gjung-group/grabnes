---
edit_url: null
---

(api-kubosubs-kubofermi)=
# `KuboFermi`

```{index} single: KuboFermi (subroutine in kubosubs)
```

```fortran
subroutine KuboFermi(npts)
```

- **Module:** {ref}`kubosubs <api-kubosubs>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubosubs.F90, line 677`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L677)

## Parameters

`npts` : `integer, intent(in)`
: —

## Calls

[GetUnit](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L606) ({ref}`file <api-file>`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)
