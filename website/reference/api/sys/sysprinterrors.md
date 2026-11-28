---
edit_url: null
---

(api-sys-sysprinterrors)=
# `SysPrintErrors`

```{index} single: SysPrintErrors (subroutine in sys)
```

```fortran
subroutine SysPrintErrors(final)
```

- **Module:** {ref}`sys <api-sys>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/sys.F90, line 284`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/sys.F90#L284)

Print the count of errors and warnings.

## Parameters

`final` : `logical, intent(in), optional`
: If is not the final count a message is displayed ('Current count')

## Calls

{ref}`num2string <api-format-interface-num2string>`, {ref}`SysPrint <api-sys-interface-sysprint>`

## Called by

{ref}`MIO_Finalize <api-mio-mio_finalize>`, {ref}`SysKill <api-sys-syskill>`, {ref}`SysKillFlush <api-sys-syskillflush>`
