---
edit_url: null
---

(api-grabnes_start-grabnesstart)=
# `GrabnesStart`

```{index} single: GrabnesStart (subroutine in grabnes_start)
```

```fortran
subroutine GrabnesStart(success)
```

- **Module:** {ref}`grabnes_start <api-grabnes_start>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/grabnes_start.F90, line 11`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/grabnes_start.F90#L11)

## Parameters

`success` : `logical, intent(out)`
: —

## Input keys read (2)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Run.Description`
  - Free text printed at the top of the output.
* - `Run.Prefix`
  - Prefix of the output files (\<prefix>.DOS, ...); default grabnes.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Initialize <api-mio-mio_initialize>`, [StartHeader](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/grabnes_start.F90#L53) ({ref}`grabnes_start <api-grabnes_start>`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)

## Called by

{ref}`grabnes <api-grabnes>` (program)
