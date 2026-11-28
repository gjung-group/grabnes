---
edit_url: null
---

(api-mio-mio_initialize)=
# `MIO_Initialize`

```{index} single: MIO_Initialize (subroutine in mio)
```

```fortran
subroutine MIO_Initialize(program, success)
```

- **Module:** {ref}`mio <api-mio>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/mio.F90, line 117`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mio.F90#L117)

Initialize MIO library

## Parameters

`program` : `character(len=*), intent(in)`
: Name of the program

`success` : `logical, intent(out)`
: —

## Calls

{ref}`DebugInit <api-debug-debuginit>`, {ref}`InputInit <api-input-inputinit>`, {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MPIIni <api-mpi-mpiini>`, {ref}`num2month <api-format-num2month>`, {ref}`OmpInit <api-omp-ompinit>`, {ref}`OptionsProcess <api-options-optionsprocess>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_SetDebug <api-sys-syssetdebug>` (= `SysSetDebug`), {ref}`TimerCount <api-timer-timercount>`

## Called by

{ref}`GrabnesStart <api-grabnes_start-grabnesstart>`
