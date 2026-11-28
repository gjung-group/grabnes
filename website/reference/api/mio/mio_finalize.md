---
edit_url: null
---

(api-mio-mio_finalize)=
# `MIO_Finalize`

```{index} single: MIO_Finalize (subroutine in mio)
```

```fortran
subroutine MIO_Finalize()
```

- **Module:** {ref}`mio <api-mio>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/mio.F90, line 184`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mio.F90#L184)

Finalize MIO library

## Calls

{ref}`CloseAll <api-file-closeall>`, {ref}`DebugClose <api-debug-debugclose>`, {ref}`InputClose <api-input-inputclose>`, {ref}`MemPrint <api-mem-memprint>`, {ref}`MPIFin <api-mpi-mpifin>`, {ref}`num2month <api-format-num2month>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`SysPrintErrors <api-sys-sysprinterrors>`, {ref}`MIO_Section <api-sys-syssection>` (= `SysSection`), {ref}`TimerPrint <api-timer-timerprint>`

## Called by

{ref}`GrabnesEnd <api-grabnes_end-grabnesend>`
