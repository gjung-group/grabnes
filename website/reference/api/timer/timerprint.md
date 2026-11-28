---
edit_url: null
---

(api-timer-timerprint)=
# `TimerPrint`

```{index} single: TimerPrint (subroutine in timer)
```

```fortran
subroutine TimerPrint()
```

- **Module:** {ref}`timer <api-timer>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/timer.F90, line 167`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/timer.F90#L167)

Print timers

## Calls

{ref}`MPITimerGet <api-mpitime-mpitimerget>`, {ref}`SysPrint <api-sys-interface-sysprint>`, {ref}`SysWarning <api-sys-syswarning>`, {ref}`TimerStop <api-timer-timerstop>`

## Called by

{ref}`MIO_Finalize <api-mio-mio_finalize>`
