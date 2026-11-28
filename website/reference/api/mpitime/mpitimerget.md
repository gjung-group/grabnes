---
edit_url: null
---

(api-mpitime-mpitimerget)=
# `MPITimerGet`

```{index} single: MPITimerGet (subroutine in mpitime)
```

```fortran
subroutine MPITimerGet(t, nc)
```

- **Module:** {ref}`mpitime <api-mpitime>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/MPI/mpitime.F90, line 78`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/MPI/mpitime.F90#L78)

Get the time

## Parameters

`t` : `real(dp), intent(out)`
: —

`nc` : `integer, intent(out)`
: —

## Called by

{ref}`TimerPrint <api-timer-timerprint>`
