---
edit_url: null
---

(api-mpitime-mpitimer)=
# `MPITimer`

```{index} single: MPITimer (subroutine in mpitime)
```

```fortran
subroutine MPITimer(action)
```

- **Module:** {ref}`mpitime <api-mpitime>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/MPI/mpitime.F90, line 42`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/MPI/mpitime.F90#L42)

Count the time for MPI modules

## Parameters

`action` : `integer, intent(in)`
: if action = 1, count else, stop

## Called by

{ref}`MPIBarrier <api-mpi-mpibarrier>`, {ref}`MPIIni <api-mpi-mpiini>`
