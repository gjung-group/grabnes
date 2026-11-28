---
edit_url: null
---

(api-mpi-mpigetnodes)=
# `MPIGetNodes`

```{index} single: MPIGetNodes (subroutine in mpi)
```

```fortran
subroutine MPIGetNodes(n, tn)
```

- **Module:** {ref}`mpi <api-mpi>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/MPI/mpi.f90, line 79`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/MPI/mpi.f90#L79)

Return the rank of the process and the total number of processes

## Parameters

`n` : `integer, intent(out)`
: —

`tn` : `integer, intent(out)`
: —
