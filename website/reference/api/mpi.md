---
edit_url: null
---

(api-mpi)=
# `mpi`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/MIO/MPI/mpi.f90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/MPI/mpi.f90#L13)
- **Group:** {ref}`MIO library <api-group-mio>`
- **Uses:** `mpi_c`, `mpi_d`, `mpi_i`, {ref}`mpi_inc <api-mpi_inc>`, `mpi_l`, `mpi_r`, `mpi_s`, `mpi_z`, {ref}`mpitime <api-mpitime>`
- **Used by:** {ref}`debug <api-debug>`, {ref}`file <api-file>`, {ref}`input <api-input>`, {ref}`mio <api-mio>`, {ref}`sys <api-sys>`, {ref}`timer <api-timer>`
- **Public interface:** 5 procedures

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`MPIBarrier <api-mpi-mpibarrier>`
  - subroutine
  - —
* - {ref}`MPIFin <api-mpi-mpifin>`
  - subroutine
  - Terminate MPI
* - {ref}`MPIGetNodes <api-mpi-mpigetnodes>`
  - subroutine
  - Return the rank of the process and the total number of processes
* - {ref}`MPIIni <api-mpi-mpiini>`
  - subroutine
  - Initialize MPI
* - {ref}`MPIIsPNode <api-mpi-mpiispnode>`
  - function
  - True if is the principal node
```
