---
edit_url: null
---

(api-mpi-mpibarrier)=
# `MPIBarrier`

```{index} single: MPIBarrier (subroutine in mpi)
```

```fortran
subroutine MPIBarrier()
```

- **Module:** {ref}`mpi <api-mpi>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/MPI/mpi.f90, line 108`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/MPI/mpi.f90#L108)

## Calls

{ref}`MPITimer <api-mpitime-mpitimer>`

## Called by

{ref}`fastNN <api-neigh-fastnn>`, {ref}`fastNNnotsquareNotRectangle <api-neigh-fastnnnotsquarenotrectangle>`, [HamUpdateNormals](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4578) ({ref}`ham <api-ham>`), {ref}`KuboCn <api-kubosubs-kubocn>`, {ref}`KuboFrac <api-kubosubs-kubofrac>`, {ref}`KuboRecursion <api-kubosubs-kuborecursion>`, [KuboUpdate\_d](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L621) ({ref}`kubosubs <api-kubosubs>`), [KuboUpdate\_z](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L648) ({ref}`kubosubs <api-kubosubs>`), {ref}`NeighList <api-neigh-neighlist>`, [NeighListOld](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L1023) ({ref}`neigh <api-neigh>`)
