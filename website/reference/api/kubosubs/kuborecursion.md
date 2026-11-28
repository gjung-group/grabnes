---
edit_url: null
---

(api-kubosubs-kuborecursion)=
# `KuboRecursion`

```{index} single: KuboRecursion (subroutine in kubosubs)
```

```fortran
subroutine KuboRecursion(Psi, Psin, Psinm1, nRecurs, a, b, H, H0, hopp, NList, wr, name, store, overlapWith, overlaps)
```

- **Module:** {ref}`kubosubs <api-kubosubs>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubosubs.F90, line 26`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L26)

Lanczos recursion from the normalised state Psi. The optional arguments do not change the recursion: store(:,n) receives the n-th Lanczos vector (store(:,1) = Psi), and overlaps(n) = \<q\_n|overlapWith>.

## Parameters

`Psi` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`Psin` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`Psinm1` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`nRecurs` : `integer, intent(in)`
: —

`a` : `real(dp), dimension(nRecurs), intent(out)`
: —

`b` : `real(dp), dimension(nRecurs), intent(out)`
: —

`H` : `complex(dp), dimension(inode1:), intent(out)`
: —

`H0` : `real(dp), dimension(inode1:), intent(in)`
: —

`hopp` : `complex(dp), dimension(:,inode1:), intent(in)`
: —

`NList` : `integer, dimension(1:,inode1:), intent(in)`
: —

`wr` : `logical, intent(in)`
: —

`name` : `character(len=*), intent(in), optional`
: —

`store` : `complex(dp), dimension(:,:), intent(out), optional`
: —

`overlapWith` : `complex(dp), dimension(inode1:), intent(in), optional`
: —

`overlaps` : `complex(dp), dimension(nRecurs), intent(out), optional`
: —

## Calls

[GetUnit](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L606) ({ref}`file <api-file>`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`MPIBarrier <api-mpi-mpibarrier>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`)

## Called by

{ref}`KuboDOS <api-kubo-kubodos>`, {ref}`KuboTEvol <api-kubo-kubotevol>`
