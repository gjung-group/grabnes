---
edit_url: null
---

(api-neigh-neighlist)=
# `NeighList`

```{index} single: NeighList (subroutine in neigh)
```

```fortran
subroutine NeighList()
```

- **Module:** {ref}`neigh <api-neigh>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/neigh.F90, line 524`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L524)

## Input keys read (9)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Structure.InterlayerDistance`
  - Here (NeighList): sets the interlayer neighbour radius sqrt((1.1 aBN)^2/3 + d^2).
* - `Neigh.LayerNeighbors`
  - NeighList (Neigh.fastNNnotsquare .false.): interlayer neighbours per atom reserved in the list.
* - `Neigh.LayerDistFactor`
  - NeighList with Neigh.BLInPlaneInteractionRadius: in-plane cutoff 1.1 aG/sqrt(3) times this factor.
* - `Neigh.CompareAll`
  - NeighList: compare every pair of atoms (with the 3 x 3 neighbouring cells). Without it NeighList is known to give wrong bands; a warning is printed.
* - `Neigh.CompareAll.UltraSmall`
  - NeighList: as Neigh.CompareAll with the 7 x 7 neighbouring cells, for cells of fewer than ~36 atoms.
* - `Neigh.BLInPlaneInteractionRadius`
  - NeighList: select neighbours by their in-plane distance (z ignored) within the in-plane cutoff.
* - `Neigh.SafetyCheck`
  - NeighList: list the atoms with fewer than Neigh.ExpectedNumber neighbours.
* - `Neigh.ExpectedNumber`
  - Neigh.SafetyCheck: minimum expected number of neighbours per atom.
* - `Neigh.Print`
  - Print the neighbour list (atom, neighbour, slot, image cell) in the output; atoms without 3 entries are marked.
```

## Calls

{ref}`AtomsSetCart <api-atoms-atomssetcart>`, {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`MPIBarrier <api-mpi-mpibarrier>`, [NeighCountsPrint](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L2708) ({ref}`neigh <api-neigh>`), [NeighIndxCell](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L1321) ({ref}`neigh <api-neigh>`), {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_StringComp <api-string-interface-stringcomp>` (= `StringComp`), {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_PrintDebug <api-sys-sysprintdebug>` (= `SysPrintDebug`), {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

{ref}`HamInit <api-ham-haminit>`
