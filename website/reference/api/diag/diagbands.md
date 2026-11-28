---
edit_url: null
---

(api-diag-diagbands)=
# `DiagBands`

```{index} single: DiagBands (subroutine in diag)
```

```fortran
subroutine DiagBands()
```

- **Module:** {ref}`diag <api-diag>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/diag.F90, line 602`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L602)

## Input keys read (13)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Bands.NumPoints`
  - k-points on the first leg of Bands.Path; the other legs get points in proportion to their length (all the same with Bands.UseSameNumberOfPoints).
* - `Bands.SparseNeig`
  - Diag.SparseSolver: number of eigenvalues computed per k-point.
* - `Bands.UseDifferentLatticeVectors`
  - Set ucell(1,2) and rcell(2,1) to zero before the path is built (rectangular-cell variant; see DiagBands).
* - `Bands.UseSameNumberOfPoints`
  - Every leg of Bands.Path gets Bands.NumPoints points, whatever its length.
* - `Output.KeepWaveFunction`
  - Dense diagonalisation that also writes the eigenvectors. Refused together with TAPW.Use or Diag.SparseSolver (it replaced those paths silently before).
* - `Diag.SparseSolver`
  - Bands by sparse ARPACK diagonalisation (DiagHamSparse): Bands.SparseNeig levels nearest E = 0, or nearest Diag.SparseSetShift with Diag.SparseUseShift. Cannot be combined with TAPW.Use.
* - `TAPW.Use`
  - Bands in the valley-projected plane-wave basis (TAPW): H is projected on plane waves around one graphene valley (Diag.N\_G, Diag.MoireAngle, Diag.UseKprimeValley); far fewer states than atoms.
* - `TAPW.UseDenseMatrix`
  - TAPW: project a dense N x N Hamiltonian (memory 16 N^2 bytes) instead of the sparse one; the sparse path is the one for large cells.
* - `TAPW.BlockMatrixFree`
  - Spinful block TAPW without the dense 2N x 2N block (matrix-free H\*X): much faster, little memory.
* - `Zeeman.Term`
  - Zeeman shift of the two spin sectors (needs Zeeman.SpinPolarized .true.; no effect, with a warning, in a one-spin run).
* - `Zeeman.PseudoTerm`
  - Pseudo-Zeeman term: not implemented, no effect (a warning is printed).
* - `TAPW.aG`
  - TAPW: graphene lattice constant of the plane-wave basis \[Angstrom\]; default 2.46019, the cells often use 2.46.
* - `Diag.SparseSetShift`
  - Here: copied for the TAPW path, which does not use it at present.
```

## Calls

[add\_high\_symmetry\_refinement](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L15864) ({ref}`diag <api-diag>`), [BuildBlockHamiltonianOnly](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L16414) ({ref}`diag <api-diag>`), [CalculateChernTAPW](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L12066) ({ref}`diag <api-diag>`), [DiagBlockHamiltonian](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L16592) ({ref}`diag <api-diag>`), [DiagCheckThreads](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L584) ({ref}`diag <api-diag>`), [DiagH0TAPW](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L5902) ({ref}`diag <api-diag>`), [DiagH0TAPW\_withBlockH](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L7289) ({ref}`diag <api-diag>`), [DiagHam](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L5128) ({ref}`diag <api-diag>`), [DiagHamSparse](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L5180) ({ref}`diag <api-diag>`), [DiagHamWF](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L10710) ({ref}`diag <api-diag>`), {ref}`MIO_InputBlock <api-input-interface-inputblock>` (= `InputBlock`), {ref}`MIO_InputFindBlock <api-input-interface-inputfindblock>` (= `InputFindBlock`), {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_PrintDebug <api-sys-sysprintdebug>` (= `SysPrintDebug`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`), {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

[CalcDiag](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L678) ({ref}`calc <api-calc>`)
