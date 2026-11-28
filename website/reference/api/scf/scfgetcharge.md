---
edit_url: null
---

(api-scf-scfgetcharge)=
# `SCFGetCharge`

```{index} single: SCFGetCharge (subroutine in scf)
```

```fortran
subroutine SCFGetCharge()
```

- **Module:** {ref}`scf <api-scf>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/scf.F90, line 20`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/scf.F90#L20)

## Input keys read (8)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Diag.KGrid`
  - Number of k-points along each reciprocal lattice vector (three integers) of the k-grid; here of the Hubbard SCF.
* - `Hubbard.SCFMaxIterations`
  - Hubbard SCF: maximum number of iterations.
* - `Hubbard.SCFTolerance`
  - Hubbard SCF: convergence threshold on the change of the charges.
* - `Hubbard.SCFMix`
  - Hubbard SCF: linear mixing, charge = (1 - mix) new + mix old.
* - `Hubbard.MOCharge`
  - Hubbard SCF: spin imbalance dq added to (and subtracted from) the starting charges of the ordered pattern.
* - `Hubbard.AFOrdering`
  - Hubbard SCF: start from an antiferromagnetic (sublattice) spin pattern.
* - `Hubbard.AFOrderingLayer`
  - Hubbard SCF: start from a pattern antiferromagnetic between layers (if Hubbard.AFOrdering is off).
* - `Hubbard.FOrderingLayer`
  - Hubbard SCF: start from a layer-ferromagnetic pattern (if the two orderings above are off).
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, [SCFHam](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/scf.F90#L216) ({ref}`scf <api-scf>`), [SCFHam0](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/scf.F90#L179) ({ref}`scf <api-scf>`), {ref}`SCFInit <api-scf-scfinit>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Section <api-sys-syssection>` (= `SysSection`)

## Called by

[CalcDiag](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L678) ({ref}`calc <api-calc>`)
