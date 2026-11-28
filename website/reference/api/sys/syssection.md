---
edit_url: null
---

(api-sys-syssection)=
# `SysSection`

```{index} single: SysSection (subroutine in sys)
```

```fortran
subroutine SysSection(title, module)
```

- **Module:** {ref}`sys <api-sys>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/sys.F90, line 363`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/sys.F90#L363)

Print a section header: a blank line, then one line of '=' that carries the title and the module. One header per stage of a run (system set-up, neighbour search, Hamiltonian, band structure, ...), so that a reader can find the stage in the output at a glance.

## Parameters

`title` : `character(len=*), intent(in)`
: Title of the section

`module` : `character(len=*), intent(in)`
: Module that called this subroutine

## Calls

{ref}`SysPrint <api-sys-interface-sysprint>`

## Called by

{ref}`AtomsPos <api-atoms-atomspos>`, [CalcDiag](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L678) ({ref}`calc <api-calc>`), [CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`), {ref}`CellGet <api-cell-cellget>`, {ref}`HaldPhaseValue <api-magf-haldphasevalue>`, {ref}`HamInit <api-ham-haminit>`, {ref}`MagfValue <api-magf-magfvalue>`, {ref}`MIO_Finalize <api-mio-mio_finalize>`, {ref}`moireBLShiftValue <api-moireblshift-moireblshiftvalue>`, {ref}`SCFGetCharge <api-scf-scfgetcharge>`, [StartHeader](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/grabnes_start.F90#L53) ({ref}`grabnes_start <api-grabnes_start>`), {ref}`TBInit <api-tbpar-tbinit>`
