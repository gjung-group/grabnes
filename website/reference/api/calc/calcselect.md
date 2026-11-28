---
edit_url: null
---

(api-calc-calcselect)=
# `CalcSelect`

```{index} single: CalcSelect (subroutine in calc)
```

```fortran
subroutine CalcSelect()
```

- **Module:** {ref}`calc <api-calc>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/calc.F90, line 13`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L13)

## Input keys read (6)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `moireBLShift.Calc`
  - Sweep over rigid layer shifts (Moire.BLShiftDefinedValues or Moire.BLShiftRange): the on-site terms and the Kubo / Diag calculations are repeated for every shift.
* - `Kubo.Calc`
  - Run the recursion calculation (CalcKubo): DOS or local DOS, then the time evolution unless Calculate.OnlyDOS.
* - `Diag.Calc`
  - Run the diagonalisation calculations (CalcDiag): bands, DOS, spectral function, selected by the Calculate.\* keys.
* - `Output.ReadDataFiles`
  - Take the hoppings from the file generate.s (one 're im' pair per line, units of g0) instead of computing them; the Peierls phases of a magnetic field are added after reading.
* - `Tunn.Calc`
  - Write the interlayer sums T\_AA, T\_AB = -sum\_j t\_ij exp(i k.d\_ij) from each Species-1 atom of layer 1 to layer 2, at k = Tunn.BigKVecX/Y/Z (CalcTunn; files \<prefix>.TunnAAR, TunnAAI, TunnABR, TunnABI).
* - `MagField.FrankMagneticField`
  - Alternative Peierls gauge of HamHopMagneticField (ham.F90); read here but not used in this routine.
```

## Calls

{ref}`AtomsSetCart <api-atoms-atomssetcart>`, [CalcDiag](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L678) ({ref}`calc <api-calc>`), [CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`), [CalcTunn](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L595) ({ref}`calc <api-calc>`), {ref}`HaldPhaseInit <api-magf-haldphaseinit>`, {ref}`HaldPhaseValue <api-magf-haldphasevalue>`, {ref}`HamHopping <api-ham-hamhopping>`, {ref}`HamInit <api-ham-haminit>`, {ref}`HamOnSite <api-ham-hamonsite>`, {ref}`HamPeierlsPeriodicGauge <api-ham-hampeierlsperiodicgauge>`, {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MagfInit <api-magf-magfinit>`, {ref}`MagfValue <api-magf-magfvalue>`, {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`moireBLShiftInit <api-moireblshift-moireblshiftinit>`, {ref}`moireBLShiftValue <api-moireblshift-moireblshiftvalue>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_PrintDebug <api-sys-sysprintdebug>` (= `SysPrintDebug`), {ref}`TBInit <api-tbpar-tbinit>`

## Called by

{ref}`grabnes <api-grabnes>` (program)
