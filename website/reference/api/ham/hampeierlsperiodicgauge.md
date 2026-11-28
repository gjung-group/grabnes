---
edit_url: null
---

(api-ham-hampeierlsperiodicgauge)=
# `HamPeierlsPeriodicGauge`

```{index} single: HamPeierlsPeriodicGauge (subroutine in ham)
```

```fortran
subroutine HamPeierlsPeriodicGauge(flux)
```

- **Module:** {ref}`ham <api-ham>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/ham.F90, line 12834`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L12834)

Peierls phases of the periodic gauge (Cresti) for a field of flux = B pi / flux quantum \[1/Ang^2\], applied to every hopping, with the corner term exp(i pi m n1 n2). Used by HamHopping and by the Output.ReadDataFiles path of CalcSelect, so that both get the same phases. Positions must be Cartesian.

## Parameters

`flux` : `real(dp), intent(in)`
: —

## Input keys read (1)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `MagField.LegacyCornerPhase`
  - Peierls phases without the exp(i pi m n1 n2) term of bonds across a cell corner, as before 2026-10-22 (wrong for odd m). The run otherwise refuses a non-integer number of flux quanta.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

{ref}`CalcSelect <api-calc-calcselect>`, {ref}`HamHopping <api-ham-hamhopping>`
