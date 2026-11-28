---
edit_url: null
---

(api-ham-haminputenergyev)=
# `HamInputEnergyEV`

```{index} single: HamInputEnergyEV (subroutine in ham)
```

```fortran
subroutine HamInputEnergyEV(key, value, defaultG0, since)
```

- **Module:** {ref}`ham <api-ham>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/ham.F90, line 13260`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L13260)

Reads an energy that the input file gives in eV and returns it in units of g0. If the input file does not set the key, value = defaultG0, the former default (already in units of g0), so such runs are unchanged. These keys were taken in units of g0 before 2026-10-14 (or before the optional date since); when the input sets one, a warning says so.

## Parameters

`key` : `character(len=*), intent(in)`
: —

`value` : `real(dp), intent(out)`
: —

`defaultG0` : `real(dp), intent(in)`
: —

`since` : `character(len=*), intent(in), optional`
: —

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_InputSearchLabel <api-input-inputsearchlabel>` (= `InputSearchLabel`), {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

[DiagSpectralEpsilon](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2103) ({ref}`diag <api-diag>`), [DiagSpectralWindow](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L2088) ({ref}`diag <api-diag>`), {ref}`HamHopping <api-ham-hamhopping>`, {ref}`HamOnSite <api-ham-hamonsite>`
