---
edit_url: null
---

(api-input-inputsearchlabel)=
# `InputSearchLabel`

```{index} single: InputSearchLabel (function in input)
```

```fortran
function InputSearchLabel(label, str, lineid) result(found)
```

- **Module:** {ref}`input <api-input>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/MIO/input.F90, line 1451`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1451)

Search a key in the input file, taking alternative names into account.

A key can have two names (InputAddAlias): the Section.Name form and the former name. The solver may ask for either, and the input file may contain either. If the file contains both, the Section.Name form is used and the conflict is reported at the end of the run. Labels of blocks (starting with '&') have no alternative names.

## Parameters

`label` : `character(len=*), intent(in)`
: —

`str` : `character(len=*), intent(out), optional`
: —

`lineid` : `integer, intent(out), optional`
: —

## Returns

`found` : `logical`
: —

## Calls

[InputLower](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1543) ({ref}`input <api-input>`), [InputMarkAsked](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1617) ({ref}`input <api-input>`), [InputSearchLabelRaw](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1733) ({ref}`input <api-input>`)

## Called by

{ref}`CellGet <api-cell-cellget>`, {ref}`DiagSpectralFunctionKGridInequivalent <api-diag-diagspectralfunctionkgridinequivalent>`, {ref}`DiagSpectralFunctionKGridInequivalent_v2 <api-diag-diagspectralfunctionkgridinequivalent_v2>`, {ref}`DiagSpectralFunctionKGridInequivalentEnergyCut <api-diag-diagspectralfunctionkgridinequivalentenergycut>`, {ref}`DiagSpectralFunctionKGridInequivalentEnergyCut_v2 <api-diag-diagspectralfunctionkgridinequivalentenergycut_v2>`, {ref}`DiagSpectralFunctionKGridInequivalentEnergyCutNickDale <api-diag-diagspectralfunctionkgridinequivalentenergycutnickdale>`, {ref}`HamInputEnergyEV <api-ham-haminputenergyev>`, [HamInputHBN](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L13105) ({ref}`ham <api-ham>`), {ref}`HamOnSite <api-ham-hamonsite>`, [InputBl\_d](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1155) ({ref}`input <api-input>`), [InputBl\_i](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L946) ({ref}`input <api-input>`), [InputBl\_l](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1264) ({ref}`input <api-input>`), [InputBl\_r](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1048) ({ref}`input <api-input>`), [InputBl\_s](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1366) ({ref}`input <api-input>`), [InputFindBl](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L756) ({ref}`input <api-input>`), [InputFindBl\_i](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L806) ({ref}`input <api-input>`), [InputFindBl\_s](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L868) ({ref}`input <api-input>`), [InputPar\_d](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L490) ({ref}`input <api-input>`), [InputPar\_i](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L343) ({ref}`input <api-input>`), [InputPar\_i\_variable](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L408) ({ref}`input <api-input>`), [InputPar\_l](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L641) ({ref}`input <api-input>`), [InputPar\_r](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L570) ({ref}`input <api-input>`), [InputPar\_s](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L709) ({ref}`input <api-input>`), {ref}`MagfInit <api-magf-magfinit>`, {ref}`TBInit <api-tbpar-tbinit>`
