---
edit_url: null
---

(api-timer-timerstop)=
# `TimerStop`

```{index} single: TimerStop (subroutine in timer)
```

```fortran
subroutine TimerStop(module)
```

- **Module:** {ref}`timer <api-timer>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/timer.F90, line 121`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/timer.F90#L121)

Stop the timer for a specified module

## Parameters

`module` : `character(len=*), intent(in)`
: Name of the module

## Calls

{ref}`SysKill <api-sys-syskill>`

## Called by

[InputBl\_d](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1155) ({ref}`input <api-input>`), [InputBl\_i](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L946) ({ref}`input <api-input>`), [InputBl\_l](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1264) ({ref}`input <api-input>`), [InputBl\_r](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1048) ({ref}`input <api-input>`), [InputBl\_s](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1366) ({ref}`input <api-input>`), {ref}`InputClose <api-input-inputclose>`, [InputFindBl](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L756) ({ref}`input <api-input>`), [InputFindBl\_i](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L806) ({ref}`input <api-input>`), [InputFindBl\_s](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L868) ({ref}`input <api-input>`), {ref}`InputInit <api-input-inputinit>`, [InputPar\_d](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L490) ({ref}`input <api-input>`), [InputPar\_i](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L343) ({ref}`input <api-input>`), [InputPar\_i\_variable](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L408) ({ref}`input <api-input>`), [InputPar\_l](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L641) ({ref}`input <api-input>`), [InputPar\_r](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L570) ({ref}`input <api-input>`), [InputPar\_s](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L709) ({ref}`input <api-input>`), {ref}`TimerPrint <api-timer-timerprint>`
