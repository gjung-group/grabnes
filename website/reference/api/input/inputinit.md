---
edit_url: null
---

(api-input-inputinit)=
# `InputInit`

```{index} single: InputInit (subroutine in input)
```

```fortran
subroutine InputInit(success)
```

- **Module:** {ref}`input <api-input>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/input.F90, line 95`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L95)

Initialize the input information

## Parameters

`success` : `logical, intent(out)`
: —

## Calls

{ref}`DebugPrint <api-debug-debugprint>`, [GetFH](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L634) ({ref}`file <api-file>`), [GetUnit](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L606) ({ref}`file <api-file>`), [InputAliasInit](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1585) ({ref}`input <api-input>`), [InputCheckDuplicates](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L208) ({ref}`input <api-input>`), [InputProcess](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1827) ({ref}`input <api-input>`), {ref}`OptionsInputName <api-options-optionsinputname>`, [SetDelete](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L584) ({ref}`file <api-file>`), {ref}`SysWarning <api-sys-syswarning>`, {ref}`TimerCount <api-timer-timercount>`, {ref}`TimerStop <api-timer-timerstop>`

## Called by

{ref}`MIO_Initialize <api-mio-mio_initialize>`
