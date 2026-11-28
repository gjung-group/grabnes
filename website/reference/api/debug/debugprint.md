---
edit_url: null
---

(api-debug-debugprint)=
# `DebugPrint`

```{index} single: DebugPrint (subroutine in debug)
```

```fortran
subroutine DebugPrint(subroutine, i)
```

- **Module:** {ref}`debug <api-debug>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/debug.F90, line 62`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/debug.F90#L62)

Write debug information. All debug information is written to file debug.log. This routine prints a report on the routines that are being used to follow track of the execution of the program.

## Parameters

`subroutine` : `character(len=*), intent(in)`
: Calling subroutine

`i` : `integer, intent(in)`
: 0 - Entering routine otherwise - Routine finished

## Called by

{ref}`InputInit <api-input-inputinit>`
