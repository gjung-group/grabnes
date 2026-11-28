---
edit_url: null
---

(api-sys)=
# `sys`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/MIO/sys.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/sys.F90#L13)
- **Group:** {ref}`MIO library <api-group-mio>`
- **Uses:** {ref}`format <api-format>`, {ref}`io <api-io>`, {ref}`mpi <api-mpi>`
- **Used by:** {ref}`file <api-file>`, {ref}`input <api-input>`, {ref}`mem <api-mem>`, {ref}`mio <api-mio>`, {ref}`parser <api-parser>`, {ref}`timer <api-timer>`
- **Public interface:** 9 procedures, 1 generic interface

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`SysDebugging <api-sys-sysdebugging>`
  - function
  - True if debug output is on for this module (or for every module).
* - {ref}`SysKill <api-sys-syskill>`
  - subroutine
  - Display an error notification and terminates the program.
* - {ref}`SysKillFlush <api-sys-syskillflush>`
  - subroutine
  - If there are errors, the program stops
* - {ref}`SysPrintDebug <api-sys-sysprintdebug>`
  - subroutine
  - Print a debug line, "module: DEBUG: message", only when debug output is on for the module (Output.Debug, Output.DebugModules) or when 'force' is true (a module's own debug key, e.g.
* - {ref}`SysPrintErrors <api-sys-sysprinterrors>`
  - subroutine
  - Print the count of errors and warnings.
* - {ref}`SysSection <api-sys-syssection>`
  - subroutine
  - Print a section header: a blank line, then one line of '=' that carries the title and the module.
* - {ref}`SysSetDebug <api-sys-syssetdebug>`
  - subroutine
  - Switch debug output on or off, for every module (no module given) or for one module.
* - {ref}`SysSubsection <api-sys-syssubsection>`
  - subroutine
  - Print a sub-section title inside a section: "module: --- title ---".
* - {ref}`SysWarning <api-sys-syswarning>`
  - subroutine
  - Display a warning notification
```

## Generic interfaces

(api-sys-interface-sysprint)=
### `SysPrint`

Specific procedures: `SysPrint_mod`, `SysPrint__`.

Defined in [`lanczosKuboCode/Src/MIO/sys.F90, line 49`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/sys.F90#L49).
