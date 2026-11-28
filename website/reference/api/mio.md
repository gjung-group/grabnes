---
edit_url: null
---

(api-mio)=
# `mio`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/MIO/mio.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mio.F90#L13)
- **Group:** {ref}`MIO library <api-group-mio>`
- **Uses:** {ref}`debug <api-debug>`, {ref}`file <api-file>`, {ref}`format <api-format>`, {ref}`input <api-input>`, {ref}`io <api-io>`, {ref}`mem <api-mem>`, {ref}`mpi <api-mpi>`, {ref}`omp <api-omp>`, {ref}`options <api-options>`, {ref}`parser <api-parser>`, {ref}`precision <api-precision>`, {ref}`string <api-string>`, {ref}`sys <api-sys>`, {ref}`timer <api-timer>`
- **Used by:** {ref}`atoms <api-atoms>`, {ref}`calc <api-calc>`, {ref}`cell <api-cell>`, {ref}`constants <api-constants>`, {ref}`diag <api-diag>`, {ref}`gauss <api-gauss>`, {ref}`grabnes_end <api-grabnes_end>`, {ref}`grabnes_start <api-grabnes_start>`, {ref}`ham <api-ham>`, {ref}`hybrid <api-hybrid>`, {ref}`interface <api-interface>`, {ref}`kubo <api-kubo>`, {ref}`kuboarrays <api-kuboarrays>`, {ref}`kubosubs <api-kubosubs>`, {ref}`magf <api-magf>`, {ref}`moireBLShift <api-moireblshift>`, {ref}`neigh <api-neigh>`, {ref}`parallel <api-parallel>`, {ref}`random <api-random>`, {ref}`scf <api-scf>`, {ref}`tbpar <api-tbpar>`
- **Public interface:** 2 procedures, 7 generic interfaces

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`MIO_Finalize <api-mio-mio_finalize>`
  - subroutine
  - Finalize MIO library
* - {ref}`MIO_Initialize <api-mio-mio_initialize>`
  - subroutine
  - Initialize MIO library
```

## Generic interfaces

(api-mio-interface-mio_debug)=
### `MIO_Debug`

Specific procedures: `DebugPrint`.

Defined in [`lanczosKuboCode/Src/MIO/mio.F90, line 46`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mio.F90#L46).

(api-mio-interface-mio_filecloseall)=
### `MIO_FileCloseAll`

Specific procedures: `CloseAll`.

Defined in [`lanczosKuboCode/Src/MIO/mio.F90, line 54`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mio.F90#L54).

(api-mio-interface-mio_timercount)=
### `MIO_TimerCount`

Specific procedures: `TimerCount`.

Defined in [`lanczosKuboCode/Src/MIO/mio.F90, line 96`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mio.F90#L96).

(api-mio-interface-mio_timerstop)=
### `MIO_TimerStop`

Specific procedures: `TimerStop`.

Defined in [`lanczosKuboCode/Src/MIO/mio.F90, line 101`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mio.F90#L101).

(api-mio-interface-num2bytes)=
### `num2bytes`

Specific procedures: `num2bytes`.

Defined in [`lanczosKuboCode/Src/MIO/mio.F90, line 70`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mio.F90#L70).

(api-mio-interface-num2month)=
### `num2month`

Specific procedures: `num2month`.

Defined in [`lanczosKuboCode/Src/MIO/mio.F90, line 65`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mio.F90#L65).

(api-mio-interface-num2str)=
### `num2str`

Specific procedures: `num2stri`, `num2strd`, `num2strz`, `num2strr`, `num2strc`.

Defined in [`lanczosKuboCode/Src/MIO/mio.F90, line 60`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mio.F90#L60).
