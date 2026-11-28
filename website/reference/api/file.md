---
edit_url: null
---

(api-file)=
# `file`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/MIO/file.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L13)
- **Group:** {ref}`MIO library <api-group-mio>`
- **Uses:** {ref}`format <api-format>`, {ref}`io <api-io>`, {ref}`mpi <api-mpi>`, {ref}`sys <api-sys>`
- **Used by:** {ref}`debug <api-debug>`, {ref}`input <api-input>`, {ref}`mio <api-mio>`
- **Public interface:** 1 procedure, 1 derived type

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`CloseAll <api-file-closeall>`
  - subroutine
  - Close all files
```

## Derived types

(api-file-type-cl_file)=
### `cl_file`

```{list-table}
:header-rows: 1

* - Component
  - Declaration
  - Default
  - Description
* - `name`
  - `character(len=maxfnlen)`
  - `'unknown'`
  - —
* - `unit`
  - `integer`
  - —
  - —
* - `form`
  - `character(len=11)`
  - `'FORMATTED'`
  - —
* - `status`
  - `character(len=7)`
  - `'UNKNOWN'`
  - —
* - `access`
  - `character(len=11)`
  - `'SEQUENTIAL'`
  - —
* - `position`
  - `character(len=6)`
  - `'ASIS'`
  - —
* - `action`
  - `character(len=9)`
  - `'READWRITE'`
  - —
* - `delete`
  - `character(len=6)`
  - `'KEEP'`
  - —
* - `recl`
  - `integer`
  - `maxlinel`
  - —
* - `iostat`
  - `integer`
  - `0`
  - —
* - `opened`
  - `logical`
  - `.false.`
  - —
* - `init`
  - `logical`
  - `.false.`
  - —
* - `prev`
  - `class(cl_file), pointer`
  - `NULL()`
  - —
* - `next`
  - `class(cl_file), pointer`
  - `NULL()`
  - —
* - `serial`
  - `logical`
  - `.false.`
  - —
* - `comm`
  - `integer`
  - `MPI_COMM_WORLD`
  - —
* - `amode`
  - `integer`
  - —
  - —
* - `disp`
  - `integer`
  - —
  - —
* - `etype`
  - `integer`
  - —
  - —
* - `filetype`
  - `integer`
  - —
  - —
* - `datarep`
  - `character(len=MPI_MAX_DATAREP_STRING)`
  - —
  - —
```

Defined in [`lanczosKuboCode/Src/MIO/file.F90, line 23`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L23).
