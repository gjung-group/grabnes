---
edit_url: null
---

(api-input)=
# `input`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/MIO/input.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L13)
- **Group:** {ref}`MIO library <api-group-mio>`
- **Uses:** {ref}`debug <api-debug>`, {ref}`file <api-file>`, {ref}`io <api-io>`, {ref}`mem <api-mem>`, {ref}`mpi <api-mpi>`, `omp_lib`, {ref}`options <api-options>`, {ref}`parser <api-parser>`, {ref}`precision <api-precision>`, {ref}`string <api-string>`, {ref}`sys <api-sys>`, {ref}`timer <api-timer>`
- **Used by:** {ref}`mio <api-mio>`
- **Public interface:** 4 procedures, 4 generic interfaces

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`InputAddAlias <api-input-inputaddalias>`
  - subroutine
  - Declare that the key newName (Section.Name form) and the key oldName are the same key.
* - {ref}`InputClose <api-input-inputclose>`
  - subroutine
  - Close input module
* - {ref}`InputInit <api-input-inputinit>`
  - subroutine
  - Initialize the input information
* - {ref}`InputSearchLabel <api-input-inputsearchlabel>`
  - function
  - Search a key in the input file, taking alternative names into account.
```

## Generic interfaces

(api-input-interface-inputblock)=
### `InputBlock`

Get elements from a block

**Arguments (from the source comment)**

`character label`
: Label of the block

`integer column`
: Column number

`* value`
: Elements obtained from the block

Specific procedures: `InputBl_i`, `InputBl_r`, `InputBl_d`, `InputBl_l`, `InputBl_s`, `InputBl_im`, `InputBl_rm`, `InputBl_dm`, `InputBl_lm`, `InputBl_imi`, `InputBl_rmi`, `InputBl_dmi`, `InputBl_lmi`.

Defined in [`lanczosKuboCode/Src/MIO/input.F90, line 72`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L72).

(api-input-interface-inputfindblock)=
### `InputFindBlock`

Finds a block

**Arguments (from the source comment)**

`character label`
: Label of the block

`logical found`
: True if the block is found

Specific procedures: `InputFindBl`, `InputFindBl_i`, `InputFindBl_s`.

Defined in [`lanczosKuboCode/Src/MIO/input.F90, line 68`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L68).

(api-input-interface-inputparameter)=
### `InputParameter`

Returns the parameter with the specified label from input file

**Arguments (from the source comment)**

`character label`
: Label of the parameter

`* def(:)`
: Default value

`character(*) defstr`
: Default string if it's required

`* value`
: Value for the parameter

`character(*) str`
: Optional string found after the values

Specific procedures: `InputPar_i`, `InputPar_d`, `InputPar_r`, `InputPar_l`, `InputPar_s`, `InputPar_is`, `InputPar_ds`, `InputPar_rs`, `InputPar_ls`.

Defined in [`lanczosKuboCode/Src/MIO/input.F90, line 58`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L58).

(api-input-interface-inputparametervariable)=
### `InputParameterVariable`

Specific procedures: `InputPar_i_variable`.

Defined in [`lanczosKuboCode/Src/MIO/input.F90, line 64`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L64).
