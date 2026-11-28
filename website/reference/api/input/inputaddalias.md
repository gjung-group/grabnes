---
edit_url: null
---

(api-input-inputaddalias)=
# `InputAddAlias`

```{index} single: InputAddAlias (subroutine in input)
```

```fortran
subroutine InputAddAlias(newName, oldName)
```

- **Module:** {ref}`input <api-input>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/input.F90, line 1518`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1518)

Declare that the key newName (Section.Name form) and the key oldName are the same key.

## Parameters

`newName` : `character(len=*), intent(in)`
: —

`oldName` : `character(len=*), intent(in)`
: —

## Calls

[InputLower](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1543) ({ref}`input <api-input>`), {ref}`SysKill <api-sys-syskill>`

## Called by

[InputAliasInit](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/input.F90#L1585) ({ref}`input <api-input>`)
