---
edit_url: null
---

(api-options-optionsinputname)=
# `OptionsInputName`

```{index} single: OptionsInputName (function in options)
```

```fortran
function OptionsInputName() result(inname)
```

- **Module:** {ref}`options <api-options>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/MIO/options.f90, line 63`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/options.f90#L63)

Returns the input file name found in input options

## Returns

`inname` : `character(len=maxfnlen)`
: —

## Called by

{ref}`InputInit <api-input-inputinit>`
