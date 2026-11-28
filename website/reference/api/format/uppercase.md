---
edit_url: null
---

(api-format-uppercase)=
# `uppercase`

```{index} single: uppercase (function in format)
```

```fortran
elemental function uppercase(str) result(lcstr)
```

- **Module:** {ref}`format <api-format>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/MIO/format.f90, line 283`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/format.f90#L283)

Converts a string to uppercase only

## Parameters

`str` : `character(len=*), intent(in)`
: String to convert

## Returns

`lcstr` : `character(len=maxlinel)`
: Converted string
