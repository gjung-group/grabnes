---
edit_url: null
---

(api-format-num2strz)=
# `num2strz`

```{index} single: num2strz (function in format)
```

```fortran
elemental function num2strz(num, dec) result(string)
```

- **Module:** {ref}`format <api-format>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/MIO/format.f90, line 142`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/format.f90#L142)

## Parameters

`num` : `complex(dp), intent(in)`
: —

`dec` : `integer, intent(in), optional`
: —

## Returns

`string` : `character(len=strlen)`
: —

## Calls

{ref}`num2string <api-format-interface-num2string>`
