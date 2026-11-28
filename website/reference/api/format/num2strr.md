---
edit_url: null
---

(api-format-num2strr)=
# `num2strr`

```{index} single: num2strr (function in format)
```

```fortran
elemental function num2strr(num, dec) result(string)
```

- **Module:** {ref}`format <api-format>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/MIO/format.f90, line 63`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/format.f90#L63)

## Parameters

`num` : `real(sp), intent(in)`
: —

`dec` : `integer, intent(in), optional`
: —

## Returns

`string` : `character(len=strlen)`
: —

## Called by

{ref}`num2bytes <api-format-num2bytes>`
