---
edit_url: null
---

(api-format-num2strd)=
# `num2strd`

```{index} single: num2strd (function in format)
```

```fortran
elemental function num2strd(num, dec) result(string)
```

- **Module:** {ref}`format <api-format>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/MIO/format.f90, line 105`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/format.f90#L105)

## Parameters

`num` : `real(dp), intent(in)`
: —

`dec` : `integer, intent(in), optional`
: —

## Returns

`string` : `character(len=strlen)`
: —
