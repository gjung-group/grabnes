---
edit_url: null
---

(api-format-num2bytes)=
# `num2bytes`

```{index} single: num2bytes (function in format)
```

```fortran
elemental function num2bytes(num) result(string)
```

- **Module:** {ref}`format <api-format>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/MIO/format.f90, line 200`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/format.f90#L200)

Converts a number of bytes to a string with the proper representation

## Parameters

`num` : `integer(8), intent(in)`
: Number of bytes

## Returns

`string` : `character(len=30)`
: —

## Calls

{ref}`num2stri <api-format-num2stri>`, {ref}`num2strr <api-format-num2strr>`

## Called by

{ref}`MemPrint <api-mem-memprint>`
