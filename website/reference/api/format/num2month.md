---
edit_url: null
---

(api-format-num2month)=
# `num2month`

```{index} single: num2month (function in format)
```

```fortran
elemental function num2month(num) result(string)
```

- **Module:** {ref}`format <api-format>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/MIO/format.f90, line 169`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/format.f90#L169)

Converts a number to a string with the corresponding month

## Parameters

`num` : `integer, intent(in)`
: The number of month

## Returns

`string` : `character(len=9)`
: —

## Called by

{ref}`MIO_Finalize <api-mio-mio_finalize>`, {ref}`MIO_Initialize <api-mio-mio_initialize>`
