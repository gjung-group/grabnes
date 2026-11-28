---
edit_url: null
---

(api-format-lowercase)=
# `lowercase`

```{index} single: lowercase (function in format)
```

```fortran
elemental function lowercase(str) result(lcstr)
```

- **Module:** {ref}`format <api-format>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/MIO/format.f90, line 244`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/format.f90#L244)

Converts a string to lowercase only

## Parameters

`str` : `character(len=*), intent(in)`
: String to convert

## Returns

`lcstr` : `character(len=maxlinel)`
: Converted string

## Called by

[ParserD](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/parser.f90#L417) ({ref}`parser <api-parser>`), [ParserInt](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/parser.f90#L336) ({ref}`parser <api-parser>`), [ParserL](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/parser.f90#L463) ({ref}`parser <api-parser>`), [ParserR](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/parser.f90#L386) ({ref}`parser <api-parser>`)
