---
edit_url: null
---

(api-format)=
# `format`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/MIO/format.f90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/format.f90#L13)
- **Group:** {ref}`MIO library <api-group-mio>`
- **Uses:** {ref}`io <api-io>`, {ref}`precision <api-precision>`
- **Used by:** {ref}`file <api-file>`, {ref}`mem <api-mem>`, {ref}`mio <api-mio>`, {ref}`parser <api-parser>`, {ref}`sys <api-sys>`
- **Public interface:** 9 procedures, 1 generic interface

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`lowercase <api-format-lowercase>`
  - function
  - Converts a string to lowercase only
* - {ref}`num2bytes <api-format-num2bytes>`
  - function
  - Converts a number of bytes to a string with the proper representation
* - {ref}`num2month <api-format-num2month>`
  - function
  - Converts a number to a string with the corresponding month
* - {ref}`num2strc <api-format-num2strc>`
  - function
  - —
* - {ref}`num2strd <api-format-num2strd>`
  - function
  - —
* - {ref}`num2stri <api-format-num2stri>`
  - function
  - —
* - {ref}`num2strr <api-format-num2strr>`
  - function
  - —
* - {ref}`num2strz <api-format-num2strz>`
  - function
  - —
* - {ref}`uppercase <api-format-uppercase>`
  - function
  - Converts a string to uppercase only
```

## Generic interfaces

(api-format-interface-num2string)=
### `num2string`

Converts a number to a character string. Integer, real and complex version are available (in single and double precision for real and complex). The number of decimal digits can be chosen for real and complex conversions, the precision should be between 1 and 15.

**Arguments (from the source comment)**

`* num`
: The number to be converted

`integer dec`
: Decimal precision for real and complex numbers

Specific procedures: {ref}`num2stri <api-format-num2stri>`, {ref}`num2strd <api-format-num2strd>`, {ref}`num2strz <api-format-num2strz>`, {ref}`num2strr <api-format-num2strr>`, {ref}`num2strc <api-format-num2strc>`.

Defined in [`lanczosKuboCode/Src/MIO/format.f90, line 33`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/format.f90#L33).
