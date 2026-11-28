---
edit_url: null
---

(api-string)=
# `string`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/MIO/string.f90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/string.f90#L13)
- **Group:** {ref}`MIO library <api-group-mio>`
- **Uses:** {ref}`io <api-io>`
- **Used by:** {ref}`input <api-input>`, {ref}`mio <api-mio>`, {ref}`parser <api-parser>`
- **Public interface:** 1 generic interface

## Generic interfaces

(api-string-interface-stringcomp)=
### `StringComp`

Compare two strings and check whether they are equal

**Arguments (from the source comment)**

`character string1, string2`
: Strings to be compared

`logical equal`
: True if the strings are equal

Specific procedures: `StringComp_s`, `StringComp_v`.

Defined in [`lanczosKuboCode/Src/MIO/string.f90, line 23`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/string.f90#L23).
