---
edit_url: null
---

(api-parser)=
# `parser`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/MIO/parser.f90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/parser.f90#L13)
- **Group:** {ref}`MIO library <api-group-mio>`
- **Uses:** {ref}`format <api-format>`, {ref}`io <api-io>`, {ref}`precision <api-precision>`, {ref}`string <api-string>`, {ref}`sys <api-sys>`
- **Used by:** {ref}`input <api-input>`, {ref}`mio <api-mio>`
- **Public interface:** 2 generic interfaces

## Generic interfaces

(api-parser-interface-parsercheck)=
### `ParserCheck`

Check if some value within a string is correct for its type

**Arguments (from the source comment)**

`character string`
: String where to look for

`integer enum`
: Number of the element inside string

`* value`
: Value obtained from the string

`integer err`
: Number of error

Specific procedures: `ParserCheck_i`, `ParserCheck_r`, `ParserCheck_d`, `ParserCheck_s`, `ParserCheck_l`.

Defined in [`lanczosKuboCode/Src/MIO/parser.f90, line 22`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/parser.f90#L22).

(api-parser-interface-parserfindvalue)=
### `ParserFindValue`

Find a value in a string. The expression is of the form: 'name = value'

**Arguments (from the source comment)**

`character name`
: Name of the value to find

`character string`
: String where to look for the value

`* value`
: The value found

`logical found`
: True if the value is found

Specific procedures: `ParserFindValue_i`, `ParserFindValue_d`.

Defined in [`lanczosKuboCode/Src/MIO/parser.f90, line 27`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/parser.f90#L27).
