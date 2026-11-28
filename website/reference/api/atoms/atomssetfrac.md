---
edit_url: null
---

(api-atoms-atomssetfrac)=
# `AtomsSetFrac`

```{index} single: AtomsSetFrac (subroutine in atoms)
```

```fortran
subroutine AtomsSetFrac()
```

- **Module:** {ref}`atoms <api-atoms>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/atoms.F90, line 1790`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/atoms.F90#L1790)

## Calls

{ref}`MIO_Debug <api-mio-interface-mio_debug>`

## Called by

{ref}`AtomsPos <api-atoms-atomspos>`, [CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`), [DiagH0TAPW](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L5902) ({ref}`diag <api-diag>`), [DiagH0TAPW\_withBlockH](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L7289) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKu](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11161) ({ref}`diag <api-diag>`), {ref}`HamHopping <api-ham-hamhopping>`, {ref}`HamInit <api-ham-haminit>`, {ref}`HamOnSite <api-ham-hamonsite>`
