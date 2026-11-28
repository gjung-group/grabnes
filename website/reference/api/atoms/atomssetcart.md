---
edit_url: null
---

(api-atoms-atomssetcart)=
# `AtomsSetCart`

```{index} single: AtomsSetCart (subroutine in atoms)
```

```fortran
subroutine AtomsSetCart()
```

- **Module:** {ref}`atoms <api-atoms>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/atoms.F90, line 1762`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/atoms.F90#L1762)

## Calls

{ref}`MIO_Debug <api-mio-interface-mio_debug>`

## Called by

[ApplySpinFlipSOC](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L16170) ({ref}`diag <api-diag>`), {ref}`AtomsPos <api-atoms-atomspos>`, [CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`), {ref}`CalcSelect <api-calc-calcselect>`, [CalcTunn](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L595) ({ref}`calc <api-calc>`), [DiagH0TAPW](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L5902) ({ref}`diag <api-diag>`), [DiagH0TAPW\_withBlockH](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L7289) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalent](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11412) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalentLee](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11693) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalentMoreOrbitals](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11542) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalentNishi](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11866) ({ref}`diag <api-diag>`), [DiagSpectralWeightWeiKuInequivalentOld](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L11286) ({ref}`diag <api-diag>`), {ref}`fastNNnotsquareNotRectangle <api-neigh-fastnnnotsquarenotrectangle>`, {ref}`GaussHeight <api-gauss-gaussheight>`, {ref}`GaussPot <api-gauss-gausspot>`, {ref}`GaussPotDefinedPositions <api-gauss-gausspotdefinedpositions>`, [HamHopMoireFourier](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5009) ({ref}`ham <api-ham>`), {ref}`HamHopping <api-ham-hamhopping>`, {ref}`HamInit <api-ham-haminit>`, {ref}`HamOnSite <api-ham-hamonsite>`, [HamOnSiteMoireFourier](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4978) ({ref}`ham <api-ham>`), {ref}`InterfacePot2 <api-interface-interfacepot2>`, {ref}`NeighList <api-neigh-neighlist>`, [NeighListOld](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L1023) ({ref}`neigh <api-neigh>`)
