---
edit_url: null
---

(api-ham-socenabledforlayer)=
# `SOCEnabledForLayer`

```{index} single: SOCEnabledForLayer (function in ham)
```

```fortran
function SOCEnabledForLayer(layerIndex) result(enabled)
```

- **Module:** {ref}`ham <api-ham>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/ham.F90, line 13602`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L13602)

## Parameters

`layerIndex` : `integer, intent(in)`
: —

## Returns

`enabled` : `logical`
: —

## Called by

[ApplySOCtoHamiltonian](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L16060) ({ref}`diag <api-diag>`), [ApplySpinFlipSOC](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L16170) ({ref}`diag <api-diag>`), [compute\_block\_hamiltonian\_derivatives](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L13471) ({ref}`diag <api-diag>`), [HamSOCBondTables](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L12929) ({ref}`ham <api-ham>`), [initialize\_sparse\_matrix](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L10457) ({ref}`diag <api-diag>`), [SOCOnsiteAdd](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L16660) ({ref}`diag <api-diag>`), [transform\_dense\_hamiltonian\_tapw](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L15335) ({ref}`diag <api-diag>`), {ref}`ZeemanEnabledForLayer <api-ham-zeemanenabledforlayer>`
