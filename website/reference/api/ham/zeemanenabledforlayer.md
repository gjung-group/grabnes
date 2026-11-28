---
edit_url: null
---

(api-ham-zeemanenabledforlayer)=
# `ZeemanEnabledForLayer`

```{index} single: ZeemanEnabledForLayer (function in ham)
```

```fortran
function ZeemanEnabledForLayer(layerIndex) result(enabled)
```

- **Module:** {ref}`ham <api-ham>`
- **Kind:** Function
- **Source:** [`lanczosKuboCode/Src/ham.F90, line 13589`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L13589)

Check if SOC should be applied to a specific layer

Layers that carry the Zeeman term: the SOC layers when SOC.LayerControl is on and Zeeman.FollowSOCLayers is .true. (default), otherwise every layer

## Parameters

`layerIndex` : `integer, intent(in)`
: Layer index to check

## Returns

`enabled` : `logical`
: True if SOC should be applied to this layer

## Calls

{ref}`SOCEnabledForLayer <api-ham-socenabledforlayer>`

## Called by

[ApplySOCtoHamiltonian](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L16060) ({ref}`diag <api-diag>`), [initialize\_sparse\_matrix](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L10457) ({ref}`diag <api-diag>`), [SOCOnsiteAdd](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L16660) ({ref}`diag <api-diag>`), [transform\_dense\_hamiltonian\_tapw](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L15335) ({ref}`diag <api-diag>`)
