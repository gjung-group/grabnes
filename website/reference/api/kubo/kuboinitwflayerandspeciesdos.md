---
edit_url: null
---

(api-kubo-kuboinitwflayerandspeciesdos)=
# `KuboInitWFLayerAndSpeciesDOS`

```{index} single: KuboInitWFLayerAndSpeciesDOS (subroutine in kubo)
```

```fortran
subroutine KuboInitWFLayerAndSpeciesDOS(Psi, layerNumber, speciesNumber, numberOfLayers, numberOfAtomsInLayer)
```

- **Module:** {ref}`kubo <api-kubo>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubo.F90, line 214`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubo.F90#L214)

## Parameters

`Psi` : `complex(dp), dimension(inode1:), intent(out)`
: —

`layerNumber` : `integer`
: —

`speciesNumber` : `integer`
: —

`numberOfLayers` : `integer`
: —

`numberOfAtomsInLayer` : `integer`
: —

## Input keys read (1)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Kubo.SeedSet`
  - Seed of the random-phase start vector of a layer-and-species DOS (Kubo.PDOSLayerAndSpecies).
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)

## Called by

[CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`)
