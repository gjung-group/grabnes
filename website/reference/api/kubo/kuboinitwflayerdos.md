---
edit_url: null
---

(api-kubo-kuboinitwflayerdos)=
# `KuboInitWFLayerDOS`

```{index} single: KuboInitWFLayerDOS (subroutine in kubo)
```

```fortran
subroutine KuboInitWFLayerDOS(Psi, layerNumber, numberOfLayers, numberOfAtomsInLayer)
```

- **Module:** {ref}`kubo <api-kubo>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubo.F90, line 140`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubo.F90#L140)

## Parameters

`Psi` : `complex(dp), dimension(inode1:), intent(out)`
: —

`layerNumber` : `integer`
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
  - Seed of the random-phase start vector of a layer DOS (Kubo.PDOSLayer); always used, not the clock.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Warning <api-sys-syswarning>` (= `SysWarning`)

## Called by

[CalcKubo](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/calc.F90#L123) ({ref}`calc <api-calc>`)
