---
edit_url: null
---

(api-magf-haldphaseinit)=
# `HaldPhaseInit`

```{index} single: HaldPhaseInit (subroutine in magf)
```

```fortran
subroutine HaldPhaseInit()
```

- **Module:** {ref}`magf <api-magf>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/magf.F90, line 148`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/magf.F90#L148)

## Input keys read (6)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Haldane.SpecifyFlux`
  - Haldane model: one phase, 2 pi Haldane.Flux / Phi\_0 (no sweep over HaldPhase.\*).
* - `Haldane.SpecifyPhase`
  - Haldane model: one phase, Haldane.Phase \[rad\] (no sweep over HaldPhase.\*).
* - `Haldane.SpecifyRange`
  - Haldane model: sweep the phase over HaldPhase.InitialInteger..FinalInteger; step m has phase m\*2\*pi/(FinalInteger - InitialInteger).
* - `HaldPhase.InitialInteger`
  - Haldane.SpecifyRange: first step of the phase sweep.
* - `HaldPhase.FinalInteger`
  - Haldane.SpecifyRange: last step of the phase sweep.
* - `HaldPhase.Step`
  - Haldane.SpecifyRange: increment of the step index.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`)

## Called by

{ref}`CalcSelect <api-calc-calcselect>`
