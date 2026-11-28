---
edit_url: null
---

(api-random-random_old-randinit)=
# `RandInit`

```{index} single: RandInit (subroutine in random)
```

```fortran
subroutine RandInit()
```

- **Module:** {ref}`random <api-random-random_old>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/random_old.F90, line 24`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/random_old.F90#L24)

## Input keys read (2)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Random.Test`
  - Old random-number module (random\_old.F90, not compiled): run a test of the generator.
* - `Random.TestNumber`
  - Random.Test: number of random numbers drawn.
```

## Calls

[GetUnit](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/file.F90#L606) ({ref}`file <api-file>`), {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`RandGen <api-random-random_old-randgen>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)
