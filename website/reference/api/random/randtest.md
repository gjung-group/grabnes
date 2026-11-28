---
edit_url: null
---

(api-random-randtest)=
# `RandTest`

```{index} single: RandTest (subroutine in random)
```

```fortran
subroutine RandTest(rng, n)
```

- **Module:** {ref}`random <api-random>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/random.F90, line 74`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/random.F90#L74)

## Parameters

`rng` : `type(rand_t), intent(inout)`
: —

`n` : `integer, intent(in)`
: —

## Calls

{ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`num2str <api-mio-interface-num2str>`, [RandNum](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/random.F90#L54) ({ref}`random <api-random>`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)
