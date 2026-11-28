---
edit_url: null
---

(api-random-randseedfrominput)=
# `RandSeedFromInput`

```{index} single: RandSeedFromInput (subroutine in random)
```

```fortran
subroutine RandSeedFromInput()
```

- **Module:** {ref}`random <api-random>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/random.F90, line 31`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/random.F90#L31)

Seed the intrinsic random-number generator for the disorder models.

With setSeed .true. the seed is seedValue, so that a disordered Hamiltonian can be reproduced; otherwise it comes from the system clock, as before.

## Input keys read (2)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Kubo.SetSeed`
  - Here: also seeds the random numbers of the disorder models (RandSeedFromInput) with Kubo.SeedValue.
* - `Kubo.SeedValue`
  - Here: seed of the random numbers of the disorder models when Kubo.SetSeed is set.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`)

## Called by

{ref}`GaussHeight <api-gauss-gaussheight>`, {ref}`GaussPot <api-gauss-gausspot>`, {ref}`GaussPotDefinedPositions <api-gauss-gausspotdefinedpositions>`, {ref}`HamOnSite <api-ham-hamonsite>`
