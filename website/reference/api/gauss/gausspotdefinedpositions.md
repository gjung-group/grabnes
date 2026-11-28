---
edit_url: null
---

(api-gauss-gausspotdefinedpositions)=
# `GaussPotDefinedPositions`

```{index} single: GaussPotDefinedPositions (subroutine in gauss)
```

```fortran
subroutine GaussPotDefinedPositions(H)
```

- **Module:** {ref}`gauss <api-gauss>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/gauss.F90, line 145`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/gauss.F90#L145)

## Parameters

`H` : `real(dp), dimension(inode1:), intent(inout)`
: —

## Input keys read (2)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Strain.BubbleGaussPotRange`
  - Strain.BubbleGaussian: width sigma of the Gaussian potential at each bubble \[Angstrom\].
* - `Strain.BubbleGaussPotStrength`
  - Strain.BubbleGaussian: amplitude of the Gaussian potential at each bubble \[units of g0\].
```

## Calls

{ref}`AtomsSetCart <api-atoms-atomssetcart>`, {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`RandSeedFromInput <api-random-randseedfrominput>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`)

## Called by

{ref}`HamHopping <api-ham-hamhopping>`
