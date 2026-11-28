---
edit_url: null
---

(api-gauss-gausspot)=
# `GaussPot`

```{index} single: GaussPot (subroutine in gauss)
```

```fortran
subroutine GaussPot(H, strength)
```

- **Module:** {ref}`gauss <api-gauss>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/gauss.F90, line 18`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/gauss.F90#L18)

## Parameters

`H` : `real(dp), dimension(inode1:), intent(inout)`
: —

`strength` : `real(dp), intent(in)`
: —

## Input keys read (2)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Disorder.GaussPotPercentage`
  - Gaussian disorder: percentage of the atoms that carry a Gaussian impurity potential.
* - `Disorder.GaussPotRange`
  - Gaussian disorder: width sigma of each impurity, exp(-r^2/(2 sigma^2)) \[Angstrom\]; random amplitudes in \[-w/2, w/2\].
```

## Calls

{ref}`AtomsSetCart <api-atoms-atomssetcart>`, {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`RandSeedFromInput <api-random-randseedfrominput>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`)

## Called by

{ref}`HamOnSite <api-ham-hamonsite>`
