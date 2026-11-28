---
edit_url: null
---

(api-gauss-gaussheight)=
# `GaussHeight`

```{index} single: GaussHeight (subroutine in gauss)
```

```fortran
subroutine GaussHeight()
```

- **Module:** {ref}`gauss <api-gauss>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/gauss.F90, line 217`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/gauss.F90#L217)

## Input keys read (3)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Strain.GaussHeightPercentage`
  - RandomStrain (refused at present): percentage of atoms that carry a Gaussian bump.
* - `Strain.GaussHeightRange`
  - RandomStrain (refused at present): width sigma of each bump \[Angstrom\].
* - `Strain.GaussHeightStrength`
  - RandomStrain (refused at present): height amplitude of the bumps (added to z) \[Angstrom\].
```

## Calls

{ref}`AtomsSetCart <api-atoms-atomssetcart>`, {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`RandSeedFromInput <api-random-randseedfrominput>`, {ref}`MIO_StringComp <api-string-interface-stringcomp>` (= `StringComp`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`)

## Called by

{ref}`HamInit <api-ham-haminit>`
