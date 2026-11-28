---
edit_url: null
---

(api-hybrid-hybridgen)=
# `HybridGen`

```{index} single: HybridGen (subroutine in hybrid)
```

```fortran
subroutine HybridGen(nAt, Rat, Species, in1, in2)
```

- **Module:** {ref}`hybrid <api-hybrid>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/hybrid.F90, line 13`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/hybrid.F90#L13)

## Parameters

`nAt` : `integer, intent(in)`
: —

`Rat` : `real(dp), dimension(3,nAt), intent(in)`
: —

`Species` : `integer, dimension(nAt), intent(inout)`
: —

`in1` : `integer, intent(in)`
: —

`in2` : `integer, intent(in)`
: —

## Input keys read (5)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Hybrid.Radius`
  - TypeOfSystem Hybrid (refused at present): radius of the hBN islands in graphene \[Angstrom\].
* - `Hybrid.Percentage`
  - TypeOfSystem Hybrid (refused at present): percentage of the area covered by hBN islands.
* - `Structure.SuperCell`
  - Here (TypeOfSystem Hybrid): sets the size of the region where islands are placed.
* - `Structure.CellSize`
  - Here (TypeOfSystem Hybrid): sets the size of the region where islands are placed.
* - `Hybrid.RandomSites`
  - TypeOfSystem Hybrid (refused at present): place the islands at random sites.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)

## Called by

{ref}`AtomsPos <api-atoms-atomspos>`
