---
edit_url: null
---

(api-tbpar-tbinit)=
# `TBInit`

```{index} single: TBInit (subroutine in tbpar)
```

```fortran
subroutine TBInit()
```

- **Module:** {ref}`tbpar <api-tbpar>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/tbpar.f90, line 20`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/tbpar.f90#L20)

## Input keys read (20)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `TB.NeighLevels`
  - Number of neighbour shells of the hopping table (at least 1).
* - `TB.Hopping`
  - Nearest-neighbour hopping g0 \[eV\], the internal energy unit: every energy is divided by it. Default -3.72 aG + 12.14 (2.9887 eV at aG 2.46).
* - `TB.BNHopping`
  - Nearest-neighbour B-N hopping of hBN \[eV\] (default: TB.Hopping).
* - `TB.CBHopping`
  - Nearest-neighbour C-B hopping \[eV\].
* - `TB.CNHopping`
  - Nearest-neighbour C-N hopping \[eV\].
* - `TB.OnSiteC`
  - On-site energy of carbon \[eV\].
* - `TB.OnSiteC1`
  - On-site energy of carbon sublattice 1 (Species 1) \[eV\].
* - `TB.OnSiteC2`
  - On-site energy of carbon sublattice 2 (Species 2) \[eV\].
* - `TB.OnSiteC1_LB`
  - On-site energy of carbon sublattice 1 of the bottom layer \[eV\], for the stack models that use it.
* - `TB.OnSiteC2_LB`
  - On-site energy of carbon sublattice 2 of the bottom layer \[eV\], for the stack models that use it.
* - `TB.OnSiteC1_LT`
  - On-site energy of carbon sublattice 1 of the top layer \[eV\], for the stack models that use it.
* - `TB.OnSiteC2_LT`
  - On-site energy of carbon sublattice 2 of the top layer \[eV\], for the stack models that use it.
* - `TB.OnSiteB`
  - On-site energy of boron \[eV\].
* - `TB.OnSiteN`
  - On-site energy of nitrogen \[eV\].
* - `TB.InterLayerHopping`
  - Interlayer hopping t exp(-d/dIntLay) of the simplest interlayer branch of HamHopping, d = in-plane distance \[eV\].
* - `Haldane.T2`
  - Here: read to leave the second-neighbour table entries at zero; the Haldane term itself is added in ham.F90.
* - `Run.TypeOfSystem`
  - Here: Ribbons, Hybrid and ReadXYZ structures get the hopping table of mixed C/B/N systems.
* - `Hubbard.U_C`
  - Hubbard U of carbon \[eV\] (two-spin runs with Hubbard.EnableSCF).
* - `Hubbard.U_B`
  - Hubbard U of boron \[eV\] (two-spin runs with Hubbard.EnableSCF).
* - `Hubbard.U_N`
  - Hubbard U of nitrogen \[eV\] (two-spin runs with Hubbard.EnableSCF).
```

## Calls

{ref}`MIO_InputBlock <api-input-interface-inputblock>` (= `InputBlock`), {ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_InputSearchLabel <api-input-inputsearchlabel>` (= `InputSearchLabel`), {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_StringComp <api-string-interface-stringcomp>` (= `StringComp`), {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Section <api-sys-syssection>` (= `SysSection`)

## Called by

{ref}`CalcSelect <api-calc-calcselect>`
