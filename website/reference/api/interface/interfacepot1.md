---
edit_url: null
---

(api-interface-interfacepot1)=
# `InterfacePot1`

```{index} single: InterfacePot1 (subroutine in interface)
```

```fortran
subroutine InterfacePot1()
```

- **Module:** {ref}`interface <api-interface>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/interface.F90, line 25`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/interface.F90#L25)

## Input keys read (13)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Interlayer.InterfaceClambda`
  - Edge potential A exp(-r/lambda)/r of a graphene/hBN interface: decay length lambda on carbon \[Angstrom\]. InterfacePot1, which reads it, is not called at present.
* - `Interlayer.InterfaceBNlambda`
  - Edge potential: decay length lambda on B and N \[Angstrom\] (InterfacePot1, not called).
* - `Interlayer.InterfaceClambda_B`
  - Edge potential: lambda on carbon from a B edge \[Angstrom\] (InterfacePot1, not called).
* - `Interlayer.InterfaceClambda_N`
  - Edge potential: lambda on carbon from an N edge \[Angstrom\] (InterfacePot1, not called).
* - `Interlayer.InterfaceBNlambda_B`
  - Edge potential: lambda on BN from a B edge \[Angstrom\] (InterfacePot1, not called).
* - `Interlayer.InterfaceBNlambda_N`
  - Edge potential: lambda on BN from an N edge \[Angstrom\] (InterfacePot1, not called).
* - `Interlayer.InterfaceAmpN`
  - Edge potential: amplitude A of N edges \[eV Angstrom\] (InterfacePot1, not called).
* - `Interlayer.InterfaceAmpB`
  - Edge potential: amplitude A of B edges \[eV Angstrom\] (InterfacePot1, not called).
* - `Interlayer.InterfaceAmpCN`
  - Edge potential: amplitude of N edges on carbon \[eV Angstrom\] (InterfacePot1, not called).
* - `Interlayer.InterfaceAmpCB`
  - Edge potential: amplitude of B edges on carbon \[eV Angstrom\] (InterfacePot1, not called).
* - `Interface.EdgeHopping`
  - Edge hopping correction at a graphene/hBN interface (InterfacePot1, not called at present).
* - `Interface.EdgeHoppingAmp`
  - Interface.EdgeHopping: amplitude (InterfacePot1, not called at present).
* - `Interface.EdgeHoppingDamp`
  - Interface.EdgeHopping: damping (InterfacePot1, not called at present).
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_Subsection <api-sys-syssubsection>` (= `SysSubsection`)
