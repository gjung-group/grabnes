---
edit_url: null
---

(api-cell-cellget)=
# `CellGet`

```{index} single: CellGet (subroutine in cell)
```

```fortran
subroutine CellGet()
```

- **Module:** {ref}`cell <api-cell>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/cell.F90, line 16`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/cell.F90#L16)

## Input keys read (14)

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Structure.LatticeParameter`
  - Graphene lattice constant aG \[Angstrom\]; the cells built by the code are multiples of it.
* - `Structure.LatticeParameterBN`
  - hBN lattice constant \[Angstrom\] (cells built from Structure.MoireCellParameters set it from the indices).
* - `Run.TypeOfSystem`
  - Kind of structure: built by the code (Graphene, BoronNitride, Ribbons, TwistedBilayerBasedOnMoireCell, ...) or read from a file (ReadXYZ). It selects the cell, the atoms and parts of the Hamiltonian.
* - `Structure.CellHeight`
  - Length of the third (out-of-plane) cell vector \[Angstrom\] for the cells built by the code.
* - `Structure.SuperCell`
  - Repeat the cell N x N in the plane (the atoms are copied, no reordering: copy 0 holds the original atoms). Equivalent to an N x N k-sampling of the original cell for local quantities.
* - `Structure.SuperCellAsymmetric`
  - Repeat the cell Structure.SuperCell times along a1 and Structure.SuperCellY times along a2.
* - `Structure.SuperCellY`
  - Structure.SuperCellAsymmetric: number of copies along the second lattice vector.
* - `Structure.BasedOnMoireCellParameters`
  - Graphene / BoronNitride / TwistedBilayer cells: build the commensurate cell from the four indices of Structure.MoireCellParameters instead of CellSize x CellSize graphene cells.
* - `Structure.MoireCellParameters`
  - Four integers m1 m2 m3 m4 of the commensurate moire cell: a1 = m1 a + m2 b on the first lattice, (m3, m4) on the second; the twist angle and lattice ratio follow from them and are printed.
* - `Structure.CellSize`
  - Number of graphene cells along each lattice vector of a cell built by the code (TypeOfSystem Graphene, BoronNitride); the cell is CellSize\*SuperCell graphene cells on a side.
* - `Structure.RibbonType`
  - TypeOfSystem Ribbons: edge of the ribbons, Zigzag or Armchair.
* - `Structure.GrapheneWidth`
  - TypeOfSystem Ribbons: width of the graphene ribbon, in rows of the ribbon cell.
* - `Structure.BNWidth`
  - TypeOfSystem Ribbons: width of the hBN ribbon, in rows of the ribbon cell.
* - `Structure.XYZFile`
  - TypeOfSystem ReadXYZ: structure file. Its first three lines are the lattice vectors \[Angstrom\], the atoms follow.
```

## Calls

{ref}`MIO_InputParameter <api-input-interface-inputparameter>` (= `InputParameter`), {ref}`MIO_InputSearchLabel <api-input-inputsearchlabel>` (= `InputSearchLabel`), {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`MIO_TimerCount <api-mio-interface-mio_timercount>`, {ref}`MIO_TimerStop <api-mio-interface-mio_timerstop>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_StringComp <api-string-interface-stringcomp>` (= `StringComp`), {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`), {ref}`MIO_PrintDebug <api-sys-sysprintdebug>` (= `SysPrintDebug`), {ref}`MIO_Section <api-sys-syssection>` (= `SysSection`)

## Called by

{ref}`grabnes <api-grabnes>` (program)
