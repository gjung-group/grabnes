---
edit_url: null
---

(api-interface-interfacepot2)=
# `InterfacePot2`

```{index} single: InterfacePot2 (subroutine in interface)
```

```fortran
subroutine InterfacePot2(H, Nneigh, NList, neighCell, neighD, Species, Nradii)
```

- **Module:** {ref}`interface <api-interface>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/interface.F90, line 120`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/interface.F90#L120)

## Parameters

`H` : `real(dp), dimension(inode1:), intent(inout)`
: —

`Nneigh` : `integer, dimension(inode1:), intent(in)`
: —

`NList` : `integer, dimension(1:maxNeigh,inode1:inode2), intent(in)`
: —

`neighCell` : `integer, dimension(1:3,1:maxNeigh,inode1:inode2), intent(in)`
: —

`neighD` : `real(dp), dimension(:,:,:), intent(in)`
: —

`Species` : `integer, dimension(inode1:), intent(in)`
: —

`Nradii` : `real(dp), dimension(:,:), intent(in)`
: —

## Calls

{ref}`AtomsSetCart <api-atoms-atomssetcart>`, {ref}`MIO_Debug <api-mio-interface-mio_debug>`, {ref}`num2str <api-mio-interface-num2str>`, {ref}`MIO_Kill <api-sys-syskill>` (= `SysKill`), {ref}`MIO_Print <api-sys-interface-sysprint>` (= `SysPrint`)
