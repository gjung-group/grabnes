---
edit_url: null
---

(api-kubosubs-kubochebyshevevolve)=
# `KuboChebyshevEvolve`

```{index} single: KuboChebyshevEvolve (subroutine in kubosubs)
```

```fortran
subroutine KuboChebyshevEvolve(Z, Psin, Psinm1, c, H0, hopp, NList, ac, bc, nPol)
```

- **Module:** {ref}`kubosubs <api-kubosubs>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/kubosubs.F90, line 557`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/kubosubs.F90#L557)

One Chebyshev time step of a single state: Z \<- U(dT) Z, with the coefficients c of KuboCn. Psin and Psinm1 are work vectors (overwritten).

## Parameters

`Z` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`Psin` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`Psinm1` : `complex(dp), dimension(inode1:), intent(inout)`
: —

`c` : `complex(dp), dimension(nPol), intent(in)`
: —

`H0` : `real(dp), dimension(inode1:), intent(in)`
: —

`hopp` : `complex(dp), dimension(:,inode1:), intent(in)`
: —

`NList` : `integer, dimension(1:,inode1:), intent(in)`
: —

`ac` : `real(dp), intent(in)`
: —

`bc` : `real(dp), intent(in)`
: —

`nPol` : `integer, intent(in)`
: —

## Calls

{ref}`MIO_Debug <api-mio-interface-mio_debug>`

## Called by

{ref}`KuboEvol <api-kubosubs-kuboevol>`
