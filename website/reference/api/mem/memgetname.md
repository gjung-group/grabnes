---
edit_url: null
---

(api-mem-memgetname)=
# `MemGetName`

```{index} single: MemGetName (subroutine in mem)
```

```fortran
subroutine MemGetName(arrname, rname, name, routine)
```

- **Module:** {ref}`mem <api-mem>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/mem_in.f90, line 269`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mem_in.f90#L269)

Get names for the variable and the routine

## Parameters

`arrname` : `character(len=32), intent(out)`
: Final name of the array

`rname` : `character(len=32), intent(out)`
: Final name of the routine

`name` : `character(len=*), intent(in), optional`
: Name of the variable

`routine` : `character(len=*), intent(in), optional`
: Calling routine
