---
edit_url: null
---

(api-mem-memcount)=
# `MemCount`

```{index} single: MemCount (subroutine in mem)
```

```fortran
subroutine MemCount(sz, arrname, rname)
```

- **Module:** {ref}`mem <api-mem>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/MIO/mem_in.f90, line 55`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/MIO/mem_in.f90#L55)

Counts the memory in use

## Parameters

`sz` : `integer(8), intent(in)`
: Size in bytes of the array allocated (if sz>0) or deallocated (if sz\<0)

`arrname` : `character(len=*), intent(in)`
: Name of the array

`rname` : `character(len=*), intent(in)`
: Name of the routine that ask for memory
