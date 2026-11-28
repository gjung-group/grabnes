---
edit_url: null
---

(api-qsort-quicksort)=
# `QuickSort`

```{index} single: QuickSort (subroutine in qsort)
```

```fortran
recursive subroutine QuickSort(array)
```

- **Module:** {ref}`qsort <api-qsort>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/math/qsort.F90, line 35`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/math/qsort.F90#L35)

Sort an array by numerical ascending order Reference:

\[1\] http://www.fortran.com/qsort\_c.f95

## Parameters

`array` : `real(dp), dimension(:), intent(inout)`
: Array to be sorted
