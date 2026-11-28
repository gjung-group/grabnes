---
edit_url: null
---

(api-sparse_matrix_module-sparse_matvec)=
# `sparse_matvec`

```{index} single: sparse_matvec (subroutine in sparse_matrix_module)
```

```fortran
subroutine sparse_matvec(nn, row_ptr, col_ind, values, x, y)
```

- **Module:** {ref}`sparse_matrix_module <api-sparse_matrix_module>`
- **Kind:** Subroutine
- **Source:** [`lanczosKuboCode/Src/test_sparse_diag.F90, line 4`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/test_sparse_diag.F90#L4)

## Parameters

`nn` : `integer, intent(in)`
: —

`row_ptr` : `integer, dimension(:), intent(in)`
: —

`col_ind` : `integer, dimension(:), intent(in)`
: —

`values` : `complex(16), dimension(:), intent(in)`
: —

`x` : `complex(16), dimension(:), intent(in)`
: —

`y` : `complex(16), dimension(nn), intent(out)`
: —

## Called by

[DiagHamSparse](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/diag.F90#L5180) ({ref}`diag <api-diag>`), {ref}`test_csr <api-test_csr>` (program)
