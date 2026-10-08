# Historical build report: Apple Silicon test copy (September 2026)

> **Historical record, not current validation.** This report describes a test
> performed on 2026-09-26 on an Apple Silicon (arm64) MacBook Pro running
> macOS 26.5.2, with a **separately patched copy of the solver** that used to
> live in `grabnes_testrun/lanczosKuboCode/`. That copy has been removed. The
> report is not evidence that the canonical `lanczosKuboCode/` sources were
> tested on macOS or on Linux.
>
> Current build instructions, results, and limitations are in
> [`docs/development/cluster-build-and-validation.md`](../docs/development/cluster-build-and-validation.md);
> the harness is described in [`README.md`](README.md).

## What was tested then

- Native arm64 compilation of the test copy with Homebrew GCC 16.2,
  Open MPI 5.0.11, ARPACK 3.9.1_1, and OpenBLAS 0.3.34, using
  `-O0 -g -fcheck=all -fbacktrace` and no OpenMP directives. The configuration
  is preserved as [`config/macos-homebrew.make.sys`](config/macos-homebrew.make.sys).
- `examples/01_graphene_bands`: the program reported `0 errors, 0 warnings`
  and `generate.bands` was byte-for-byte identical to `reference/bands.dat`
  (SHA-256 `312203821daaa76b6f9891988b3a4d59f78e1e7cfc0dc0011a262d394e8c9df8`).

## What happened to the changes made in that copy

| Change in the test copy | Status in the canonical solver |
| --- | --- |
| `.eq.` between logicals replaced by `.eqv.` (`calc.F90`, `diag.F90`) | Applied |
| Unreachable self-aliasing `move_alloc` block removed (`diag.F90`) | Applied |
| `size()` of the unassociated MPI receive list guarded (`ham.F90`) | Applied, extended to `calc.F90`, and the list pointers are now initialized to `NULL()` |
| Array size taken before deallocation in the memory accounting (`MIO/mem_temp.f90`) | Applied |
| 1024-character file-name buffers in the two code generators | Applied |
| `-fallow-argument-mismatch`, `-ffree-line-length-none` | In `make.sys.example` |
| `-std=legacy` | Not needed; not used |
| `-DTIMER` given explicitly | No longer needed: the Makefile adds it as it was always meant to |
| OpenMP directives disabled, `libgomp` linked | Kept for GNU Fortran, see the validation document |
| `nonBulkSmall .true.` and `WriteDataFiles .false.` in the graphene input | Already part of `examples/01_graphene_bands/Gendata.in` |
| `WritePos .true.` in the test copy of the graphene input | Not kept; the canonical input uses `.false.` |

The test copy only ever ran the two-atom graphene example. The additional
defects found when the canonical solver was built with optimization and run on
all four examples on Linux are listed in the validation document.
