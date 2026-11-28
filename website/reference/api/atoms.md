---
edit_url: null
---

(api-atoms)=
# `atoms`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/atoms.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/atoms.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`cell <api-cell>`, {ref}`constants <api-constants>`, {ref}`hybrid <api-hybrid>`, {ref}`math <api-math>`, {ref}`mio <api-mio>`, {ref}`name <api-name>`, {ref}`parallel <api-parallel>`
- **Used by:** {ref}`calc <api-calc>`, {ref}`diag <api-diag>`, {ref}`gauss <api-gauss>`, {ref}`grabnes <api-grabnes>`, {ref}`ham <api-ham>`, {ref}`interface <api-interface>`, {ref}`kubo <api-kubo>`, {ref}`kubosubs <api-kubosubs>`, {ref}`neigh <api-neigh>`, {ref}`scf <api-scf>`
- **Public interface:** 4 procedures, 33 variables
- **Input keys:** 119 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`AtomsPos <api-atoms-atomspos>`
  - subroutine
  - —
* - {ref}`AtomsRotate <api-atoms-atomsrotate>`
  - subroutine
  - counterclockwise if angle is positive
* - {ref}`AtomsSetCart <api-atoms-atomssetcart>`
  - subroutine
  - —
* - {ref}`AtomsSetFrac <api-atoms-atomssetfrac>`
  - subroutine
  - —
```

## Module variables

```{list-table}
:header-rows: 1

* - Name
  - Declaration
  - Initial value
  - Description
* - `dIntLay`
  - `real(dp)`
  - —
  - —
* - `displacements`
  - `real(dp), dimension(:,:), pointer`
  - `NULL()`
  - —
* - `displacements_b`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `displacements_t`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `frac`
  - `logical, save`
  - `.true.`
  - —
* - `in1`
  - `integer, save`
  - —
  - For species: 1 - Carbon site A 2 - Carbon site B 3 - Boron 4 - Nitrogen
* - `in2`
  - `integer, save`
  - —
  - For species: 1 - Carbon site A 2 - Carbon site B 3 - Boron 4 - Nitrogen
* - `indxDiv`
  - `integer, dimension(:), pointer`
  - —
  - —
* - `indxNode`
  - `integer, dimension(:), pointer`
  - —
  - —
* - `inode1`
  - `integer, save`
  - —
  - For species: 1 - Carbon site A 2 - Carbon site B 3 - Boron 4 - Nitrogen
* - `inode2`
  - `integer, save`
  - —
  - For species: 1 - Carbon site A 2 - Carbon site B 3 - Boron 4 - Nitrogen
* - `interlayerDistances`
  - `real(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `layerIndex`
  - `integer, dimension(:), pointer`
  - —
  - —
* - `nAt`
  - `integer, save`
  - —
  - —
* - `nAtBN`
  - `integer, save`
  - —
  - —
* - `nAtC`
  - `integer, save`
  - —
  - —
* - `nAtC1`
  - `integer, save`
  - —
  - —
* - `nAtC2`
  - `integer, save`
  - —
  - —
* - `nAtC3`
  - `integer, save`
  - —
  - —
* - `nEl`
  - `real(dp)`
  - —
  - —
* - `phiForEffectiveModel`
  - `real(dp)`
  - —
  - —
* - `procAt`
  - `integer, dimension(:), pointer`
  - —
  - —
* - `Rat`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `RatBN`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `RatC`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `RatC1`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `RatC2`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `RatC3`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `RatInit`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `RatTemp`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `secAt`
  - `integer, dimension(:), pointer`
  - —
  - —
* - `Species`
  - `integer, dimension(:), pointer`
  - —
  - —
* - `SpeciesTemp`
  - `integer, dimension(:), pointer`
  - —
  - —
```

## Input keys read by private procedures

Procedures that are not part of the public interface but read input keys.

**`AtomsOrder`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/atoms.F90#L1890))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Structure.AtomsOrderDeactivated`
  - Keep the atoms in input order (default). .false. reorders them in spatial bins for the parallel decomposition; refused for ReadXYZ, whose per-atom files are in file order.
```
