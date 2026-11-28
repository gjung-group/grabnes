---
edit_url: null
---

(api-neigh)=
# `neigh`

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/neigh.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L1)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`atoms <api-atoms>`, {ref}`cell <api-cell>`, {ref}`constants <api-constants>`, {ref}`math <api-math>`, {ref}`mio <api-mio>`, {ref}`parallel <api-parallel>`, {ref}`tbpar <api-tbpar>`
- **Used by:** {ref}`calc <api-calc>`, {ref}`diag <api-diag>`, {ref}`gauss <api-gauss>`, {ref}`ham <api-ham>`, {ref}`interface <api-interface>`, {ref}`kubosubs <api-kubosubs>`, {ref}`scf <api-scf>`
- **Public interface:** 9 procedures, 18 variables
- **Input keys:** 16 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`fastNN <api-neigh-fastnn>`
  - subroutine
  - calculate nearest neighbors of a sample with a fast method that uses binning
* - {ref}`fastNNnotsquare <api-neigh-fastnnnotsquare>`
  - subroutine
  - —
* - {ref}`fastNNnotsquareBulk <api-neigh-fastnnnotsquarebulk>`
  - subroutine
  - —
* - {ref}`fastNNnotsquareBulkSmall <api-neigh-fastnnnotsquarebulksmall>`
  - subroutine
  - —
* - {ref}`fastNNnotsquareNotRectangle <api-neigh-fastnnnotsquarenotrectangle>`
  - subroutine
  - —
* - {ref}`fastNNnotsquareSmall <api-neigh-fastnnnotsquaresmall>`
  - subroutine
  - —
* - {ref}`NeighList <api-neigh-neighlist>`
  - subroutine
  - —
* - {ref}`NeighShellCutoff2 <api-neigh-neighshellcutoff2>`
  - function
  - Squared in-plane search radius that contains the first nshell neighbour shells of a honeycomb lattice with bond length aCC.
* - {ref}`NeighShellWarn2 <api-neigh-neighshellwarn2>`
  - function
  - Squared in-plane length beyond which an intralayer pair is closer to the search radius of NeighShellCutoff2 than a rigid lattice allows: the last quarter of the interval between the last requested shell and the radius.
```

## Module variables

```{list-table}
:header-rows: 1

* - Name
  - Declaration
  - Initial value
  - Description
* - `maxNeigh`
  - `integer`
  - `3`
  - —
* - `maxShells`
  - `integer, parameter`
  - `8`
  - Neighbour shells of the honeycomb lattice: radius in units of the bond length and number of sites. The entry after the last shell only serves to place the outermost search radius.
* - `neighCell`
  - `integer, dimension(:,:,:), pointer`
  - —
  - —
* - `NeighD`
  - `real(dp), dimension(:,:,:), pointer`
  - —
  - —
* - `NeighIntraWarn2`
  - `real(dp), save`
  - `0.0_dp`
  - Squared in-plane length above which an intralayer pair counts as "close to the search radius" in NeighSearchLayered (0 = no check). Set by HamInit.
* - `NList`
  - `integer, dimension(:,:), pointer`
  - —
  - —
* - `NList2`
  - `integer, dimension(:,:), pointer`
  - —
  - —
* - `Nneigh`
  - `integer, dimension(:), pointer`
  - —
  - —
* - `nodeSndRcv`
  - `integer, dimension(:,:), pointer`
  - `NULL()`
  - —
* - `Nradii`
  - `real(dp), dimension(:,:), pointer`
  - —
  - —
* - `rcvIndx`
  - `integer, dimension(:,:), pointer`
  - `NULL()`
  - The halo lists are only allocated for nProc>1; keep them disassociated (rather than undefined) otherwise so that associated() can be tested.
* - `rcvList`
  - `integer, dimension(:), pointer`
  - `NULL()`
  - The halo lists are only allocated for nProc>1; keep them disassociated (rather than undefined) otherwise so that associated() can be tested.
* - `rcvSz`
  - `integer`
  - —
  - —
* - `shellCount`
  - `integer, dimension(maxShells), parameter`
  - `[3, 6, 3, 6, 6, 6, 6, 3]`
  - —
* - `shellRadius`
  - `real(dp), dimension(maxShells+1), parameter`
  - `[1.0_dp, 1.7320508075688772_dp, 2.0_dp, 2.6457513110645906_dp, 3.0_dp, 3.4641016151377544_dp, 3.6055512754639891_dp, 4.0_dp, 4.3588989435406740_dp]`
  - —
* - `sndIndx`
  - `integer, dimension(:,:), pointer`
  - `NULL()`
  - —
* - `sndList`
  - `integer, dimension(:), pointer`
  - `NULL()`
  - —
* - `sndSz`
  - `integer`
  - —
  - —
```

## Input keys read by private procedures

Procedures that are not part of the public interface but read input keys.

**`NeighFinishSearch`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L319))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Output.WriteDataFiles`
  - Write the neighbour lists and the Hamiltonian to files (v, dx, dy, dz, pos, neighCell.dat, neighD.dat, \<prefix>.e, \<prefix>.s.mag, ...): a complete export of H, large for big cells.
```

**`NeighSetRadii`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L348))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Run.TypeOfSystem`
  - Here: BoronNitride uses the hBN lattice constant for every neighbour radius.
```

**`NeighReadDataFiles`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L446))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Neigh.ReadDetails`
  - With Output.ReadDataFiles: also read neighCell.dat and neighD.dat (image cell and bond vector of each pair).
```

**`NeighListOld`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/neigh.F90#L1023))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Neigh.Distance`
  - NeighListOld: first-neighbour distance \[Angstrom\]; pairs within 1.1 times it are neighbours.
```
