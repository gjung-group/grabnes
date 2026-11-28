---
edit_url: null
---

(api-ham)=
# `ham`

Tight-binding Hamiltonian module for twisted bilayer graphene and related 2D materials.

Provides routines to initialize on-site terms and compute hopping integrals for tight-binding Hamiltonians, including moiré strain, Zeeman terms, and interlayer coupling models. Parallelization via OpenMP/MPI is supported.

- **Kind:** Module
- **Source:** [`lanczosKuboCode/Src/ham.F90`](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L27)
- **Group:** {ref}`Solver <api-group-solver>`
- **Uses:** {ref}`atoms <api-atoms>`, {ref}`cell <api-cell>`, {ref}`constants <api-constants>`, {ref}`gauss <api-gauss>`, `ieee_arithmetic`, {ref}`interface <api-interface>`, {ref}`magf <api-magf>`, {ref}`math <api-math>`, {ref}`mio <api-mio>`, {ref}`moireBLShift <api-moireblshift>`, {ref}`name <api-name>`, {ref}`neigh <api-neigh>`, {ref}`parallel <api-parallel>`, {ref}`random <api-random>`, {ref}`tbpar <api-tbpar>`
- **Used by:** {ref}`calc <api-calc>`, {ref}`diag <api-diag>`, {ref}`scf <api-scf>`
- **Public interface:** 8 procedures, 32 variables
- **Input keys:** 449 reads

## Procedures

```{list-table}
:header-rows: 1

* - Procedure
  - Kind
  - Summary
* - {ref}`HaldaneEnabledForLayer <api-ham-haldaneenabledforlayer>`
  - function
  - Check if Haldane should be applied to a specific layer
* - {ref}`HamHopping <api-ham-hamhopping>`
  - subroutine
  - Assemble hopping integrals (intra/inter-layer) into hopp(:,:).
* - {ref}`HamInit <api-ham-haminit>`
  - subroutine
  - Initialize Hamiltonian parameters and neighbor lists.
* - {ref}`HamInputEnergyEV <api-ham-haminputenergyev>`
  - subroutine
  - Reads an energy that the input file gives in eV and returns it in units of g0.
* - {ref}`HamOnSite <api-ham-hamonsite>`
  - subroutine
  - Compute and add on-site energy terms to H0.
* - {ref}`HamPeierlsPeriodicGauge <api-ham-hampeierlsperiodicgauge>`
  - subroutine
  - Peierls phases of the periodic gauge (Cresti) for a field of flux = B pi / flux quantum \[1/Ang^2\], applied to every hopping, with the corner term exp(i pi m n1 n2).
* - {ref}`SOCEnabledForLayer <api-ham-socenabledforlayer>`
  - function
  - —
* - {ref}`ZeemanEnabledForLayer <api-ham-zeemanenabledforlayer>`
  - function
  - Check if SOC should be applied to a specific layer
```

## Module variables

```{list-table}
:header-rows: 1

* - Name
  - Declaration
  - Initial value
  - Description
* - `gZeeman`
  - `real(dp), save`
  - —
  - —
* - `H0`
  - `real(dp), dimension(:), pointer`
  - `NULL()`
  - Hamiltonian matrices (public interface)
* - `H0Bottom`
  - `real(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `H0Top`
  - `real(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `HABimag`
  - `real(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `HABreal`
  - `real(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `HaldaneLayerControl`
  - `logical, save`
  - —
  - Haldane layer control
* - `HaldaneLayersArray`
  - `integer, dimension(:), save, allocatable`
  - —
  - —
* - `Ho`
  - `real(dp), dimension(:), pointer`
  - `NULL()`
  - —
* - `hopp`
  - `complex(dp), dimension(:,:), pointer`
  - `NULL()`
  - —
* - `InterlayerSOCterm`
  - `logical, save`
  - `.false.`
  - —
* - `IntrinsicSOCterm`
  - `logical, save`
  - —
  - SOC terms flags
* - `IsingSOCterm`
  - `logical, save`
  - —
  - —
* - `KaneMeleSOCterm`
  - `logical, save`
  - `.false.`
  - —
* - `lambdaI`
  - `real(dp), save`
  - —
  - SOC parameters
* - `lambdaIL`
  - `real(dp), save`
  - `0.0_dp`
  - —
* - `lambdaIsing`
  - `real(dp), save`
  - —
  - —
* - `lambdaKM`
  - `real(dp), save`
  - `0.0_dp`
  - —
* - `lambdaPIA`
  - `real(dp), save`
  - —
  - —
* - `lambdaPIAA`
  - `real(dp), save`
  - `0.0_dp`
  - SOC.PIA (since November 2026): Gmitra et al. PRB 93, 155104 (2016), 2i/3 lamPIA^X (s x d)\_z on the second-neighbour graphene bonds of sublattice X, in the spinful block Hamiltonian (socBond table)
* - `lambdaPIAB`
  - `real(dp), save`
  - `0.0_dp`
  - SOC.PIA (since November 2026): Gmitra et al. PRB 93, 155104 (2016), 2i/3 lamPIA^X (s x d)\_z on the second-neighbour graphene bonds of sublattice X, in the spinful block Hamiltonian (socBond table)
* - `lambdaR`
  - `real(dp), save`
  - —
  - —
* - `nspin`
  - `integer, save`
  - —
  - —
* - `PIASOCterm`
  - `logical, save`
  - —
  - —
* - `rashbaAngle`
  - `real(dp), save`
  - `0.0_dp`
  - —
* - `RashbaConsistentPhase`
  - `logical, save`
  - `.false.`
  - Spinful bond terms of the block (2N x 2N) Hamiltonian, opt-in (2026-10-16).  All need SOC.RashbaConsistentPhase, which puts every spin-flip and spin-dependent bond term on the Bloch phase of the spin-diagonal hopping. SOC.RashbaConsistentPhase: Rashba in the Kane-Mele form 2i/3 lamR' (s x d)\_z with lamR' = 3/2 SOC.LambdaR (the grabnes normalisation is kept: the K-point splitting is 3 SOC.LambdaR), Bloch phase of the hopping. Off: the former term, whose phase carries the opposite sign of k.d (exact only at K and Gamma). SOC.RashbaAngle \[deg\]: spin rotation exp(-i phi s\_z/2) of the Rashba term; 0 = tangential, 90 = radial texture. SOC.KaneMele + SOC.LambdaKM \[eV\]: i lamKM/(3 sqrt3) nu\_ij s\_z on the second-neighbour graphene bonds, nu\_ij from HamHaldaneChirality; gap 2 lamKM in a single layer.  (IntrinsicSOCterm is a spin-independent sublattice mass.) SOC.Interlayer + SOC.LambdaInterlayer \[eV\]: the interlayer spin-flip term lambda\_4 of Parra-Martinez et al., arXiv:2610.05905, between same-sublattice sites of adjacent UNTWISTED graphene layers at in-plane distance aCC: i lam4/3 c\_ij sgn(z\_j - z\_i) (d\_ij x s)\_z, c\_ij = +1 if d\_ij is parallel to a first-neighbour bond of i, -1 if antiparallel (inversion fixes this sublattice sign); at K its matrix element is lam4 sigma\_z s\_+/- .
* - `RashbaSOCterm`
  - `logical, save`
  - —
  - —
* - `socBond`
  - `complex(dp), dimension(:,:,:,:), save, allocatable`
  - —
  - per bond (j,i), filled by HamSOCBondTables: socBond(s,s',j,i) = \<i s\|H\|NList(j,i) s'> without the Bloch phase
* - `socBondOn`
  - `logical, dimension(:,:), save, allocatable`
  - —
  - —
* - `SOCLayerControl`
  - `logical, save`
  - —
  - SOC layer control
* - `SOCLayersArray`
  - `integer, dimension(:), save, allocatable`
  - —
  - —
* - `Zterm`
  - `logical, save`
  - —
  - —
```

## Input keys read by private procedures

Procedures that are not part of the public interface but read input keys.

**`HamMeasureLayerTwist`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4124))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Stack.MeasureTwist`
  - Report the twist of every interface measured from the structure.
* - `Stack.ModelFromTwist`
  - Untwisted AB pairs found in the structure get the bilayer F2G2 intra- and interlayer model (Jung-MacDonald PRB 89 035405); AA and twisted pairs keep the two-centre model.
* - `Stack.AlignModelEnergies`
  - With Stack.ModelFromTwist: shift the AB-pair sites so their K-point energy matches the single-layer set.
* - `Stack.TwistTolerance`
  - Stack.MeasureTwist / ModelFromTwist: twists below this count as untwisted \[degrees\].
```

**`HamCheckHermiticity`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4491))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Moire.OffDiag`
  - Off-diagonal (bond) moire term of the G/hBN effective model, Cab. Evaluated at atom i: H slightly non-Hermitian unless Moire.OffDiagMidpoint.
* - `TBG.OffDiag`
  - Off-diagonal (bond) stacking-dependent term of twisted bilayer graphene.
* - `GBN.OffDiag`
  - Off-diagonal (bond) stacking-dependent term of G/hBN, from atom i's displacement.
* - `Strain.StrainedMoire`
  - Strained-moire hopping term (checked for Hermiticity here).
* - `Moire.OffDiagMidpoint`
  - Moire.OffDiag: evaluate the bond term at the bond midpoint (one shared routine for both moire slots): Hermitian H, needed for recursion runs. Plain effective model only.
```

**`HamMoireFourierSetup`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L4894))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Structure.SuperCellX`
  - Here: input cell = simulation cell / SuperCell for the Moire.Fourier wave vectors.
* - `Moire.FourierKprime`
  - Moire.Fourier: the term U is given for the K' valley.
* - `Moire.FourierGrapheneAngle`
  - Moire.Fourier: orientation of the graphene lattice \[degrees\].
* - `Moire.FourierLatticeConstant`
  - Moire.Fourier: graphene lattice constant of the Fourier components \[Angstrom\] (default aG).
```

**`HamHopValleyZeeman`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L5061))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Valley.Zeeman`
  - Valley Zeeman term lambda i nu\_ij s\_i /(3 sqrt 3) on the graphene second-neighbour bonds of each layer (= tau\_z at the Dirac points) \[eV\]; 0 = off. Added before the Peierls phase.
```

**`HamSOCBondTables`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L12929))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Stack.TwistTolerance`
  - Here: interfaces twisted by less than this count as aligned for the SOC bond tables \[degrees\].
```

**`HamDisplacementsFromPositions`** (private, [source](https://github.com/gjung-group/grabnes/blob/1e390c6ead7c835389ae89a599cb1f39e3fc691e/lanczosKuboCode/Src/ham.F90#L13132))

```{list-table}
:header-rows: 1

* - Key
  - Description
* - `Structure.DisplacementsFromPositions`
  - Here: the stacking displacements are computed from generateInit.xyz and the structure.
* - `Stack.GlobalTwist`
  - Here: refused with Structure.DisplacementsFromPositions (those displacements are already in the lab frame).
* - `Stack.EncapsulatedThreeLayers`
  - Here: graphene atoms get a bottom and a top displacement.
```
