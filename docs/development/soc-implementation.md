# Spin-orbit coupling and TAPW implementation status

> Development status, reviewed against the source on 2026-09-27. This page
> describes research functionality, not a validated public example.

## Scope

GRABNES contains intrinsic, Ising/valley-Zeeman, Rashba, and PIA spin-orbit
coupling (SOC) paths. The implementation follows the general model discussed
by Gmitra *et al.*, Phys. Rev. B **93**, 155104 (2016),
doi:10.1103/PhysRevB.93.155104.

The code is the source of truth. Input parsing and SOC parameters are in
`lanczosKuboCode/Src/ham.F90`; Hamiltonian construction, TAPW projection, and
diagonalization are in `lanczosKuboCode/Src/diag.F90`.

## Current architecture

The scalar-spin path applies diagonal SOC contributions through
`ApplySOCtoHamiltonian`. Rashba coupling requires a combined spin block, which
is assembled by `BuildBlockHamiltonianOnly` and projected through
`DiagH0TAPW_withBlockH`. PIA also has hopping contributions in the maintained
source. Layer selection is handled by `SOCEnabledForLayer` when
`SOCLayerControl` is enabled.

Relevant input switches include:

```text
IntrinsicSOCterm
IsingSOCterm
RashbaSOCterm
PIASOCterm
LambdaI
LambdaIsing
LambdaR
LambdaPIA
SOCLayerControl
SOCLayers
```

TAPW controls are read using the `Diag.` namespace; see the
[TAPW/SOC input guide](../user-guide/tapw-soc-input.md) and confirm values
against `calc.F90` when developing a new workflow.

## Implemented pieces

- Parsing and normalization of the four SOC coupling parameters.
- Layer-specific SOC selection.
- Diagonal intrinsic, Ising, and PIA contributions.
- Rashba spin-flip hopping in a `2N x 2N` block Hamiltonian.
- A block-Hamiltonian TAPW projection path.
- TAPW valley selection through `Diag.UseKprimeValley`.
- A forced block-TAPW diagnostic path through `Diag.forceBlockTAPW`.

Presence in the source does not by itself constitute physical validation.
Only workflows represented by runnable examples and reference data should be
described as supported.

## Validation status

The historical single-layer checks are retained in `tests/soc/README.md` at
the repository root. They record completed and pending manual checks, but they
do not yet provide executable regression fixtures. Before announcing SOC/TAPW
support, add small inputs and reference observables for:

1. no-SOC spin degeneracy;
2. intrinsic and Ising gap/splitting scaling;
3. Rashba spin mixing and coupling scaling;
4. PIA onsite versus hopping behavior;
5. equivalence of ordinary and forced block-TAPW paths when SOC is disabled;
6. K/K-prime sign conventions; and
7. layer-selective SOC.

Each regression should state units, expected tolerance, compiler information,
and the source of the physical reference value.

## Known risks

- The block Hamiltonian uses roughly four times the storage of a scalar-spin
  Hamiltonian before TAPW reduction.
- Several pathways are selected by nested runtime flags, so branch coverage is
  important.
- Debug output and research-only controls remain in the main source.
- The historical notes contained status statements that predated the current
  block-TAPW routines; those statements have been removed from this maintained
  summary.

## Maintenance rule

Update this page when a relevant routine or input name changes. Put design
proposals and unresolved defects in development or known-issue pages; put
stable instructions in the user guide only after an automated example covers
them.
