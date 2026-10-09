# Single-Layer Graphene SOC Benchmark Tests

> **Validation plan, not yet an automated test suite.** Historical completion
> marks below describe manual development checks. Promote a case to a supported
> regression only after committing its input, observable, tolerance, and
> reference output under this directory.

## Overview
This document provides a systematic test suite for validating all SOC terms on single-layer graphene, starting from simplest to most complex.

## Test Structure
Each test builds on the previous one, adding one SOC term at a time to isolate effects and verify correct implementation.

---

## Test 0: Baseline (No SOC)

### Input Parameters
```fortran
SpinPolarized .true.
EnableSCF .false.

! All SOC terms disabled
ZeemanTerm .false.
IntrinsicSOCterm .false.
RashbaSOCterm .false.
PIASOCterm .false.
```

### Expected Results
- ✅ Two sets of bands (spin-up and spin-down)
- ✅ Dirac point at K point: bands cross at E = 0
- ✅ No gap at Dirac point
- ✅ Identical spin-up and spin-down bands (degenerate)

### Observable
- Energy at Dirac point should be exactly 0 (no gap)
- Spin channels completely degenerate

---

## Test 1: Valley Zeeman Only (λVZ)

### Input Parameters
```fortran
SpinPolarized .true.
EnableSCF .false.

ZeemanTerm .true.
ZeemanFactor 0.002  ! Adjusted to avoid E=0 alignment in Test 3
BmagZeeman 10.0

IntrinsicSOCterm .false.
RashbaSOCterm .false.
PIASOCterm .false.
```

### Expected Results
- ✅ Spin-up bands shifted UP by `+gZeeman × BmagZeeman = +0.01`
- ✅ Spin-down bands shifted DOWN by `-gZeeman × BmagZeeman = -0.01`
- ✅ Bands still CROSS at Dirac point (no gap opened)
- ✅ Total splitting: `2 × 0.001 × 10 = 0.02`

### Observable
- **Energy difference between spin channels at any k-point**
  - Should be exactly `0.02` (in units consistent with your code)
- Spin-up band always above corresponding spin-down band
- **At Dirac point**: Spin-up crosses at E ≈ +0.01, Spin-down crosses at E ≈ -0.01
- Bands still touch/cross (gapless), just at different energy levels

### Verification
```python
# Check energy difference between channels
E_up[i] - E_dn[i] ≈ 0.02  # Should be constant across all i
```

---

## Test 2: Intrinsic SOC Only (λI)

### Input Parameters
```fortran
SpinPolarized .true.
EnableSCF .false.

ZeemanTerm .false.

IntrinsicSOCterm .true.
LambdaI 0.01

RashbaSOCterm .false.
PIASOCterm .false.
```

### Expected Results
- ✅ Gap opens at Dirac point
- ✅ Gap size ≈ `2λI = 2 × 0.01 = 0.02`
- ✅ Both spin channels have identical gap
- ✅ No net spin splitting (both channels shifted equally)

### Observable
- **Energy gap at Dirac point**
  - Valence band maximum: `-λI = -0.01`
  - Conduction band minimum: `+λI = +0.01`
  - Gap: `2λI = 0.02`

### Verification
```python
# Find band indices at Dirac point
# Check that gap = 2*lambdaI
gap = min(E_cond) - max(E_valence)
assert abs(gap - 2*lambdaI) < 0.001
```

---

## Test 3: Intrinsic + Valley Zeeman (λI + λVZ)

### Input Parameters
```fortran
SpinPolarized .true.
EnableSCF .false.

ZeemanTerm .true.
ZeemanFactor 0.002  ! Adjusted to avoid E=0 alignment in Test 3
BmagZeeman 10.0

IntrinsicSOCterm .true.
LambdaI 0.01

RashbaSOCterm .false.
PIASOCterm .false.
```

### Expected Results
- ✅ Gap opens at Dirac point for both spin channels
- ✅ Gap size: `2λI = 0.02` (same for both spin channels)
- ✅ Spin-up and spin-down gaps shifted relative to each other by `2gZeeman×Bmag = 0.02`
- ✅ Spin-resolved gap at Dirac point (different for each spin channel)

### Observable
- **Gap position in energy**
  - Spin-up valence: `-λI + shift = -0.01 + 0.01 = 0.00`
  - Spin-up conduction: `+λI + shift = +0.01 + 0.01 = 0.02`
  - Spin-down valence: `-λI - shift = -0.01 - 0.01 = -0.02`
  - Spin-down conduction: `+λI - shift = +0.01 - 0.01 = 0.00`
- **Key observation**: Spin-up valence band and spin-down conduction band both at E=0 (crossing still possible at Dirac point)
- **Gap size**: Still `2×λI = 0.02` for each spin channel

### Verification
```python
# Both bands should be gapped
# Gap magnitude different for each spin channel
gap_up = min(E_cond_up) - max(E_val_up)  # should be 0.02
gap_dn = min(E_cond_dn) - max(E_val_dn)  # should be 0.02
# But the energy positions are shifted relative to each other
```

---

## Test 4: Rashba Only (λR) - Most Complex

### Input Parameters
```fortran
SpinPolarized .true.
EnableSCF .false.

ZeemanTerm .false.
IntrinsicSOCterm .false.

RashbaSOCterm .true.
LambdaR 0.5  ! Large value for visibility (will be divided by g0)

PIASOCterm .false.

! REQUIRED for this test: Add Intrinsic SOC to open gap
! Without this, Rashba effects are subtle on gapless Dirac cone
IntrinsicSOCterm .true.
LambdaI 0.01
```

### Expected Results
- ✅ **Block Hamiltonian used** (you should see different eigenvalue output)
- ✅ Spin mixing at Dirac point
- ✅ Linear Rashba splitting away from Dirac point
- ✅ Splitting magnitude proportional to λR × k⊥
- ✅ Both spin channels show Rashba contribution

### Observable
- **Energy splitting at finite k**
  - For small k near Dirac point: `ΔE ≈ ±λR × k⊥`
  - Rashba creates characteristic "Rashba cones" in dispersion
  - Maximum splitting occurs perpendicular to specific directions

### Verification
```python
# For band structure at k away from Dirac point
# Should see Rashba-induced splitting that depends on k⊥
# This is the most complex test - requires fitting to Rashba dispersion
```

---

## Test 5: PIA Only (λPIA) - Hopping Terms Critical

### Input Parameters
```fortran
SpinPolarized .true.
EnableSCF .false.

ZeemanTerm .false.
IntrinsicSOCterm .false.
RashbaSOCterm .false.

PIASOCterm .true.
LambdaPIA 0.01
```

### Expected Results
- ✅ Onsite: adds λPIA to diagonal (spin-independent shift)
- ✅ Hopping: creates imaginary hopping amplitude
- ✅ Breaks sublattice/chiral symmetry
- ✅ Weak dispersion modification near Dirac point

### Observable
- **Onsite contribution only**: Simple energy shift, both spins shift equally
- **With hopping**: Subtle asymmetry in dispersion, can create small gap-like features
- **Test both** by comparing onsite-only vs full PIA

### Verification
Compare results with PIA hopping disabled vs enabled to see the difference.

---

## Test 6: All Terms Combined

### Input Parameters
```fortran
SpinPolarized .true.
EnableSCF .false.

ZeemanTerm .true.
ZeemanFactor 0.002  ! Adjusted to avoid E=0 alignment in Test 3
BmagZeeman 10.0

IntrinsicSOCterm .true.
LambdaI 0.01

RashbaSOCterm .true.
LambdaR 0.1  ! Larger value needed for visibility since it's divided by g0

PIASOCterm .true.
LambdaPIA 0.01
```

### Expected Results
- ✅ All effects present simultaneously
- ✅ Gap from Intrinsic SOC
- ✅ Spin splitting from Zeeman
- ✅ Rashba mixing + splitting
- ✅ PIA symmetry breaking

### Observable
- **Complex band structure** with all contributions
- Validate that all terms combine constructively
- No unphysical behavior or numerical instabilities

---

## Test 7: PIA Onsite vs Full (Critical Test)

### Purpose
Verify PIA hopping terms make a difference

### Test A: Onsite Only
```fortran
PIASOCterm .true.
LambdaPIA 0.01
! Comment out PIA hopping in code temporarily
```

### Test B: Full PIA (Onsite + Hopping)
```fortran
PIASOCterm .true.
LambdaPIA 0.01
! Enable PIA hopping in code
```

### Expected Difference
- Onsite only: Just adds constant λPIA to diagonal
- Full PIA: Adds λPIA to diagonal PLUS imaginary hopping terms
- The hopping terms create additional band structure modifications

### Verification
Compare band structures - they should differ near Dirac point if hopping is significant.

---

## Test 8: Rashba Scaling Test

### Purpose
Verify Rashba effect scales linearly with λR

### Tests
```fortran
RashbaSOCterm .true.
LambdaR 0.001  ! Test 1
LambdaR 0.002  ! Test 2
LambdaR 0.005  ! Test 3
```

### Expected
- Splitting magnitude should scale linearly with λR
- `Splitting_2 / Splitting_1 ≈ 2.0`
- `Splitting_3 / Splitting_1 ≈ 5.0`

---

## Test 9: Intrinsic SOC Gap Size

### Purpose
Verify gap scales as 2λI

### Tests
```fortran
IntrinsicSOCterm .true.
LambdaI 0.001  ! Test 1: Gap = 0.002
LambdaI 0.002  ! Test 2: Gap = 0.004
LambdaI 0.005  ! Test 3: Gap = 0.010
```

### Expected
- Gap size = 2 × λI exactly
- Linear scaling verified
- Both spin channels have identical gaps (unless Zeeman also present)

---

## Recommended Testing Order

1. **Test 0**: Baseline - verify no SOC works ✅ **COMPLETED**
2. **Test 1**: Valley Zeeman - simplest, cleanest test ✅ **COMPLETED**
3. **Test 2**: Intrinsic SOC - verify gap opening ✅ **COMPLETED**
4. **Test 3**: Combined (test 1 + 2) - verify terms add properly ✅ **COMPLETED**
5. **Test 9**: Intrinsic scaling - verify λI scaling ✅ **COMPLETED**
6. **Test 4**: Rashba - most complex ❌ **BLOCKED** - Needs Rashba routing fix in all code paths
7. **Test 8**: Rashba scaling - verify λR scaling
8. **Test 7**: PIA comparison - verify hopping matters
9. **Test 5**: PIA standalone
10. **Test 6**: All terms combined - final validation

## Testing Progress (Updated)

### Completed Tests ✅
- **Test 0**: Baseline verification - Spin degeneracy confirmed
- **Test 1**: Valley Zeeman (λVZ) - Spin splitting confirmed, ZeemanFactor adjusted to 0.002
- **Test 2**: Intrinsic SOC (λI) - Gap opening confirmed after sublattice-dependent fix
- **Test 3**: Combined Intrinsic + Zeeman - Both effects work together properly
- **Test 9**: Intrinsic SOC scaling - Linear scaling with λI verified

### Issues Found 🔧
- **Intrinsic SOC**: Required sublattice-dependent implementation (fixed)
- **Rashba routing**: Only implemented in sparseDiagSolver branch, needs expansion to all code paths
- **SOC parameter normalization**: Added g0 normalization for consistency

### Remaining Work 📋
- **Test 4**: Rashba implementation - Fix routing to all Hamiltonian builder branches
- **Test 5**: PIA standalone - After Rashba fix
- **Test 6**: All SOC terms combined - Final integration
- **Test 7**: PIA onsite vs full comparison
- **Test 8**: Rashba scaling verification

---

## Critical Checks for Each Test

### Compilation Check
- ✅ No compilation errors
- ✅ No runtime errors

### Output Structure
- ✅ Proper number of eigenvalues
- ✅ Spin channels handled correctly
- ✅ k-path output formatted correctly

### Physics Check
- ✅ Energy conservation (no unphysical energies)
- ✅ Hermiticity (eigenvalues are real)
- ✅ Symmetry properties match expectations

### Comparison Check
- ✅ Changing λ values produces expected proportional changes
- ✅ Disabling terms reverts to simpler cases
- ✅ No unexplained discontinuities or artifacts

---

## Single-Layer Graphene Input Setup

For clean tests, use single-layer graphene:

```
! Single layer graphene unit cell
nAt = 2 (carbon atoms)
! No interlayer distance
! Simple hexagonal cell with pristine graphene
```

This eliminates moiré effects and provides clean baseline for SOC testing.
