# TAPW + SOC Input Parameters

> **Experimental interface.** TAPW with SOC is not yet covered by a public
> end-to-end reference example. Confirm parameter names against
> `lanczosKuboCode/Src/calc.F90` and `ham.F90`, record the commit used, and
> validate the result before using it in production work.

## Required Parameters to Activate TAPW + SOC

### Core Flags
```
useTAPW .true.
RashbaSOCterm .true.  # (or your existing SOC flags)
```

### Essential TAPW Parameters

The following parameters have default values but can be customized:

```
# TAPW G-vector cutoff (default: 180)
# Controls the size of the TAPW basis: M = NG × Nlabel
# Larger values = more accurate but slower
Diag.N_G 5

# Moiré/twist angle in degrees (default: -33.004491598883078)
# For twisted bilayer graphene, set to your actual twist angle
# Note: Can be negative (default is negative)
Diag.MoireAngle 1.08

# Physical twist angle for unitary check (default: 1.08)
# Used if checkTAPWUnitary is enabled
Diag.PhysicalTwistAngle 1.08

# Graphene lattice constant in Angstroms (default: 2.46019)
TAPW.aG 2.46019
```

### Optional TAPW Parameters

```
# Use dense matrix instead of sparse (default: .false.)
# Recommended: .false. for large systems, .true. for small/medium
useDenseMatrixTAPW .false.

# Enable TAPW debug output (default: .false.)
Diag.TAPWDebug .true.

# Use triangular G-vector truncation (default: .false.)
Diag.TriangularTruncation .false.

# Check TAPW unitary transformation (default: .false.)
# If true, automatically calculates optimal NGrange for complete BZ
Diag.CheckTAPWUnitary .false.

# Use K' valley instead of K valley (default: .false.)
Diag.UseKprimeValley .false.

# Use rigid reference positions (default: .false.)
Diag.UseRigidPositions .false.

# G-grid rotation angle in degrees (default: 0.0)
Diag.gGridRotationAngle 0.0
```

### SOC Parameters (your existing ones)

```
# Your existing SOC flags
RashbaSOCterm .true.
IntrinsicSOCterm .false.  # if applicable
PIASOCterm .false.        # if applicable

# Layer-specific SOC control (if applicable)
SOCLayerControl .true.
SOCLayers 1 5             # space-separated layer indices

# SOC coupling strengths
lambdaR 0.1               # Rashba coupling
lambdaI 0.0               # Intrinsic SOC
lambdaPIA 0.0             # PIA SOC
```

## Minimal Working Example

For a basic TAPW + SOC calculation:

```
useTAPW .true.
RashbaSOCterm .true.
Diag.N_G 5
Diag.MoireAngle 1.08
TAPW.aG 2.46019
```

## Notes

1. **TAPW basis size**: The TAPW basis size is M = NG × Nlabel, where Nlabel is the number of unique layer×sublattice combinations. Start with small NG (e.g., 5) and increase if needed.

2. **Performance**: 
   - Small NG (~5-10): Fast, good for testing
   - Medium NG (~20-50): Balanced accuracy/speed
   - Large NG (~100-180): Very accurate but slow

3. **SOC layer control**: If you're using `SOCLayerControl`, make sure the layers you specify exist in your system.

4. **Chern calculation**: To disable Chern calculation (as requested), simply don't set:
   ```
   Diag.calculateChern .true.  # Don't include this!
   ```

## Example Full Input Block

```
# Enable TAPW + SOC
useTAPW .true.
RashbaSOCterm .true.

# TAPW parameters
Diag.N_G 5
Diag.MoireAngle 1.08
TAPW.aG 2.46019

# Optional: Debug output
Diag.TAPWDebug .true.
Diag.socDebug .true.

# SOC parameters (adjust as needed)
lambdaR 0.1

# Layer control (if applicable)
SOCLayerControl .true.
SOCLayers 1 2

# NOTE: Do NOT set Diag.calculateChern for bands-only calculation
```
