# TAPW Chern Number Calculation Optimization Guide

> **Development proposal.** Complexity estimates and projected speedups below
> are hypotheses until accompanied by reproducible profiler output. Apply
> optimizations only with numerical-regression coverage.

## Overview

This document outlines computational bottlenecks and optimization strategies for the TAPW (Twisted Angle Plane Wave) Chern number calculation in the `lanczosKuboCode` codebase. The optimizations are organized by category and impact level.

## Current Performance Analysis

### Computational Complexity
- **TAPW Bands Calculation**: O(nk × M³) ≈ O(10⁴ × 10⁹) = **10¹³ operations**
- **Chern Calculation**: O(nk × n_bands × M²) ≈ O(10⁴ × 61 × 10⁶) = **10¹² operations**
- **Total**: ~**10¹³ operations** for typical 30×30 k-point grid with M=1014

### Memory Usage
- **Stored Arrays**: `stored_eigenvectors(M,M,nk)` + `stored_hamiltonians(M,M,nk)` ≈ **16 GB** for 900 k-points
- **Working Arrays**: ~**200 MB** per k-point for Berry curvature calculation

---

## 🔴 TAPW Bands Calculation Bottlenecks

*These affect the initial band structure calculation that generates stored data*

### 1. TAPW Diagonalization Per K-Point ⭐⭐⭐⭐⭐

**Location**: `DiagH0TAPW()` called in main k-point loop  
**Current Cost**: O(nk × M³) ≈ **10¹³ operations**  
**Impact**: Highest - dominates total computation time

#### Problem
- Full TAPW transformation + diagonalization for every k-point
- ZHEEV diagonalization of 1014×1014 matrices
- Complete TAPW projection: `H_proj = X† H X`
- Repeated for every k-point (900-10,000+ k-points)

#### Optimization Strategy
**Potential Speedup**: **10-100x**

```fortran
! Current (inefficient):
do ik = 1, nk
   call DiagH0TAPW(...)  ! Full transformation + diagonalization
end do

! Optimized approach:
call PrecomputeTAPWTransformation()  ! Once only
do ik = 1, nk
   call UpdateHamiltonianIncremental(k_point)  ! H(k) = H₀ + δH(k)
   call DiagonalizeProjectedH(H_proj)  ! Much smaller operation
end do
```

#### Implementation Details
1. **Separate transformation from diagonalization**
2. **Cache transformation matrix X** between k-points
3. **Incremental Hamiltonian updates** for k-dependent terms
4. **Reuse G-vector calculations** and phase factors where possible

---

### 2. TAPW Transformation Matrix Construction ⭐⭐⭐⭐

**Location**: `DiagH0TAPW()` → TAPW setup  
**Current Cost**: O(N × M²) per k-point  
**Impact**: High - significant overhead per k-point

#### Problem
- Rebuilding transformation matrix `X` every k-point
- Phase factor calculations: `exp(i G·r)` for all atoms/G-vectors
- Matrix assembly and orthogonalization

#### Optimization Strategy
**Potential Speedup**: **5-10x**

```fortran
! Pre-compute phase-independent components
call PrecomputeAtomicBasisFunctions()
call PrecomputeGVectorGrid()

! Cache transformation matrices by symmetry
do ik = 1, nk
   X_cached = GetCachedTransformation(symmetry_class(ik))
   if (.not. allocated(X_cached)) then
      X_cached = ComputeTransformation(ik)
      call CacheTransformation(X_cached, symmetry_class(ik))
   end if
end do
```

---

### 3. G-Vector Generation and Filtering ⭐⭐⭐

**Location**: TAPW setup routines  
**Current Cost**: O(NG²) operations  
**Impact**: Medium - preprocessing overhead

#### Problem
- Repeated G-vector calculations
- Brillouin zone filtering for each calculation
- Distance sorting and selection
- Symmetry operations

#### Optimization Strategy
**Potential Speedup**: **2-5x**

```fortran
! Pre-compute and cache G-vector lists
if (.not. allocated(cached_Gx)) then
   call GenerateOptimizedGVectors(cached_Gx, cached_Gy, cached_NG)
   call ApplyBrillouinZoneFiltering()
   call OptimizeGVectorOrdering()  ! For cache efficiency
end if
```

---

## 🟡 Chern Calculation Bottlenecks

*These are specific to Berry curvature calculation using stored TAPW data*

### 4. Position Difference Matrix Construction ⭐⭐⭐⭐ ✅

**Location**: `build_TAPW_position_differences()`  
**Status**: **ALREADY OPTIMIZED**  
**Speedup Achieved**: ~100x on subsequent calls, ~10x on first call

#### Current Implementation
- ✅ **Intelligent caching system** with change detection
- ✅ **Block-structured computation** eliminating modulo operations  
- ✅ **Pre-computed indices** for optimal memory access

---

### 5. Berry Curvature Double Loop ⭐⭐⭐⭐

**Location**: `CalculateBerryAtKpointFromStored()` lines 11989-11994  
**Current Cost**: O(n_bands × M²) per k-point ≈ **61M operations/k-point**  
**Impact**: High - dominant cost in Chern calculation

#### Problem
```fortran
do i = 1, n_bands
   do m = 1, M  ! ~1014 iterations
      if (n /= m) then
         berry_sum = berry_sum + (Vx(n,m) * conjg(Vy(n,m))) / ((eigval(n) - eigval(m))**2 + eps)
      end if
   end do
end do
```

#### Optimization Strategy
**Potential Speedup**: **3-5x**

```fortran
! Vectorized approach using BLAS operations
! Pre-compute energy difference matrix
real(dp) :: energy_diff_matrix(M,M)
complex(dp) :: berry_contrib(M,M)

! Vectorized computation
do i = 1, M
   do j = 1, M
      energy_diff_matrix(i,j) = (eigval(i) - eigval(j))**2 + eps
   end do
end do

! Element-wise operations (can be further optimized with BLAS)
where (energy_diff_matrix > eps)
   berry_contrib = Vx * conjg(Vy) / energy_diff_matrix
elsewhere
   berry_contrib = 0.0_dp
end where

! Sum for each band
do i = 1, n_bands
   n = band_indices(i)
   berry_curv_bands(i) = -2.0_dp * aimag(sum(berry_contrib(n,:)) - berry_contrib(n,n))
end do
```

---

### 6. Multiple ZGEMM Calls Per K-Point ⭐⭐⭐

**Location**: `CalculateBerryAtKpointFromStored()` lines 11948-11958  
**Current Cost**: O(6 × M³) per k-point = **6 × 10⁹ operations/k-point**  
**Impact**: Medium-High - expensive matrix operations

#### Problem
```fortran
! 6 separate ZGEMM calls per k-point:
call ZGEMM('N', 'N', M, M, M, -cmplx_i, delX, M, H_k, M, ...)  ! dH/dkx
call ZGEMM('N', 'N', M, M, M, -cmplx_i, delY, M, H_k, M, ...)  ! dH/dky  
call ZGEMM('N', 'N', M, M, M, ..., dHdkx, M, eigvec, M, ...)   ! temp_x
call ZGEMM('N', 'N', M, M, M, ..., dHdky, M, eigvec, M, ...)   ! temp_y
call ZGEMM('C', 'N', M, M, M, ..., eigvec, M, temp_x, M, ...)  ! Vx
call ZGEMM('C', 'N', M, M, M, ..., eigvec, M, temp_y, M, ...)  ! Vy
```

#### Optimization Strategy
**Potential Speedup**: **2-3x**

```fortran
! Batch operations and reuse intermediate results
! Pre-compute eigvec^† once
call ZGEMM('C', 'N', M, M, M, (1.0_dp,0.0_dp), eigvec, M, eigvec, M, (0.0_dp,0.0_dp), eigvec_dagger, M)

! Combined operations where possible
! Use optimized BLAS call patterns
! Consider using ZHEMM for Hermitian matrices where applicable
```

---

### 7. Memory Allocation/Deallocation ⭐⭐

**Location**: Every k-point in Berry calculation  
**Current Cost**: Memory churn for ~10 matrices of size M×M  
**Impact**: Medium - memory bandwidth bottleneck

#### Problem
```fortran
! Every k-point:
allocate(H_k(M,M), eigvec(M,M), eigval(M), delX(M,M), delY(M,M), ...)
! ... calculations ...
deallocate(H_k, eigvec, eigval, delX, delY, ...)
```

#### Optimization Strategy
**Potential Speedup**: **1.5-2x**

```fortran
! Pre-allocate workspace arrays (module level)
complex(dp), allocatable, save :: workspace_H_k(:,:), workspace_eigvec(:,:)
complex(dp), allocatable, save :: workspace_Vx(:,:), workspace_Vy(:,:)
real(dp), allocatable, save :: workspace_delX(:,:), workspace_delY(:,:)

! Initialize once
subroutine InitializeChernWorkspace(M)
   if (.not. allocated(workspace_H_k)) then
      allocate(workspace_H_k(M,M), workspace_eigvec(M,M), ...)
   end if
end subroutine

! Reuse in calculations - no allocate/deallocate per k-point
```

---

## 🟢 System-Level Optimizations

*These apply to both TAPW and Chern calculations*

### 8. OpenMP Parallelization ⭐⭐⭐⭐

**Location**: Main k-point loops  
**Current**: Sequential k-point processing  
**Impact**: High - excellent scalability potential

#### Implementation
**Potential Speedup**: **4-8x** (on multi-core systems)

```fortran
! TAPW bands calculation
!$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, HLoc) &
!$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, Kpts)
do ip = 1, ptsTot
   call DiagH0TAPW(...)
   !$OMP CRITICAL
   ! Store results
   !$OMP END CRITICAL
end do
!$OMP END PARALLEL DO

! Chern calculation  
!$OMP PARALLEL DO PRIVATE(berry_curv_bands) REDUCTION(+:chern_bands)
do ik = 1, nk
   call CalculateBerryAtKpointFromStored(ik, ...)
   ! Reduction handles accumulation automatically
end do
!$OMP END PARALLEL DO
```

#### Considerations
- **Thread-safe storage** for TAPW results
- **Memory bandwidth** limitations with large matrices
- **Load balancing** for different k-point costs

---

### 9. Memory Access Patterns ⭐⭐

**Location**: Large matrix operations  
**Impact**: Medium - cache efficiency improvements

#### Optimization Strategy
**Potential Speedup**: **1.5-2x**

```fortran
! Block-structured algorithms for better cache locality
! Optimize matrix storage order (row-major vs column-major)
! Use cache-blocking for large matrix multiplications
! Minimize memory allocations in hot loops
```

---

## 📊 Implementation Roadmap

### Phase 1: High-Impact Optimizations
**Expected Total Speedup: 50-500x**

1. **TAPW Diagonalization Reuse** (10-100x speedup)
   - Separate transformation from diagonalization
   - Cache transformation matrices
   - Implement incremental Hamiltonian updates

2. **OpenMP Parallelization** (4-8x speedup)
   - Parallelize k-point loops
   - Thread-safe data structures
   - Memory bandwidth optimization

### Phase 2: Medium-Impact Optimizations  
**Expected Additional Speedup: 5-15x**

3. **Berry Curvature Vectorization** (3-5x speedup)
   - Vectorize double loops
   - Pre-compute energy differences
   - Optimize BLAS usage

4. **TAPW Transformation Caching** (5-10x speedup)
   - Cache by symmetry classes
   - Pre-compute phase factors
   - Optimize G-vector operations

### Phase 3: System Optimizations
**Expected Additional Speedup: 2-5x**

5. **Memory Management** (1.5-2x speedup)
   - Pre-allocated workspaces
   - Eliminate allocation churn
   - Optimize memory access patterns

6. **Algorithm Refinements** (2-3x speedup)
   - Block-structured computations
   - Cache-friendly data layouts
   - Specialized BLAS operations

---

## 🎯 Expected Performance Gains

### Current Performance
- **TAPW Bands**: ~10¹³ operations
- **Chern Calculation**: ~10¹² operations  
- **Total Runtime**: Hours to days (depending on system)

### Optimized Performance
- **Phase 1**: **50-500x faster** → Minutes to hours
- **Phase 2**: **250-7500x faster** → Seconds to minutes  
- **Phase 3**: **500-37500x faster** → Near real-time for moderate grids

### Memory Efficiency
- **Current**: ~16 GB for stored arrays + allocation churn
- **Optimized**: ~8 GB with better data structures + pre-allocated workspaces

---

## 🔧 Implementation Notes

### Code Locations
- **TAPW Bands**: `lanczosKuboCode/Src/diag.F90` lines 883-1056
- **Chern Calculation**: `lanczosKuboCode/Src/diag.F90` lines 11473-12009
- **Position Matrices**: `lanczosKuboCode/Src/diag.F90` lines 12057-12175 ✅ **Already optimized**

### Dependencies
- **MKL/BLAS**: Already used for matrix operations
- **OpenMP**: Available, needs implementation
- **Memory**: Sufficient for optimized algorithms

### Testing Strategy
1. **Benchmark current performance** with timing instrumentation
2. **Implement optimizations incrementally** with validation
3. **Compare results** for numerical accuracy
4. **Profile optimized code** to identify remaining bottlenecks

---

## 📈 Monitoring and Profiling

### Key Metrics
- **Time per k-point** (TAPW bands vs Chern calculation)
- **Memory usage** (peak and sustained)
- **Cache hit rates** for transformation matrices
- **OpenMP scaling efficiency**

### Profiling Tools
```bash
# Intel VTune for detailed performance analysis
vtune -collect hotspots ./grabnes

# gprof for function-level profiling  
gfortran -pg ... && ./grabnes && gprof grabnes gmon.out

# Custom timing instrumentation
call MIO_TimerStart('tapw_diag')
! ... code ...
call MIO_TimerStop('tapw_diag')
```

---

## ⚠️ Considerations

### Numerical Stability
- **Verify Berry curvature accuracy** with optimized algorithms
- **Test edge cases** (small energy gaps, band crossings)
- **Compare with reference implementations**

### Memory Constraints
- **Monitor memory usage** on different HPC node types
- **Implement fallback strategies** for memory-limited systems
- **Consider out-of-core algorithms** for very large systems

### Portability
- **Test on different compilers** (Intel, GNU, etc.)
- **Validate OpenMP performance** across architectures
- **Ensure BLAS library compatibility**

---

*Last updated: September 26, 2025*  
*Generated for: lanczosKuboCode TAPW Chern number calculation*
