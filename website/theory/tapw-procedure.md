(tapw)=
# TAPW band-structure procedure

```{note}
Source: *Twisted Angle Planar Wave (TAPW) Method Implementation — Band
Structure Calculation Procedure*, GRABNES contributors, pre-release technical
note.
```

## Introduction

The Twisted Angle Planar Wave (TAPW) method is implemented in the `DiagBands`
routine for calculating electronic band structures of moiré systems. This
document provides a comprehensive overview of the TAPW procedure from entry
into `DiagBands` through generation of the `generate.bands` file, excluding
Chern number calculations.

## Overview of TAPW flow

The TAPW calculation follows this sequence:

1. **Parameter Setup**: Read TAPW parameters and initialize coordinate systems
2. **G-vector Generation**: Create reciprocal space grid around reference K-point
3. **Basis Construction**: Build TAPW basis functions using plane waves
4. **Hamiltonian Projection**: Project tight-binding Hamiltonian onto TAPW basis
5. **Diagonalization**: Solve eigenvalue problem for each k-point
6. **Output Generation**: Write band structure to `generate.bands`

## Detailed implementation

### Entry point and parameter setup

The TAPW calculation begins when `useTAPW = .true.` in the `DiagBands` routine
(*TAPW Entry Point*):

```fortran
! In DiagBands routine
else if (useTAPW) then
   call MIO_Print('calling TAPW routine for '//trim(num2str(ip)),'diag')
   call DiagH0TAPW(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, 
                   maxNeigh, hopp, NList, Nneigh, neighCell, neig)
end if
E(:, is, ip) = ELoc
```

Key parameters are read from input:

- `Diag.N_G`: Number of G-vectors (`tapwNG`)
- `Diag.MoireAngle`: Rotation angle for graphene lattice

### Coordinate system construction

#### Graphene reciprocal lattice

The code constructs rotated graphene reciprocal vectors `rG`:

$$
\mathbf{a}_1 = a_G \begin{pmatrix} 1 \\ 0 \end{pmatrix}, \quad
\mathbf{a}_2 = a_G \begin{pmatrix} \cos(60°) \\ \sin(60°) \end{pmatrix}
$$

After rotation by `moireAngle`:

$$
\mathbf{a}'_i = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} \mathbf{a}_i
$$

Reciprocal vectors are computed as:

$$
\mathbf{rG} = 2\pi \cdot (\mathbf{A}^T)^{-1}
$$

```fortran
! Build graphene lattice vectors (60-degree structure)
sGlattice(1,1) = 1.0_dp * aG                    ! (aG, 0)
sGlattice(2,1) = 0.0_dp
sGlattice(1,2) = cos(60.0_dp*pi/180.0_dp) * aG  ! (aG*cos60°, aG*sin60°)
sGlattice(2,2) = sin(60.0_dp*pi/180.0_dp) * aG

! Apply rotation
rotation_angle = moireAngle * pi / 180.0_dp
sGlattice = matmul(reshape((/ cos_rot, sin_rot, -sin_rot, cos_rot /), 
                           (/2,2/)), sGlattice)

! Convert to reciprocal space: rG = 2π * inverse(lattice.T)
rG = 2.0_dp * pi * transpose(sGlattice_inv)
```

#### Reference K-point

The K-point in graphene's Brillouin zone is calculated as:

$$
\mathbf{k}_{ref} = \frac{2}{3}\mathbf{b}_1 + \frac{1}{3}\mathbf{b}_2
$$

```fortran
real(dp), parameter :: refF(2) = (/ 2.0_dp/3.0_dp, 1.0_dp/3.0_dp /)
k_ref = refF(1)*rG(:,1) + refF(2)*rG(:,2)
```

### G-vector generation

The TAPW method generates a set of reciprocal lattice vectors around the
reference K-point. Two methods are implemented:

#### Hexagonal shell method

```fortran
call generate_shifted_G_list(rcell, k_ref, NGrange, Gx, Gy, NG)
```

This generates G-vectors as:

$$
\mathbf{G}_{n_1,n_2} = n_1 \mathbf{b}_{moire,1} + n_2 \mathbf{b}_{moire,2}
$$

The vectors are:

1. Generated on a grid around `k_ref`
2. Sorted by distance from the nearest lattice point
3. Selected to keep the `N_G` closest vectors
4. Filtered for uniqueness

#### Triangular truncation method

Alternative method following Python implementation:

$$
\mathbf{G}_{frac} = \begin{pmatrix} j + k_{ref,x} \\ i + k_{ref,y} \end{pmatrix}
$$

Selected within triangular Brillouin zone boundaries.

### TAPW basis construction

#### Basis functions

The TAPW basis functions are plane waves modulated by atomic orbitals:

$$
\phi^{\alpha}_{\mathbf{G}}(\mathbf{r}) = \sum_{\mathbf{R}} e^{i(\mathbf{k} + \mathbf{G}) \cdot (\mathbf{R} + \boldsymbol{\tau}_\alpha)} \psi_\alpha(\mathbf{r} - \mathbf{R} - \boldsymbol{\tau}_\alpha)
$$

where:

- $\alpha$ labels atomic species/sublattices
- $\mathbf{G}$ are the selected reciprocal lattice vectors
- $\mathbf{R}$ are lattice translation vectors
- $\boldsymbol{\tau}_\alpha$ are atomic positions within unit cell

#### X-matrix construction

The transformation matrix $X$ connects atomic orbitals to TAPW basis:

```fortran
call build_X(X, x, y, label, Gx, Gy, N, NG, Nlabel)
```

Matrix elements:

$$
X_{i,\alpha\mathbf{G}} = e^{i\mathbf{G} \cdot \mathbf{r}_i} \delta_{\alpha,\text{species}(i)}
$$

Implemented as:

```fortran
do i = 1, N
   do iG = 1, NG
      do iLabel = 1, Nlabel
         if (label(i) == iLabel) then
            idx = (iLabel-1)*NG + iG
            phase = Gx(iG)*x(i) + Gy(iG)*y(i)
            X(i, idx) = cmplx(cos(phase), sin(phase), dp)
         end if
      end do
   end do
end do
```

### Hamiltonian projection

#### Sparse matrix construction

The tight-binding Hamiltonian is first constructed in sparse CSR format:

```fortran
call initialize_sparse_matrix(N, maxN, H0, hopp, NList, Nneigh, 
                             neighCell, ns, is, KLoc, cell, 
                             row_ptr, col_ind, values, sigma)
```

Matrix elements include:

$$
H_{ij} = \epsilon_i \delta_{ij} + \sum_{\mathbf{R}} t_{ij}(\mathbf{R}) e^{i\mathbf{k} \cdot \mathbf{R}}
$$

where $t_{ij}(\mathbf{R})$ are hopping parameters and $\epsilon_i$ are on-site
energies.

#### TAPW projection

The projected Hamiltonian in TAPW space:

$$
H^{TAPW} = X^\dagger H X
$$

```fortran
call transform_sparse_hamiltonian(N, M, row_ptr, col_ind, values, X, Hproj)
```

Matrix dimensions: $M = N_G \times N_{label}$ where $N_{label}$ typically
equals 2 (sublattices A and B).

### Eigenvalue solution

The projected Hamiltonian is diagonalized using LAPACK:

$$
H^{TAPW} \psi^{TAPW}_n = E_n \psi^{TAPW}_n
$$

```fortran
allocate(eigvals(M))
lwork = 2*M
allocate(ZWorkLoc(2*M), DWorkLoc(3*M))

call ZHEEV('V','U', M, Hproj, M, eigvals, ZWorkLoc, lwork, DWorkLoc, info)

if (info /= 0) then
   print *, "Diagonalization failed: ZHEEV info =", info
   stop
end if
```

### K-point loop and progress tracking

The TAPW calculation is performed for each k-point in the path:

```fortran
do ip = 1, ptsTot
   ! Progress tracking for TAPW (show every 10%)
   if (useTAPW .and. modulo(ip, max(1, int(ptsTot/10))) == 0) then
      call MIO_Print('TAPW progress: '//trim(num2str(nint(100.0_dp*ip/ptsTot)))//
                     '% ('//trim(num2str(ip))//'/'//trim(num2str(ptsTot))//
                     ' k-points)', 'diag')
   end if
   
   do is = 1, nspin
      KptsLoc = Kpts(:, ip)
      call DiagH0TAPW(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, 
                      maxNeigh, hopp, NList, Nneigh, neighCell, neig)
      E(:, is, ip) = ELoc
   end do
end do
```

## Output generation

### Band structure file

After all k-points are calculated, the results are written to
`generate.bands`:

```fortran
! Write TAPW matrix size header
if (useTAPW) then
   write(uu,'(A,I0)') '# TAPW matrix size: M = ', M_tapw
end if

! Write eigenvalues for each k-point
do ip = 1, ptsTot
   do is = 1, nspin
      do i = 1, nAt
         write(uu, *) E(i, is, ip)
      end do
   end do
end do
```

### Debug output files

Several debug files are generated:

- `brillouin_zones_debug.dat`: Brillouin zone vertices and K-points
- `moire_bz_debug.dat`: Moiré BZ information
- `g_vectors_debug.dat`: Generated G-vectors
- `kpath_debug_*`: K-path coordinates

## Key equations summary

### Core TAPW equations

$$
\begin{aligned}
\text{Basis functions:} \quad & \phi^{\alpha}_{\mathbf{G}}(\mathbf{r}) = \sum_{\mathbf{R}} e^{i(\mathbf{k} + \mathbf{G}) \cdot (\mathbf{R} + \boldsymbol{\tau}_\alpha)} \psi_\alpha(\mathbf{r} - \mathbf{R} - \boldsymbol{\tau}_\alpha) \\
\text{Transformation matrix:} \quad & X_{i,\alpha\mathbf{G}} = e^{i\mathbf{G} \cdot \mathbf{r}_i} \delta_{\alpha,\text{species}(i)} \\
\text{Projected Hamiltonian:} \quad & H^{TAPW} = X^\dagger H X \\
\text{Eigenvalue equation:} \quad & H^{TAPW} \psi^{TAPW}_n = E_n \psi^{TAPW}_n
\end{aligned}
$$

### Coordinate systems

$$
\begin{aligned}
\text{Graphene lattice:} \quad & \mathbf{a}_1 = a_G(1,0), \quad \mathbf{a}_2 = a_G(\cos 60°, \sin 60°) \\
\text{Reciprocal lattice:} \quad & \mathbf{rG} = 2\pi \cdot (\mathbf{A}^T)^{-1} \\
\text{K-point:} \quad & \mathbf{k}_{ref} = \frac{2}{3}\mathbf{b}_1 + \frac{1}{3}\mathbf{b}_2 \\
\text{G-vectors:} \quad & \mathbf{G}_{n_1,n_2} = n_1 \mathbf{b}_{1} + n_2 \mathbf{b}_{2}
\end{aligned}
$$

## Implementation parameters

| Parameter | Variable | Description |
| --- | --- | --- |
| `Diag.N_G` | `tapwNG` | Number of G-vectors (default: 180) |
| `Diag.MoireAngle` | `moireAngle` | Rotation angle in degrees |
| `refF` | `(2/3, 1/3)` | Fractional K-point coordinates |
| `Nlabel` | 2 | Number of sublattices (A, B) |
| `M` | `NG × Nlabel` | TAPW matrix dimension |

## Computational complexity

- **G-vector generation**: $O(N_G \log N_G)$ for sorting
- **X-matrix construction**: $O(N \times N_G \times N_{label})$
- **Hamiltonian projection**: $O(N_{sparse} \times M)$ where $N_{sparse}$ is
  number of non-zero elements
- **Diagonalization**: $O(M^3)$ per k-point
- **Total per k-point**: $O(M^3)$ dominated by diagonalization

For typical parameters ($N_G = 180$, $N_{label} = 2$), $M = 360$, making
diagonalization the computational bottleneck.

## Verification and debugging

The implementation includes comprehensive verification:

- Direct-reciprocal lattice duality checks
- G-vector uniqueness and ordering verification
- Hamiltonian Hermiticity validation
- Matrix norm consistency checks
- Coordinate system orientation verification

## Conclusion

The TAPW implementation provides a systematic approach to calculating
electronic band structures in moiré systems by:

1. Constructing appropriate reciprocal space grids
2. Building plane wave basis functions
3. Projecting tight-binding Hamiltonians onto reduced subspaces
4. Efficiently solving eigenvalue problems
5. Generating standard band structure output

The method effectively reduces the computational complexity from the full
supercell size to a manageable projected space while maintaining the essential
physics of the moiré system.
