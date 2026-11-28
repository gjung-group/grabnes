# Electronic structure observables

## Density of states: Lanczos recursion vs. exact diagonalization

The density of states (DOS) can be computed in two ways, and the choice
between them depends on the system size.

**Exact diagonalization** is the more accurate option, but due to its poor
scaling with system size one is limited to systems containing about 10,000
atoms on a low-memory machine and 50,000 atoms on a high-memory machine. The
DOS is accumulated over a Monkhorst–Pack grid (`KGrid`) with a Lorentzian
broadening `Epsilon`, so the k-point sampling must be dense enough to
properly populate the DOS.

**Lanczos recursion** easily deals with systems containing 20 million atoms
on a low-memory machine, with the precision controlled by the number of
recursion steps included in the continued fraction. The method tridiagonalizes
the original Hamiltonian $H$ into a smaller Hamiltonian $\hat{H}$ on which
the operator $\delta(E-\hat{H})$ corresponding to the DOS can be efficiently
evaluated numerically:

$$
\delta(E-\hat{H}) =
\cfrac{1}{E + i \eta - a_1 - \cfrac{b_1^2}{\dots \cfrac{b_2^2}{\dots \cfrac{\dots}{E+ i \eta - a_{N-1} - \cfrac{b_{N-1}^2}{{E+i \eta - a_N - b_N^2 \times Term}}}}}}
$$

where $a_N$ are the diagonal recursion coefficients and $b_N$ the off-diagonal
coefficients in the tridiagonal representation, while $\eta$ is a numerical
broadening in the Lorentzian basis function
$\frac{1}{E+i\eta-\hat{H}}$. The expression is terminated with `Term` chosen
as 1, but more advanced methods exist in the literature. The coefficients are
computed recursively, so the computational cost scales with the number of
recursion coefficients. See
[this reference](https://www.sciencedirect.com/science/article/pii/S0370157320304245)
for further details.

```{note}
When using Lanczos recursion for the DOS, a large enough supercell is needed
to sample enough states in real space and thus properly populate the DOS —
analogous to the requirement of adding sufficient reciprocal-space k-points
when using exact diagonalization. Larger systems give a much smoother curve
for a fixed number of recursion steps and energy broadening.
```

The two approaches give similar results at the same energy resolution. Note
that the Lanczos DOS output is renormalized by the number of atoms: to match
it with the diagonalization DOS, multiply the Lanczos DOS by two times the
number of atoms and divide by the first nearest-neighbor hopping (the
`TB.Hopping` value of the input file); the energy axis is also given in units
of the first nearest-neighbor hopping.

## Band structures

With `Calculate.Bands .true.` the code diagonalizes the Hamiltonian on a path
through the high-symmetry points of the (super)cell Brillouin zone, defined
in reduced coordinates via the `&begin Bands.Path ... &end Bands.Path` code
block. The result is written to `generate.bands`.

## Spectral functions and band unfolding

The code calculates the band-unfolded spectral band structure
(`Calculate.spectral`), which projects the supercell spectral weight back
onto the folded Brillouin zone of the constituent layers. This is the natural
quantity to compare with ARPES experiments: GRABNES simulations have been
shown to be in good agreement with nano-ARPES data on twisted bilayer
graphene, see the
[Nature Physics publication](https://www.nature.com/articles/s41567-020-0974-x)
and [](../benchmarks/index.md). The local density of states and the
probability density based on the wavefunction distribution are also
available.

## Magnetic field and Landau levels

A perpendicular magnetic field enters the Hamiltonian through the standard
Peierls phase substitution, see for instance
[A. Cresti, Phys. Rev. B 103, 045402 (2021)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.045402).
Because the Peierls phases destroy the simple translational periodicity that
exact diagonalization relies on, magnetic-field calculations are performed
with Lanczos recursion on large supercells.

The field is controlled either physically (`MagField <B>` in Tesla) or via
the integer that relates the flux quantum to the simulation cell size
(`MagField.Integer <n>`). Scanning the field with
`MagField.InitialInteger`, `MagField.FinalInteger` and `MagField.IntegerStep`
produces the data for Hofstadter butterfly or Wannier diagrams.

Combining a magnetic field with on-site Anderson disorder
(`Anderson .true.`, `AndersonAmp <W>`) illustrates the physics of Landau
level broadening: the Landau levels become wider due to the presence of more
localized disorder states in the tails around the center of the Landau level
where the extended states reside. If the disorder is increased further, the
Landau levels cannot form anymore and a regular Anderson insulator develops
at all energies, where all states are localized.

## Transport

The Kubo–Greenwood conductivity is computed with the Lanczos machinery
(`Kubo.Calc .true.`), reusing the recursion coefficients of the DOS
calculation.
