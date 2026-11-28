# Citing GRABNES

This page is quoted from the
[README](https://github.com/gjung-group/grabnes/blob/main/README.md) and
[`CITATION.cff`](https://github.com/gjung-group/grabnes/blob/main/CITATION.cff)
of the code repository.

## Citation

Please cite GRABNES by its authors, the repository, and the exact version or
Git commit you used. The metadata are in [`CITATION.cff`](https://github.com/gjung-group/grabnes/blob/main/CITATION.cff), which
GitHub displays under "Cite this repository":

> J. Jung, R. Martinez-Gordillo, and N. Leconte, *GRABNES: Graphene and Boron
> Nitride Electronic Structure*, https://github.com/gjung-group/grabnes,
> commit `<hash>`.

A software paper and an archived release with a DOI do not exist yet; this
section will name them when they do.

## Publications to cite with the code

Please also cite the publications behind the methods and models you use.

Methods of the code (real-space lattice models, Lanczos recursion):

- N. Leconte, A. Ferreira, and J. Jung, *Efficient multiscale lattice
  simulations of strained and disordered graphene*, Semicond. Semimet. **95**,
  35-99 (2016), [doi:10.1016/bs.semsem.2016.04.002](https://doi.org/10.1016/bs.semsem.2016.04.002).

Twisted and multilayer graphene:

- N. Leconte, S. Javvaji, J. An, A. Samudrala, and J. Jung, *Relaxation effects
  in twisted bilayer graphene: A multiscale approach*, Phys. Rev. B **106**,
  115410 (2022), [doi:10.1103/PhysRevB.106.115410](https://doi.org/10.1103/PhysRevB.106.115410)
  (twisted-bilayer model; example 08).
- N. Leconte, Y. Park, J. An, A. Samudrala, and J. Jung, *Electronic structure
  of lattice relaxed alternating twist tNG-multilayer graphene: from few layers
  to bulk AT-graphite*, 2D Mater. **9**, 044002 (2022),
  [doi:10.1088/2053-1583/ac8a00](https://doi.org/10.1088/2053-1583/ac8a00)
  (alternating-twist multilayers; example 09).
- K. Shin, N. Leconte, J. Jung, and H. Min, *Supermoire reconstruction and
  domain-resolved valley topology in helical multilayer graphene*, 2D Mater.
  **13**, 045021 (2026), [doi:10.1088/2053-1583/aeb140](https://doi.org/10.1088/2053-1583/aeb140)
  (helical multilayers).
- P. Moon and M. Koshino, *Energy spectrum and quantum Hall effect in twisted
  bilayer graphene*, Phys. Rev. B **85**, 195458 (2012),
  [doi:10.1103/PhysRevB.85.195458](https://doi.org/10.1103/PhysRevB.85.195458)
  (two-center interlayer hopping; examples 03 and 04).

Graphene on hexagonal boron nitride:

- J. Jung, A. Raoux, Z. Qiao, and A. H. MacDonald, *Ab initio theory of moire
  superlattice bands in layered two-dimensional materials*, Phys. Rev. B **89**,
  205414 (2014), [doi:10.1103/PhysRevB.89.205414](https://doi.org/10.1103/PhysRevB.89.205414)
  (effective moire-potential model; examples 05 to 07).
- N. Leconte and J. Jung, *Commensurate and incommensurate double moire
  interference in graphene encapsulated by hexagonal boron nitride*, 2D Mater.
  **7**, 031005 (2020), [doi:10.1088/2053-1583/ab891a](https://doi.org/10.1088/2053-1583/ab891a)
  (graphene encapsulated by hBN, double moire).
- J. An, N. Leconte, S. Javvaji, Y. Park, and J. Jung, *Atomistic substrate
  relaxation effects in the band gaps of graphene on hexagonal boron nitride*,
  [arXiv:2602.18893](https://arxiv.org/abs/2602.18893) (2026) (atomistic
  graphene/hBN with a relaxed substrate).
- W. Shi, S. Kahn, N. Leconte, T. Taniguchi, K. Watanabe, M. Crommie, J. Jung,
  and A. Zettl, *High-order fractal quantum oscillations in graphene/BN
  superlattices in the extreme doping limit*, Phys. Rev. Lett. **130**, 186204
  (2023), [doi:10.1103/PhysRevLett.130.186204](https://doi.org/10.1103/PhysRevLett.130.186204)
  (graphene/hBN in a magnetic field).
- Z. Wang, P. M. Perez-Piskunow, C. P. Y. Wong, M. Holwill, H. Sun, W. Fu,
  J. Hu, Y. Zheng, X. Zhang, T. Taniguchi, K. Watanabe, A. Ariando, G. Yu,
  L. Lin, K. E. J. Goh, S. Roche, J. Jung, K. S. Novoselov, and N. Leconte,
  *Beyond the flat-band paradigm: open orbits as a route to correlated states
  in dispersive moire minibands*, Proc. Natl. Acad. Sci. USA, accepted (2026),
  [arXiv:2610.06110](https://arxiv.org/abs/2610.06110) (graphene/hBN
  minibands and their orbits in a magnetic field).

## Authors and maintenance

The authors of GRABNES are Jeil Jung, Rafael Martinez-Gordillo, and Nicolas
Leconte. The code grew out of a tight-binding and Kubo transport program
written by Rafael Martinez-Gordillo, whose input/output and math libraries it
still uses, and was developed further in the group of Jeil Jung.

Nicolas Leconte is the current developer and the maintainer of this
repository; questions and reports should be addressed to him through the
issue tracker.

## License

GRABNES is free software, distributed under the terms of the GNU General
Public License, version 3 or (at your option) any later version. The full
text is in [`LICENSE`](https://github.com/gjung-group/grabnes/blob/main/LICENSE). The program comes without any warranty.

Some files were written by other authors and keep their own notices (GNU GPL,
GNU LGPL, two-clause BSD); they are listed in
[`THIRD_PARTY_LICENSES.md`](https://github.com/gjung-group/grabnes/blob/main/THIRD_PARTY_LICENSES.md). The reasoning behind
the choice of license is recorded in
[](../docs/development/licensing-audit.md).
