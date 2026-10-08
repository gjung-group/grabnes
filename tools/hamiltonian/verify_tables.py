#!/usr/bin/env python3
"""Independent check of the tight-binding Hamiltonian assembled by GRABNES.

Run a band calculation with ``WriteDataFiles .true.`` and pass the run
directory. The script reads the tables written by the solver

    v, neighCell.dat, neighD.dat   neighbor list, lattice translations, displacements
    <prefix>.s.mag                 hopping of every neighbor entry (units of g0, sign reversed)
    <prefix>.pos, .cell, .e        positions, lattice vectors, on-site energies and layers
    <prefix>.bands                 eigenvalues along the band path

and compares them with a Hamiltonian that is built here from the atomic
positions alone: every pair (i, m, R) inside the search radii is enumerated
exhaustively and its matrix element is evaluated from the model parameters
given on the command line. Nothing of the solver's neighbor list, lattice
translations, or hopping values enters this second Hamiltonian.

Bloch convention of GRABNES (``DiagHam``): with R the lattice translation of
neighbor entry j of atom i (the neighbor sits at r_m + R),

    H[m, i](k) = sum_j H_j exp(-i k.R),      H_j = -hopp(j, i) * g0 .

The phases only contain lattice vectors, so H(k + G) = H(k) element by element.
Moving an atom by a lattice vector T multiplies its row and column by
exp(+-i k.T): a unitary change of basis that leaves the eigenvalues unchanged.
Both properties are tested, as is the convention with phases exp(i k.d) of the
full displacement d, which has the same spectrum.

Checks, with the residual of each one printed:

1. stored displacement = r_m + R - r_i for every entry;
2. solver entries = exhaustive enumeration (missing, extra, duplicates);
3. every entry (i, m, R) has the partner (m, i, -R) with the conjugate value;
4. every stored matrix element = the independent model value;
5. H(k) Hermitian, and equal to the independent H(k), at Gamma, K, M, generic
   and seeded random k-points;
6. invariance under reciprocal lattice vectors and under moving atoms to other
   unit cells;
7. eigenvalues of the independent H(k) = bands written by the solver.

Requires NumPy. Exit status 1 if a residual exceeds its limit.
"""

import argparse
import re
import sys

import numpy as np

# Neighbor shells of the honeycomb lattice in units of the bond length; the
# ninth value only places the outermost search radius. Kept independent of the
# table in neigh.F90 on purpose.
SHELL_RADIUS = np.sqrt(np.array([1.0, 3.0, 4.0, 7.0, 9.0, 12.0, 13.0, 16.0, 19.0]))
DZ_INTRA, DZ_INTER, DZ_SECOND = 1.5, 4.5, 7.5


def load(run, prefix):
    lines = open(f"{run}/v").read().split("\n")
    neighbors = []
    for i in range(0, len(lines) - 1, 2):
        if not lines[i].strip():
            break
        row = [int(x) - 1 for x in lines[i + 1].split()]
        assert len(row) == int(lines[i]), "inconsistent neighbor file v"
        neighbors.append(row)

    def rows(name, convert=float):
        return [[convert(x) for x in line.split()] for line in open(f"{run}/{name}") if line.strip()]

    pair = re.compile(r"\(\s*([-+0-9.Ee]+)\s*,\s*([-+0-9.Ee]+)\s*\)")
    hop = [[complex(float(a), float(b)) for a, b in pair.findall(line)]
           for line in open(f"{run}/{prefix}.s.mag") if line.strip()]
    images = np.array(rows("neighCell.dat", int))
    disp = np.array(rows("neighD.dat"))
    cell = np.array(rows(f"{prefix}.cell"))
    onsite = rows(f"{prefix}.e")
    counts = [len(row) for row in neighbors]
    assert [len(h) for h in hop] == counts and len(images) == sum(counts) == len(disp), "inconsistent tables"
    owner = np.repeat(np.arange(len(neighbors)), counts)
    return dict(
        owner=owner, neighbor=np.concatenate([np.array(r, int) for r in neighbors]),
        image=images, disp=disp, hop=np.concatenate([np.array(h) for h in hop]),
        pos=np.array(rows(f"{prefix}.pos")), ucell=cell[:3].T, rcell=cell[3:].T,
        h0=np.array([r[0] for r in onsite]), species=np.array([int(r[1]) for r in onsite]),
        layer=np.array([int(r[2]) for r in onsite]),
    )


def read_bands(path):
    lines = open(path).read().split("\n")
    bands, spins, points = map(int, lines[3].split())
    values = np.array([float(x) for line in lines[4:] for x in line.split()])
    values = values.reshape(points, bands * spins + 1)
    return values[:, 0], values[:, 1:]


def path_k_points(rcell, corners, path_length):
    """k-points of the band file. The solver spaces them uniformly on each
    segment; the six-decimal path coordinate of the file is only used to count
    the points per segment, so that the k-points are exact."""
    corners = [rcell @ np.array(c) for c in corners]
    ends = np.cumsum([np.linalg.norm(b - a) for a, b in zip(corners[:-1], corners[1:])])
    x = np.asarray(path_length)
    result = [corners[0]]
    first = 1
    for a, b, end in zip(corners[:-1], corners[1:], ends):
        last = int(np.searchsorted(x, end + 1e-5))      # points up to this corner
        steps = last - first
        for j in range(1, steps + 1):
            result.append(a + (b - a) * j / steps)
        first = last
    assert len(result) == len(x), "band path does not match the given corners"
    return result


class Model:
    """Hamiltonian defined by the command-line parameters only."""

    def __init__(self, args, bond):
        self.bond = bond
        self.two_center_shells = args.koshino_intralayer
        self.shell_energy = args.intralayer or []
        n = self.two_center_shells or len(self.shell_energy) or args.structure_only
        self.rc_intra = 0.5 * (SHELL_RADIUS[n - 1] + SHELL_RADIUS[n]) * bond
        self.shell_edges = 0.5 * (SHELL_RADIUS[:n] + SHELL_RADIUS[1:n + 1]) * bond
        self.rc_inter = args.interlayer_cutoff
        self.vpi, self.vsigma = -args.vpppi0, args.vppsigma0
        self.d0, self.delta = args.interlayer_distance, args.decay
        self.second = args.second_layer
        self.periodic_z = args.periodic_z

    def classify(self, d):
        """0: not a pair, 1: intralayer, 2: interlayer (arrays over displacements)."""
        rho = np.hypot(d[..., 0], d[..., 1])
        adz = abs(d[..., 2])
        intra = (adz < DZ_INTRA) & (rho < self.rc_intra)
        top = DZ_SECOND if self.second else DZ_INTER
        inter = (adz > DZ_INTRA) & (adz < top) & (rho < self.rc_inter)
        if not self.second:
            inter &= adz < DZ_INTER
        return np.where(intra, 1, np.where(inter, 2, 0))

    def element(self, d, kind):
        """Matrix element (eV) of pairs with displacement d and class kind."""
        rho = np.hypot(d[..., 0], d[..., 1])
        dist = np.linalg.norm(d, axis=-1)
        safe = np.where(dist > 0, dist, 1.0)
        cz = (d[..., 2] / safe) ** 2
        inter = (self.vpi * np.exp(-(dist - self.bond) / self.delta) * (1 - cz)
                 + self.vsigma * np.exp(-(dist - self.d0) / self.delta) * cz)
        if self.two_center_shells:
            intra = inter          # KoshinoIntralayer: the same two-center form within a layer
        else:
            shell = np.searchsorted(self.shell_edges, rho)
            intra = np.array(self.shell_energy + [0.0])[np.minimum(shell, len(self.shell_energy))]
        return np.where(kind == 1, intra, np.where(kind == 2, inter, 0.0))

    def enumerate(self, pos, ucell):
        """All pairs (i, m, R != self) inside the search radii: arrays i, m, R, d, H."""
        rmax = max(self.rc_intra, self.rc_inter)
        area = abs(ucell[0, 0] * ucell[1, 1] - ucell[0, 1] * ucell[1, 0])
        frac = np.linalg.solve(ucell[:2, :2], pos[:, :2].T).T
        spread = frac.max(axis=0) - frac.min(axis=0)
        reach = [int(np.ceil(rmax * np.linalg.norm(ucell[:2, 1 - a]) / area + spread[a])) + 1 for a in (0, 1)]
        reach_z = 0
        if self.periodic_z:
            top = (DZ_SECOND if self.second else DZ_INTER) + np.ptp(pos[:, 2])
            reach_z = int(np.ceil(top / ucell[2, 2])) + 1
        out = [[], [], [], [], []]
        for a in range(-reach[0], reach[0] + 1):
            for b in range(-reach[1], reach[1] + 1):
                for c in range(-reach_z, reach_z + 1):
                    image = np.array([a, b, c])
                    d = pos[None, :, :] + ucell @ image - pos[:, None, :]
                    kind = self.classify(d)
                    if a == 0 and b == 0 and c == 0:
                        np.fill_diagonal(kind, 0)
                    i, m = np.nonzero(kind)
                    out[0].append(i)
                    out[1].append(m)
                    out[2].append(np.tile(image, (len(i), 1)))
                    out[3].append(d[i, m])
                    out[4].append(self.element(d[i, m], kind[i, m]))
        return (np.concatenate(out[0]), np.concatenate(out[1]), np.concatenate(out[2]),
                np.concatenate(out[3]), np.concatenate(out[4]))


def bloch(n, owner, neighbor, value, phase_vector, k, onsite):
    """H[m, i](k) = sum value * exp(-i k.phase_vector): the convention of DiagHam."""
    h = np.zeros((n, n), complex)
    np.add.at(h, (neighbor, owner), value * np.exp(-1j * (phase_vector @ k)))
    return h + np.diag(onsite)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_dir")
    parser.add_argument("--prefix", default="generate")
    parser.add_argument("--lattice-parameter", type=float, default=2.46)
    parser.add_argument("--g0", type=float, help="energy unit in eV (default: 12.14 - 3.72 a)")
    parser.add_argument("--intralayer", type=lambda s: [float(x) for x in s.split(",")],
                        help="intralayer matrix elements in eV, one per neighbor shell, e.g. -2.9888,0,0")
    parser.add_argument("--onsite", default="", metavar="SPECIES:E,...",
                        help="on-site energies of the model in eV by species number, e.g. 3:3.09,4:-1.89 "
                             "(default: zero)")
    parser.add_argument("--structure-only", type=int, metavar="SHELLS",
                        help="no model: only check the neighbor geometry for this many intralayer "
                             "shells, the translations, the reverse partners and Hermiticity of the "
                             "stored values, and list the stored elements by species and distance")
    parser.add_argument("--koshino-intralayer", type=int, metavar="SHELLS",
                        help="KoshinoIntralayer: intralayer elements from the two-center form with "
                             "vpppi0, on this many neighbor shells (instead of --intralayer)")
    parser.add_argument("--vpppi0", type=float, default=3.5)
    parser.add_argument("--vppsigma0", type=float, default=0.48)
    parser.add_argument("--interlayer-distance", type=float, default=3.34)
    parser.add_argument("--decay", type=float, help="decay length in A (default: 0.184 a)")
    parser.add_argument("--interlayer-cutoff", type=float,
                        help="in-plane interlayer search radius in A (default: 1.1 a_cc * 6.2); "
                             "1.0 corresponds to Neigh.LayerNeighbors 0")
    parser.add_argument("--periodic-z", action="store_true", help="system periodic along z (Bulk, nonBulkSmall)")
    parser.add_argument("--second-layer", action="store_true", help="addSecondLayerInteractions")
    parser.add_argument("--path", default="2/3,1/3;0,0;1/2,0;2/3,1/3",
                        help="corners of the band path in reciprocal coordinates")
    parser.add_argument("--seed", type=int, default=20261008)
    parser.add_argument("--tol-matrix", type=float, default=1e-9,
                        help="limit for matrix-element and Hermiticity residuals (eV)")
    parser.add_argument("--tol-bands", type=float, default=1e-6,
                        help="limit for the comparison with the six-decimal band file (eV)")
    args = parser.parse_args()

    if sum(x is not None for x in (args.intralayer, args.koshino_intralayer, args.structure_only)) != 1:
        parser.error("give exactly one of --intralayer, --koshino-intralayer and --structure-only")
    a = args.lattice_parameter
    bond = a / np.sqrt(3.0)
    g0 = args.g0 if args.g0 is not None else 12.14 - 3.72 * a
    if args.decay is None:
        args.decay = 0.184 * a
    if args.interlayer_cutoff is None:
        args.interlayer_cutoff = bond * 1.1 * 6.2
    t = load(args.run_dir, args.prefix)
    pos, ucell, rcell = t["pos"], t["ucell"], t["rcell"]
    n = len(pos)
    model = Model(args, bond)
    failures = []

    def report(label, value, limit=None, unit="eV"):
        flag = ""
        if limit is not None and not value <= limit:
            flag = f"   <-- above the limit {limit:g}"
            failures.append(label)
        print(f"  {label:<62s} {value:10.3e} {unit}{flag}")

    print(f"atoms: {n}   neighbor entries: {len(t['owner'])}   g0 = {g0:.5f} eV")
    print(f"cell: |A1| = {np.linalg.norm(ucell[:, 0]):.5f} A, |A2| = {np.linalg.norm(ucell[:, 1]):.5f} A; "
          f"intralayer radius {model.rc_intra:.4f} A, "
          f"interlayer radius {model.rc_inter:.4f} A")

    solver_value = -t["hop"] * g0                      # H_j of every entry (eV)
    species_energy = {int(a_): float(b_) for a_, b_ in (item.split(":") for item in args.onsite.split(",") if item)}
    model_onsite = np.array([species_energy.get(int(sp_), 0.0) for sp_ in t["species"]])
    lattice = t["image"] @ ucell.T                     # R of every entry

    print("1. lattice translations")
    expected_disp = pos[t["neighbor"]] + lattice - pos[t["owner"]]
    report("max |stored displacement - (r_m + R - r_i)|", abs(expected_disp - t["disp"]).max(), 1e-9, "A")

    print("2. completeness of the neighbor list")
    ei, em, eimage, edisp, evalue = model.enumerate(pos, ucell)
    solver_keys = {}
    duplicates = 0
    for index, key in enumerate(zip(t["owner"], t["neighbor"], *t["image"].T)):
        duplicates += key in solver_keys
        solver_keys[key] = index
    enum_keys = {key: index for index, key in enumerate(zip(ei, em, *eimage.T))}
    missing = [enum_keys[key] for key in enum_keys if key not in solver_keys]
    extra = [solver_keys[key] for key in solver_keys if key not in enum_keys]
    print(f"  pairs inside the search radii (exhaustive enumeration): {len(enum_keys)}; "
          f"in the solver list: {len(solver_keys)}")
    report("pairs missing from the solver list", len(missing), 0, "")
    report("solver entries outside the search radii", len(extra), 0, "")
    report("duplicate solver entries", duplicates, 0, "")
    if missing:
        report("largest |H| among the missing pairs", abs(evalue[missing]).max())
    if extra:
        report("largest |H| among the extra entries", abs(solver_value[extra]).max())

    print("3. reverse partners (H_ij(R) = conj(H_ji(-R)))")
    unpaired, worst = 0, 0.0
    for (i, m, a1, a2, a3), index in solver_keys.items():
        partner = solver_keys.get((m, i, -a1, -a2, -a3))
        if partner is None:
            unpaired += 1
        else:
            worst = max(worst, abs(solver_value[index] - np.conj(solver_value[partner])))
    report("entries without the reverse entry", unpaired, 0, "")
    report("max |H_ij(R) - conj(H_ji(-R))|", worst, args.tol_matrix)

    common = [(solver_keys[key], enum_keys[key]) for key in solver_keys if key in enum_keys]
    si = np.array([c[0] for c in common])
    ee = np.array([c[1] for c in common])
    intra = abs(edisp[ee][:, 2]) < DZ_INTRA
    if args.structure_only:
        print("4. stored matrix elements by species pair and distance (no model given)")
        sp = t["species"]
        pairs = np.stack([np.minimum(sp[t["owner"]], sp[t["neighbor"]]), np.maximum(sp[t["owner"]], sp[t["neighbor"]])], 1)
        dist = np.round(np.linalg.norm(t["disp"], axis=1), 2)
        in_plane = abs(t["disp"][:, 2]) < DZ_INTRA
        print("     species  in-plane  distance (A)  entries   H mean (eV)    H min..max (eV)")
        shown = 0
        for a_, b_ in sorted(set(map(tuple, pairs))):
            for flag in (True, False):
                sel = (pairs[:, 0] == a_) & (pairs[:, 1] == b_) & (in_plane == flag)
                for r in np.unique(dist[sel])[: (99 if flag else 3)]:
                    m = sel & (dist == r)
                    v = solver_value[m].real
                    print(f"     {a_}-{b_}      {'yes' if flag else 'no ':3s}   {r:10.2f}   {m.sum():7d}   {v.mean():+11.5f}    {v.min():+.5f}..{v.max():+.5f}")
                    shown += 1
        report("largest |Im H| among the stored elements", abs(solver_value.imag).max())
        on = t["h0"] * g0
        for sp_ in np.unique(sp):
            print(f"     on-site energy, species {sp_}: {on[sp == sp_].min():+.5f}..{on[sp == sp_].max():+.5f} eV")
    else:
        print("4. matrix elements against the independent model")
        diff = abs(solver_value[si] - evalue[ee])
        intra = abs(edisp[ee][:, 2]) < DZ_INTRA
        if intra.any():
            report("intralayer entries: max |H_solver - H_model|", diff[intra].max(), args.tol_matrix)
        if (~intra).any():
            report("interlayer entries: max |H_solver - H_model|", diff[~intra].max(), args.tol_matrix)
        report("on-site energies: max |E_solver - E_model|", abs(t["h0"] * g0 - model_onsite).max(), args.tol_matrix)
        rho = np.round(np.hypot(edisp[ee][:, 0], edisp[ee][:, 1])[intra], 2)
        print("  intralayer shells:  distance (A)   neighbors/atom   H_solver (eV)   H_model (eV)")
        for radius in np.unique(rho):
            sel = np.flatnonzero(intra)[rho == radius]
            print(f"                      {radius:10.2f}   {len(sel) / n:14.2f}   {solver_value[si][sel].real.mean():+12.5f}   "
                  f"{evalue[ee][sel].mean():+11.5f}")
        if (~intra).any():
            print(f"  interlayer: {(~intra).sum() / n:.2f} neighbors/atom, largest |H| = "
                  f"{abs(solver_value[si][~intra]).max():.5f} eV, smallest |H| = {abs(solver_value[si][~intra]).min():.2e} eV")

    print("5. H(k): Hermiticity and agreement with the independent model")
    onsite = t["h0"] * g0
    corners = [[float(eval(x, {"__builtins__": {}})) for x in c.split(",")] + [0.0] for c in args.path.split(";")]
    rng = np.random.default_rng(args.seed)
    named = {"Gamma": [0, 0, 0], "K": [2 / 3, 1 / 3, 0], "M": [0.5, 0, 0],
             "generic (0.137, 0.291)": [0.137, 0.291, 0], "generic (-0.412, 0.073)": [-0.412, 0.073, 0]}
    for index in range(4):
        f = rng.uniform(-1, 1, 2)
        named[f"random {index + 1} ({f[0]:+.3f}, {f[1]:+.3f})"] = [f[0], f[1], 0]
    e_lattice = eimage @ ucell.T

    def solver_h(k):
        return bloch(n, t["owner"], t["neighbor"], solver_value, lattice, k, onsite)

    if args.structure_only:
        worst_herm = worst_g = 0.0
        k0 = rcell @ np.array([0.137, 0.291, 0.0])
        for label, f in named.items():
            hs = solver_h(rcell @ np.array(f, float))
            worst_herm = max(worst_herm, abs(hs - hs.conj().T).max())
        for g in ([1, 0, 0], [0, 1, 0], [-2, 3, 0]):
            worst_g = max(worst_g, abs(solver_h(k0 + rcell @ np.array(g, float)) - solver_h(k0)).max())
        report(f"largest Hermiticity residual of H(k) at {len(named)} k-points", worst_herm, args.tol_matrix)
        report("H(k+G) - H(k), element by element", worst_g, args.tol_matrix)
        x, bands = read_bands(f"{args.run_dir}/{args.prefix}.bands")
        worst_band = 0.0
        for kp, energies in zip(path_k_points(rcell, corners, x), bands):
            hs = solver_h(kp)
            worst_band = max(worst_band, abs(np.linalg.eigvalsh(0.5 * (hs + hs.conj().T)) - energies).max())
        report("eigenvalues of the stored tables against the band file", worst_band, args.tol_bands)
        if failures:
            print("MISMATCH: " + "; ".join(failures))
            return 1
        print("all residuals within their limits (structure only; no independent model)")
        return 0

    def model_h(k):
        return bloch(n, ei, em, evalue, e_lattice, k, model_onsite)

    worst_herm = worst_h = worst_e = 0.0
    print("     k-point                      |H-H^+|_max   |H_solver-H_model|_max   max |dE|")
    for label, f in named.items():
        k = rcell @ np.array(f, float)
        hs, hm = solver_h(k), model_h(k)
        herm = abs(hs - hs.conj().T).max()
        dh = abs(hs - hm).max()
        de = abs(np.linalg.eigvalsh(0.5 * (hs + hs.conj().T)) - np.linalg.eigvalsh(hm)).max()
        worst_herm, worst_h, worst_e = max(worst_herm, herm), max(worst_h, dh), max(worst_e, de)
        print(f"     {label:<28s} {herm:11.2e}   {dh:20.2e}   {de:9.2e}")
    report("largest Hermiticity residual", worst_herm, args.tol_matrix)
    report("largest |H_solver(k) - H_model(k)|", worst_h, args.tol_matrix)
    report("largest eigenvalue difference", worst_e, args.tol_matrix)

    print("6. invariances")
    k = rcell @ np.array([0.137, 0.291, 0.0])
    reference = np.linalg.eigvalsh(model_h(k))
    worst_g_h = worst_g_e = 0.0
    for g in ([1, 0, 0], [0, 1, 0], [-2, 3, 0]):
        kg = k + rcell @ np.array(g, float)
        worst_g_h = max(worst_g_h, abs(solver_h(kg) - solver_h(k)).max())
        worst_g_e = max(worst_g_e, abs(np.linalg.eigvalsh(model_h(kg)) - reference).max())
    report("H_solver(k+G) - H_solver(k), element by element", worst_g_h, args.tol_matrix)
    report("eigenvalues of the model at k+G", worst_g_e, args.tol_matrix)
    # atoms moved to other unit cells: new positions, new enumeration
    shift = rng.integers(-2, 3, size=(n, 2))
    moved = pos + np.hstack([shift, np.zeros((n, 1), int)]) @ ucell.T
    mi, mm, mimage, _, mvalue = model.enumerate(moved, ucell)
    hm_moved = bloch(n, mi, mm, mvalue, mimage @ ucell.T, k, model_onsite)
    gauge = np.exp(-1j * (moved - pos) @ k)
    report("atoms moved by lattice vectors: number of pairs changes by", abs(len(mi) - len(ei)), 0, "")
    report("  eigenvalues", abs(np.linalg.eigvalsh(hm_moved) - reference).max(), args.tol_matrix)
    report("  |H' - D^+ H D| with D = diag(exp(-i k.T))",
           abs(hm_moved - gauge.conj()[:, None] * model_h(k) * gauge[None, :]).max(), args.tol_matrix)
    # convention with the full displacement in the phase
    h_atomic = bloch(n, ei, em, evalue, edisp, k, model_onsite)
    report("phases exp(-i k.d) instead of exp(-i k.R): eigenvalues",
           abs(np.linalg.eigvalsh(h_atomic) - reference).max(), args.tol_matrix)

    print("7. eigenvalues of the independent model against the solver's band file")
    x, bands = read_bands(f"{args.run_dir}/{args.prefix}.bands")
    worst_band = 0.0
    for kp, energies in zip(path_k_points(rcell, corners, x), bands):
        worst_band = max(worst_band, abs(np.linalg.eigvalsh(model_h(kp)) - energies).max())
    report(f"max |E_model - E_solver| over {len(x)} k-points x {bands.shape[1]} bands", worst_band, args.tol_bands)

    if failures:
        print("MISMATCH: " + "; ".join(failures))
        return 1
    print("all residuals within their limits")
    return 0


if __name__ == "__main__":
    sys.exit(main())
