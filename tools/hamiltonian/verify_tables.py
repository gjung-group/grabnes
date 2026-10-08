#!/usr/bin/env python3
"""Independent check of the tight-binding Hamiltonian assembled by GRABNES.

Run a band calculation with ``WriteDataFiles .true.`` and point this script at
the run directory. From the tables the solver writes (neighbor list ``v``,
displacements ``dx``/``dy``/``dz``, hoppings ``<prefix>.s.mag``, positions,
cell, on-site energies) it

1. lists the intralayer neighbor shells and the matrix element on each;
2. checks the interlayer elements against the two-center formula
   H = Vpi exp(-(d-a_cc)/delta) (1 - (dz/d)^2) + Vsigma exp(-(d-d0)/delta) (dz/d)^2;
3. reports how complete and how symmetric the neighbor list is;
4. rebuilds H(k) from the tables, once with the lattice image of every
   neighbor taken from its displacement vector and once with the image rule of
   ``neigh.F90`` (matching by distance only), and compares both with the bands
   the solver wrote;
5. builds H(k) from scratch on the solver's atomic positions, with a complete
   and Hermitian neighbor set, and compares again.

Step 4 with the solver's rule validates the diagonalization and the reading of
the tables; step 5 measures the total effect of the neighbor-search
approximations. Requires NumPy. The band path is assumed to be the one of the
public examples (K - Gamma - M - K').
"""

import argparse
import collections
import re
import sys

import numpy as np


def load(run, prefix):
    lines = open(f"{run}/v").read().split("\n")
    count, neighbors = [], []
    i = 0
    while i < len(lines) and lines[i].strip():
        count.append(int(lines[i]))
        neighbors.append([int(x) - 1 for x in lines[i + 1].split()])
        i += 2

    def rows(name):
        return [[float(x) for x in line.split()] for line in open(f"{run}/{name}") if line.strip()]

    pair = re.compile(r"\(\s*([-+0-9.Ee]+)\s*,\s*([-+0-9.Ee]+)\s*\)")
    hop = [
        [complex(float(a), float(b)) for a, b in pair.findall(line)]
        for line in open(f"{run}/{prefix}.s.mag")
        if line.strip()
    ]
    disp = [np.array(list(zip(x, y, z))) for x, y, z in zip(rows("dx"), rows("dy"), rows("dz"))]
    cell = np.array(rows(f"{prefix}.cell"))
    onsite = rows(f"{prefix}.e")
    table = dict(
        count=count, neighbors=neighbors, hop=hop, disp=disp,
        pos=np.array(rows(f"{prefix}.pos")), ucell=cell[:3].T, rcell=cell[3:].T,
        h0=np.array([r[0] for r in onsite]), layer=np.array([int(r[2]) for r in onsite]),
    )
    assert all(len(h) == c == len(d) for h, c, d in zip(hop, count, disp)), "inconsistent tables"
    return table


def read_bands(path):
    lines = open(path).read().split("\n")
    bands, spins, points = map(int, lines[3].split())
    values = np.array([float(x) for line in lines[4:] for x in line.split()])
    values = values.reshape(points, bands * spins + 1)
    return values[:, 0], values[:, 1:]


def k_points(rcell, corners, path_length):
    corners = [rcell @ np.array(c) for c in corners]
    segments = [np.linalg.norm(b - a) for a, b in zip(corners[:-1], corners[1:])]
    result = []
    for x in path_length:
        i = 0
        while i < len(segments) - 1 and x > segments[i] + 1e-7:
            x -= segments[i]
            i += 1
        result.append(corners[i] + (corners[i + 1] - corners[i]) * min(x / segments[i], 1.0))
    return result


def lattice_images(t, rule):
    """Integer lattice vector of every neighbor entry, and the entries where the
    distance-only rule of neigh.F90 (last match in a -5..5 scan) picks another one."""
    ucell, pos = t["ucell"], t["pos"]
    images, different = [], 0
    scan = [np.array([ix, iy, 0]) for ix in range(-5, 6) for iy in range(-5, 6)]
    for i, (nbrs, disp) in enumerate(zip(t["neighbors"], t["disp"])):
        row = []
        for m, d in zip(nbrs, disp):
            true = np.round(np.linalg.solve(ucell, d - (pos[m] - pos[i]))).astype(int)
            chosen = true
            if rule == "solver":
                dist = np.linalg.norm(d)
                chosen = np.zeros(3, int)
                for n in scan:
                    if i == m and not n.any():
                        continue
                    if abs(np.linalg.norm(pos[i] - pos[m] - ucell @ n) - dist) < 0.1:
                        chosen = n
                different += not np.array_equal(chosen, true)
            row.append(chosen)
        images.append(row)
    return images, different


def h_from_tables(t, images, k, g0):
    n = len(t["count"])
    h = np.zeros((n, n), complex)
    for i in range(n):
        for m, hop, image in zip(t["neighbors"][i], t["hop"][i], images[i]):
            h[m, i] -= hop * np.exp(-1j * k @ (t["ucell"] @ image))
    lower = np.tril(h, -1)  # ZHEEV('N','L') only references the lower triangle
    return (lower + lower.conj().T + np.diag(t["h0"])) * g0


def h_independent(t, k, shells, two_center, cutoff):
    pos, ucell, layer = t["pos"], t["ucell"], t["layer"]
    h = np.zeros((len(pos), len(pos)), complex)
    same = layer[:, None] == layer[None, :]
    reach = int(np.ceil(cutoff / np.linalg.norm(ucell[:, 0]))) + 1
    for a in range(-reach, reach + 1):
        for b in range(-reach, reach + 1):
            shift = ucell @ np.array([a, b, 0])
            d = pos[None, :, :] + shift - pos[:, None, :]
            rho = np.hypot(d[..., 0], d[..., 1])
            dist = np.linalg.norm(d, axis=2)
            phase = np.exp(1j * k @ shift)
            for radius, energy in shells.items():
                h += np.where(same & (abs(rho - radius) < 0.02), energy, 0.0) * phase
            inter = (~same) & (rho < cutoff)
            cz = np.where(inter, (d[..., 2] / np.where(dist > 0, dist, 1.0)) ** 2, 0.0)
            v = (two_center["vpi"] * np.exp(-(dist - two_center["acc"]) / two_center["delta"]) * (1 - cz)
                 + two_center["vsigma"] * np.exp(-(dist - two_center["d0"]) / two_center["delta"]) * cz)
            h += np.where(inter, v, 0.0) * phase
    return h + np.diag(t["h0"])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_dir")
    parser.add_argument("--prefix", default="generate")
    parser.add_argument("--lattice-parameter", type=float, default=2.46)
    parser.add_argument("--g0", type=float, help="energy unit in eV (default: 12.14 - 3.72 a)")
    parser.add_argument("--vpppi0", type=float, default=3.5)
    parser.add_argument("--vppsigma0", type=float, default=0.48)
    parser.add_argument("--interlayer-distance", type=float, default=3.34)
    parser.add_argument("--layer-dist-factor", type=float, default=6.2)
    parser.add_argument("--max-table-dev", type=float,
                        help="fail if the tables (solver image rule) deviate more from the bands (eV)")
    parser.add_argument("--max-model-dev", type=float,
                        help="fail if the independent model deviates more from the bands (eV)")
    args = parser.parse_args()

    a = args.lattice_parameter
    acc = a / np.sqrt(3.0)
    g0 = args.g0 if args.g0 is not None else 12.14 - 3.72 * a
    t = load(args.run_dir, args.prefix)
    n = len(t["count"])
    print(f"atoms: {n}   energy unit g0 = {g0:.5f} eV   largest |on-site| = {abs(t['h0']).max() * g0:.3e} eV")

    shells = collections.defaultdict(list)
    inter = []
    for hops, disp in zip(t["hop"], t["disp"]):
        for hop, d in zip(hops, disp):
            if abs(d[2]) < 1e-6:
                shells[round(float(np.hypot(d[0], d[1])), 2)].append(-hop.real * g0)
            else:
                inter.append((np.hypot(d[0], d[1]), d[2], -hop.real * g0))
    print("intralayer shells:  distance (A)  neighbors/atom  H_ij (eV)")
    model_shells = {}
    for radius in sorted(shells):
        values = np.array(shells[radius])
        print(f"                    {radius:10.2f}  {len(values) / n:14.2f}  {values.mean():+.5f}"
              + ("" if np.ptp(values) < 1e-9 else f"  (spread {np.ptp(values):.1e})"))
        if abs(values.mean()) > 1e-12:
            model_shells[radius] = values.mean()

    two_center = dict(vpi=-args.vpppi0, vsigma=args.vppsigma0, acc=acc,
                      d0=args.interlayer_distance, delta=0.184 * a)
    cutoff = acc * 1.1 * args.layer_dist_factor
    if inter:
        inter = np.array(inter)
        dist = np.hypot(inter[:, 0], inter[:, 1])
        cz = (inter[:, 1] / dist) ** 2
        formula = (two_center["vpi"] * np.exp(-(dist - acc) / two_center["delta"]) * (1 - cz)
                   + two_center["vsigma"] * np.exp(-(dist - two_center["d0"]) / two_center["delta"]) * cz)
        print(f"interlayer: {len(inter) / n:.2f} neighbors/atom, in-plane distance up to {inter[:, 0].max():.2f} A, "
              f"largest |H_ij| = {abs(inter[:, 2]).max():.4f} eV")
        print(f"  deviation from the two-center formula: {abs(formula - inter[:, 2]).max():.2e} eV")
        pos, ucell, layer = t["pos"], t["ucell"], t["layer"]
        expected = 0
        reach = int(np.ceil(cutoff / np.linalg.norm(ucell[:, 0]))) + 1
        for a_ in range(-reach, reach + 1):
            for b_ in range(-reach, reach + 1):
                d = pos[None, :, :] + ucell @ np.array([a_, b_, 0]) - pos[:, None, :]
                expected += int(((layer[:, None] != layer[None, :]) & (np.hypot(d[..., 0], d[..., 1]) < cutoff)).sum())
        print(f"  pairs inside the nominal in-plane cutoff of {cutoff:.2f} A: {expected}; "
              f"in the solver list: {len(inter)} ({100 * len(inter) / expected:.1f} %)")

    true_images, _ = lattice_images(t, "true")
    solver_images, different = lattice_images(t, "solver")
    entries = {}
    for i in range(n):
        for m, image in zip(t["neighbors"][i], true_images[i]):
            entries[(i, m, image[0], image[1])] = True
    unpaired = sum((m, i, -a_, -b_) not in entries for (i, m, a_, b_) in entries)
    print(f"neighbor entries: {len(entries)}; without the reverse entry: {unpaired}; "
          f"given another lattice image by the distance-only rule: {different}")

    path, bands = read_bands(f"{args.run_dir}/{args.prefix}.bands")
    ks = k_points(t["rcell"], ([2 / 3, 1 / 3, 0], [0, 0, 0], [0.5, 0, 0], [2 / 3, 1 / 3, 0]), path)

    def deviation(build):
        return max(abs(np.linalg.eigvalsh(build(k)) - e).max() for k, e in zip(ks, bands))

    dev_solver = deviation(lambda k: h_from_tables(t, solver_images, k, g0))
    dev_true = deviation(lambda k: h_from_tables(t, true_images, k, g0))
    dev_model = deviation(lambda k: h_independent(t, k, model_shells, two_center, cutoff))
    print("largest band deviation from the solver output (eV):")
    print(f"  tables, lattice images as chosen by the solver : {dev_solver:.2e}")
    print(f"  tables, lattice images from the displacements  : {dev_true:.2e}")
    print(f"  independent model, complete neighbor set       : {dev_model:.2e}")

    status = 0
    if args.max_table_dev is not None and dev_solver > args.max_table_dev:
        print(f"MISMATCH: table reconstruction deviates by more than {args.max_table_dev:g} eV")
        status = 1
    if args.max_model_dev is not None and dev_model > args.max_model_dev:
        print(f"MISMATCH: independent model deviates by more than {args.max_model_dev:g} eV")
        status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
