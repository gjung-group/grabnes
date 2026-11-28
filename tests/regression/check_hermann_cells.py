#!/usr/bin/env python3
"""Structures built from Hermann's four indices (Structure.MoireCellParameters m1 m2 m3 m4): independent checks.

For every case the solver builds the cell, writes its tables (WriteDataFiles; ten recursion steps follow) and the
script checks, from generate.pos, generate.cell, the neighbour tables and job.out alone:

  atoms     the number of atoms of every layer is 2 N(pair) (N = m^2 + m n + n^2 of the pair the layer is built from)
  cell      the cell vectors are A1 = aG (p1 a1 + p2 a2) of the pair built first and A2 = A1 turned by +60 deg
  lattice   every atom of a layer sits on the honeycomb lattice that the indices define (layer 2 = layer 1 scaled by
            |alpha/beta| and turned by arg(alpha/beta), alpha = m1 + m2 tau, beta = m3 + m4 tau, tau = exp(i pi/3)),
            on the sublattice of its Species
  cosets    the atoms of each sublattice of each layer are a bijection onto the lattice modulo the superlattice:
            N(pair) distinct classes, none missing, no duplicate position
  periodic  A1 and A2 are lattice vectors of every layer, and every atom has exactly three in-layer neighbours at the
            bond length of its layer in the solver's own neighbour list (periodic boundary conditions)
  csl       the shortest common vector of the two lattices, found by brute force from the positions' lattices (no
            index arithmetic): |A1| when the solver says "primitive", |A1|/sqrt(n) when it says "holds n primitive
            cells"; and the printed signed angle and primitive indices agree with tools/structure/hermann.py
  same      index sets that describe the same cell (unit multiples) give the same set of positions modulo the cell

Refusals (zero pairs, TwistedBilayer with BasedOnMoireCellParameters) must stop with the solver's message.

    python3 check_hermann_cells.py --bin /path/to/grabnes [--work DIR] [--launcher "mpirun -np 1"] [--only NAME]

With --tbg-sweep NATOMS (not part of the regression suite: hours) the cases are instead every primitive commensurate
twisted bilayer with at most NATOMS atoms and an angle in (0, 60) deg, one cell per angle from the signed-index search
of tools/structure/hermann.py; each is checked as above and against Mele (2010) and Lopes dos Santos et al. (2012):
the same set of angles (two-way), the same atom count and area for every angle. --jobs runs several cases at once.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools", "structure"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools", "hamiltonian"))
import hermann                                   # noqa: E402
from verify_tables import load                   # noqa: E402

AG = 2.46
TAU = complex(0.5, np.sqrt(3.0) / 2.0)
RUN = """Prefix generate
Kubo.Calc .true.
Diag.Calc .false.
Calculate.OnlyDOS .true.
Calculate.Bands .false.
RecursionNumber 10
NumberofEnergyPoints 11
Epsilon 0.05
EnergyMin 0-9.0
EnergyMax 9.0
setSeed .true.
seedValue 123456
WriteDataFiles .true.
WritePos .false.
SuperCell 1
Neigh.fastNNnotsquare .true.
TB.NeighLevels 1
CellHeight 40.0
InterlayerDistance 3.35
"""
BILAYER = "TypeOfBL Koshino\ntwoLayers .true.\ntwoLayersZ1 18.325\ntwoLayersZ2 21.675\nNeigh.LayerNeighbors 60\n" \
          "Neigh.LayerDistFactor 2.0\n"

# name: (TypeOfSystem, indices, extra input, layers from the bottom: 'f' = pair built first, 's' = second)
CASES = {
    # twisted bilayer graphene: the table of the task, its primitive cells in both handednesses, signed indices
    "tbg_2_1_1_2": ("TwistedBilayerBasedOnMoireCell", (2, 1, 1, 2), BILAYER, "sf"),
    "tbg_4_1_1_4": ("TwistedBilayerBasedOnMoireCell", (4, 1, 1, 4), BILAYER, "sf"),
    "tbg_3_m1_2_1": ("TwistedBilayerBasedOnMoireCell", (3, -1, 2, 1), BILAYER, "sf"),
    "tbg_2_1_3_m1": ("TwistedBilayerBasedOnMoireCell", (2, 1, 3, -1), BILAYER, "sf"),
    "tbg_5_2_2_5": ("TwistedBilayerBasedOnMoireCell", (5, 2, 2, 5), BILAYER, "sf"),
    "tbg_4_m1_3_1": ("TwistedBilayerBasedOnMoireCell", (4, -1, 3, 1), BILAYER, "sf"),
    "tbg_3_1_4_m1": ("TwistedBilayerBasedOnMoireCell", (3, 1, 4, -1), BILAYER, "sf"),
    "tbg_7_1_1_7": ("TwistedBilayerBasedOnMoireCell", (7, 1, 1, 7), BILAYER, "sf"),
    "tbg_5_m2_3_2": ("TwistedBilayerBasedOnMoireCell", (5, -2, 3, 2), BILAYER, "sf"),
    "tbg_7_4_4_7": ("TwistedBilayerBasedOnMoireCell", (7, 4, 4, 7), BILAYER, "sf"),
    "tbg_6_m1_5_1": ("TwistedBilayerBasedOnMoireCell", (6, -1, 5, 1), BILAYER, "sf"),
    # the primitive indices job.out prints for 4 1 1 4 and 7 1 1 7, written with a bare minus sign (*_bare)
    "tbg_1_2_m1_3_bare": ("TwistedBilayerBasedOnMoireCell", (1, 2, -1, 3), BILAYER, "sf"),
    "tbg_2_3_m2_5_bare": ("TwistedBilayerBasedOnMoireCell", (2, 3, -2, 5), BILAYER, "sf"),
    # equivalent representations: unit multiples (same cell), layer exchange and mirror (opposite twist)
    "tbg_m2_m1_m1_m2": ("TwistedBilayerBasedOnMoireCell", (-2, -1, -1, -2), BILAYER, "sf"),
    "tbg_m1_3_m2_3": ("TwistedBilayerBasedOnMoireCell", (-1, 3, -2, 3), BILAYER, "sf"),
    "tbg_1_4_4_1": ("TwistedBilayerBasedOnMoireCell", (1, 4, 4, 1), BILAYER, "sf"),
    "tbg_5_m1_5_m4": ("TwistedBilayerBasedOnMoireCell", (5, -1, 5, -4), BILAYER, "sf"),
    "tbg_2_1_m2_3": ("TwistedBilayerBasedOnMoireCell", (2, 1, -2, 3), BILAYER, "sf"),     # layer 2 turned 60 deg more
    # zero indices, non-coprime indices, a common Eisenstein factor without a common integer factor
    "tbg_5_0_0_5": ("TwistedBilayerBasedOnMoireCell", (5, 0, 0, 5), BILAYER, "sf"),
    "tbg_1_0_1_0": ("TwistedBilayerBasedOnMoireCell", (1, 0, 1, 0), BILAYER, "sf"),
    "tbg_4_2_2_4": ("TwistedBilayerBasedOnMoireCell", (4, 2, 2, 4), BILAYER, "sf"),
    "tbg_5_6_m1_10": ("TwistedBilayerBasedOnMoireCell", (5, 6, -1, 10), BILAYER, "sf"),
    # first cell vector along y (the fractional coordinates divided by its x component)
    "tbg_1_1_1_m2": ("TwistedBilayerBasedOnMoireCell", (1, 1, 1, -2), BILAYER, "sf"),
    "tbg_11_2_7_m14": ("TwistedBilayerBasedOnMoireCell", (11, 2, 7, -14), BILAYER, "sf"),
    # a larger cell, magic-angle neighbourhood
    "tbg_17_16_16_17": ("TwistedBilayerBasedOnMoireCell", (17, 16, 16, 17), BILAYER, "sf"),
    # lattice-mismatched pairs (Graphene_Over_BN: the BN lattice constant follows from the indices)
    "gbn_5_0_4_0": ("Graphene_Over_BN", (5, 0, 4, 0), BILAYER, "sf"),
    "gbn_10_0_8_0": ("Graphene_Over_BN", (10, 0, 8, 0), BILAYER, "sf"),
    "gbn_56_0_55_0": ("Graphene_Over_BN", (56, 0, 55, 0), BILAYER, "sf"),
    "gbn_4_0_m3_3": ("Graphene_Over_BN", (4, 0, -3, 3), BILAYER, "sf"),
    "gbn_9_m4_4_3": ("Graphene_Over_BN", (9, -4, 4, 3), BILAYER, "sf"),
    "gbn_3_0_2_m4": ("Graphene_Over_BN", (3, 0, 2, -4), BILAYER, "sf"),
    # graphene on hBN (aBN = 2.504 Ang to 0.2 %) twisted by 0.6 and 5 deg: primitive cells of hermann.search_hermann
    "gbn_26_26_26_25": ("Graphene_Over_BN", (26, 26, 26, 25), BILAYER, "sf"),
    "gbn_m2_31_1_29": ("Graphene_Over_BN", (-2, 31, 1, 29), BILAYER, "sf"),
    "bgbn_5_0_4_0": ("BilayerGraphene_Over_BN", (5, 0, 4, 0), BILAYER, "sff"),
    "bgbn_9_m4_4_3": ("BilayerGraphene_Over_BN", (9, -4, 4, 3), BILAYER, "sff"),
    "encap_3_m1_2_1": ("MoireEncapsulatedBilayerBasedOnMoireCell", (3, -1, 2, 1), BILAYER, "sf"),
    "graphene_3_m1_2_1": ("Graphene", (3, -1, 2, 1), "BasedOnMoireCellParameters .true.\n", "f"),
    "graphene_5_0_0_5": ("Graphene", (5, 0, 0, 5), "BasedOnMoireCellParameters .true.\n", "f"),
    # refused
    "refuse_zero_first": ("TwistedBilayerBasedOnMoireCell", (0, 0, 1, 2), BILAYER, None),
    "refuse_zero_second": ("TwistedBilayerBasedOnMoireCell", (1, 2, 0, 0), BILAYER, None),
    "refuse_twistedbilayer_basedon": ("TwistedBilayer", (2, 1, 1, 2),
                                      BILAYER + "BasedOnMoireCellParameters .true.\nCellSize 6\n", None),
}
# index sets that describe the same structure (same positions modulo the cell)
SAME = [("tbg_2_1_1_2", "tbg_m2_m1_m1_m2"), ("tbg_2_1_1_2", "tbg_m1_3_m2_3"),
        # theta + 60 deg: the layers turn about a hexagon centre, which maps each honeycomb layer onto itself
        ("tbg_2_1_1_2", "tbg_2_1_m2_3")]
# (x, y): y is the mirror image (y -> -y) of x: theta and 60 deg - theta (21.79 / 38.21 ...) with this stacking, the
# layer exchange and the conjugate indices (m1+m2, -m2, m3+m4, -m4)
MIRROR = [("tbg_2_1_1_2", "tbg_3_m1_2_1"), ("tbg_3_1_4_m1", "tbg_4_m1_3_1"), ("tbg_4_1_1_4", "tbg_1_4_4_1"),
          ("tbg_4_1_1_4", "tbg_5_m1_5_m4")]
# (primitive, non-primitive, n): the non-primitive structure is n copies of the primitive one (positions modulo the
# primitive cell); the primitive cells are the same-handedness ones ((2a+b)/3, (b-a)/3, (a+2b)/3, (a-b)/3)
COPIES = [("tbg_3_m1_2_1", "tbg_4_1_1_4", 3), ("tbg_4_m1_3_1", "tbg_5_2_2_5", 3), ("tbg_5_m2_3_2", "tbg_7_1_1_7", 3),
          ("tbg_6_m1_5_1", "tbg_7_4_4_7", 3), ("tbg_1_2_m1_3_bare", "tbg_4_1_1_4", 3),
          ("tbg_2_3_m2_5_bare", "tbg_7_1_1_7", 3), ("tbg_2_1_1_2", "tbg_4_2_2_4", 4), ("gbn_5_0_4_0", "gbn_10_0_8_0", 4)]
REFUSAL = {"refuse_zero_first": "must both be nonzero", "refuse_zero_second": "must both be nonzero",
           "refuse_twistedbilayer_basedon": "not implemented for TypeOfSystem TwistedBilayer"}


def int_in(i):
    return f"0-{-i}" if i < 0 else str(i)


def run(a, name, typ, m, extra):
    d = os.path.join(a.work, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    fmt = str if name.endswith("_bare") else int_in     # -1, or the 0-1 form older binaries need
    text = RUN + f"TypeOfSystem {typ}\nMoireCellParameters {' '.join(fmt(x) for x in m)}\n" + extra
    open(os.path.join(d, "Gendata.in"), "w").write(text)
    with open(os.path.join(d, "job.out"), "w") as out:
        rc = subprocess.run(a.launcher.split() + [os.path.abspath(a.bin), "Gendata.in"], cwd=d, stdout=out,
                            stderr=subprocess.STDOUT, timeout=1800,
                            env=dict(os.environ, OMP_NUM_THREADS="2", OMP_STACKSIZE="512M")).returncode
    return d, rc, open(os.path.join(d, "job.out")).read()


def order(m):
    """The pair built first (layer 1) and second, as AtomsPos orders them: the pairs are exchanged when
    N(m3, m4) >= N(m1, m2)."""
    al, be = (m[0], m[1]), (m[2], m[3])
    return (be, al) if hermann.norm(be) >= hermann.norm(al) else (al, be)


def cz(x):
    return x[0] + x[1] * TAU


def lattice_basis(z):
    """2x2 Cartesian basis (columns a1, a2) of the honeycomb Bravais lattice z * aG * Z[tau]."""
    return np.array([[(z * AG).real, (z * AG * TAU).real], [(z * AG).imag, (z * AG * TAU).imag]])


def coset_keys(n, pair):
    """Classes of the integer lattice points n (rows) modulo the superlattice of pair (columns (p1,p2), (-p2,p1+p2))."""
    p1, p2 = pair
    m = np.array([[p1, -p2], [p2, p1 + p2]], float)
    det = int(round(np.linalg.det(m)))
    f = np.linalg.solve(m, n.T).T * det                       # integer rows: det * M^-1 n
    f = np.rint(f).astype(np.int64) % abs(det)
    return set(map(tuple, f)), len(f)


def csl_shortest(z2, rmax):
    """Shortest nonzero vector common to the lattices aG Z[tau] and z2 aG Z[tau], by brute force over layer 1."""
    b1, b2 = lattice_basis(1.0), lattice_basis(z2)
    inv2 = np.linalg.inv(b2)
    nmax = int(rmax / AG * 2 / np.sqrt(3)) + 2
    i, j = np.meshgrid(np.arange(-nmax, nmax + 1), np.arange(-nmax, nmax + 1))
    v = b1 @ np.vstack([i.ravel(), j.ravel()])
    r = np.hypot(v[0], v[1])
    keep = (r > 1e-9) & (r <= rmax * (1 + 1e-9))
    v, r = v[:, keep], r[keep]
    s = inv2 @ v
    common = np.all(np.abs(s - np.rint(s)) < 1e-6, axis=0)
    return r[common].min() if common.any() else None


def joined_log(log):
    """job.out as one line: module prefixes and the line wrapping of MIO_Print removed."""
    return " ".join(" ".join(re.sub(r"^\s*[a-z]+:\s?", "", l) for l in log.split("\n")).split())


def check_case(a, name):
    typ, m, extra, layers = CASES[name]
    d, rc, log = run(a, name, typ, m, extra)
    joined = joined_log(log)
    problems = []
    if layers is None:
        if rc == 0 or REFUSAL[name] not in joined:
            problems.append(f"expected a refusal containing '{REFUSAL[name]}' (status {rc})")
        else:
            print(f"  refused: {REFUSAL[name]}")
        return problems, None
    if rc != 0:
        return [f"the solver ended with status {rc} (see {d}/job.out)"], None
    t = load(d, "generate")
    pos, ucell, species = t["pos"], t["ucell"], np.array(t["species"])
    first, second = order(m)
    zf = 1.0
    zs = cz(first) / cz(second)                          # lattice of the pair built second, relative to layer 1
    want = {"f": (first, zf), "s": (second, zs)}
    # cell
    a1, a2 = AG * cz(first), AG * cz(first) * TAU
    cell_err = max(abs(ucell[0, 0] - a1.real), abs(ucell[1, 0] - a1.imag), abs(ucell[0, 1] - a2.real),
                   abs(ucell[1, 1] - a2.imag))
    if cell_err > 1e-6:
        problems.append(f"cell vectors differ from aG (first pair) by {cell_err:.2e} Ang")
    # layers, from the bottom
    zval = np.round(pos[:, 2], 3)
    levels = sorted(set(zval))
    if len(levels) != len(layers):
        return problems + [f"{len(levels)} atomic planes, expected {len(layers)}"], None
    keys_all = {}
    for lev, tag in zip(levels, layers):
        sel = np.where(zval == lev)[0]
        pair, z = want[tag]
        nexp = 2 * hermann.norm(pair)
        if len(sel) != nexp:
            problems.append(f"plane z={lev}: {len(sel)} atoms, expected 2 N{pair} = {nexp}")
            continue
        s = np.linalg.solve(lattice_basis(z), pos[sel, :2].T).T
        sub = []
        bad = 0
        for row in s:
            for k, b in enumerate((1 / 3, 2 / 3)):
                if np.all(np.abs((row - b) - np.rint(row - b)) < 1e-6):
                    sub.append(k)
                    break
            else:
                bad += 1
                sub.append(-1)
        sub = np.array(sub)
        if bad:
            problems.append(f"plane z={lev}: {bad} atoms off the honeycomb lattice of {pair}")
            continue
        sp = species[sel]
        if not (np.all((sp[sub == 0] % 2) == 1) and np.all((sp[sub == 1] % 2) == 0)):
            problems.append(f"plane z={lev}: Species do not follow the sublattice")
        for k, b in enumerate((1 / 3, 2 / 3)):
            n = np.rint(s[sub == k] - b)
            keys, cnt = coset_keys(n, pair)
            if cnt != hermann.norm(pair) or len(keys) != cnt:
                problems.append(f"plane z={lev} sublattice {k}: {cnt} atoms, {len(keys)} distinct classes "
                                f"modulo the cell, expected {hermann.norm(pair)}")
        # A1, A2 lattice vectors of this layer
        sa = np.linalg.solve(lattice_basis(z), np.array([[a1.real, a2.real], [a1.imag, a2.imag]]))
        if np.abs(sa - np.rint(sa)).max() > 1e-6:
            problems.append(f"plane z={lev}: the cell vectors are not lattice vectors of the layer")
        keys_all[lev] = pos[sel]
    # in-layer coordination from the solver's neighbour list (periodic boundary conditions)
    dz = np.abs(t["disp"][:, 2])
    dl = np.hypot(t["disp"][:, 0], t["disp"][:, 1])
    bond = {lev: AG * abs(want[tag][1]) / np.sqrt(3) for lev, tag in zip(levels, layers)}
    bl = np.array([bond[zval[i]] for i in t["owner"]])
    intra = (dz < 0.5) & (np.abs(dl - bl) < 1e-4)
    count = np.bincount(t["owner"][intra], minlength=len(pos))
    if np.any(count != 3):
        problems.append(f"{np.sum(count != 3)} atoms without exactly three in-layer neighbours at the bond length")
    # coincidence lattice by brute force, against the solver's verdict
    if "s" in layers and "f" in layers:
        lmin = csl_shortest(zs, abs(a1) + 1e-6)
        d_ = hermann.describe(m)
        ratio = round((abs(a1) / lmin) ** 2) if lmin else 0
        said = 1 if "Hermann cell: primitive" in joined else None
        mm = re.search(r"it holds (\d+) primitive cells", joined)
        if mm:
            said = int(mm.group(1))
        if said != ratio or ratio != d_["mult"]:
            problems.append(f"coincidence lattice: shortest common vector gives {ratio} primitive cells, the solver "
                            f"says {said}, hermann.py {d_['mult']}")
        th = re.search(r"Hermann cell: the \(m3, m4\) lattice is the \(m1, m2\) lattice turned by\s+(\S+) deg", joined)
        if not th or abs(float(th.group(1)) - d_["theta"]) > 1e-5:
            problems.append(f"printed signed angle {th.group(1) if th else None}, expected {d_['theta']:.6f}")
        if d_["mult"] > 1:
            pm = re.search(r"Primitive cell with the same twist:\s+(-?\d+) (-?\d+) (-?\d+) (-?\d+)", joined)
            if not pm or not hermann.same_cell(tuple(int(x) for x in pm.groups()), d_["prim"]):
                problems.append("printed primitive indices missing or not the primitive cell")
        print(f"  {typ} {m}: {len(pos)} atoms, (m3, m4) lattice turned by {d_['theta']:+.4f} deg, |A1| {abs(a1):.4f} Ang, "
              f"{'primitive' if ratio == 1 else f'{ratio} primitive cells'} (brute-force coincidence lattice)")
    else:
        print(f"  {typ} {m}: {len(pos)} atoms, one layer")
    return problems, (pos, ucell)


def folded(x, cell):
    """Positions folded into the cell (fractional, rounded), with z, as a list of keys."""
    f = (np.linalg.inv(cell[:2, :2]) @ x[:, :2].T).T
    f = f - np.floor(f + 1e-7)
    f[np.abs(f - 1) < 1e-6] = 0.0
    f[np.abs(f) < 1e-6] = 0.0
    return list(map(tuple, np.round(np.column_stack([f, x[:, 2]]), 5)))


def copies_of(p, q, n):
    """q is n copies of p: q folded into the cell of p hits every position of p exactly n times."""
    from collections import Counter
    (pp, up), (qq, _) = p, q
    kp, kq = Counter(folded(pp, up)), Counter(folded(qq, up))
    return len(qq) == n * len(pp) and set(kp) == set(kq) and all(c == 1 for c in kp.values()) \
        and all(c == n for c in kq.values())


def tbg_sweep(a):
    """Every primitive TBG cell up to a.tbg_sweep atoms through the solver, against Mele and Lopes dos Santos."""
    from concurrent.futures import ThreadPoolExecutor
    table = hermann.tbg_table(a.tbg_sweep)
    mele, lopes = hermann.mele_angles(a.tbg_sweep), hermann.lopes_angles(a.tbg_sweep)
    names = []
    for th, (nat, m) in sorted(table.items()):
        name = "sweep_" + "_".join(int_in(x).replace("0-", "m") for x in m)
        CASES[name] = ("TwistedBilayerBasedOnMoireCell", m, BILAYER, "sf")
        names.append((name, th, nat, m))

    def one(item):
        name, th, nat, m = item
        try:
            problems, _ = check_case(a, name)
        except Exception as exc:
            problems = [f"{type(exc).__name__}: {exc}"]
        if not a.keep:
            shutil.rmtree(os.path.join(a.work, name), ignore_errors=True)
        return name, th, nat, m, problems
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        results = list(ex.map(one, names))
    failures = 0
    for name, th, nat, m, problems in results:
        if problems:
            failures += 1
            print(f"FAIL {name} {th:.6f} deg: " + "; ".join(problems))
    mk = {round(k, 7): v for k, v in mele.items()}
    lk = {round(k, 7): v for k, v in lopes.items()}
    hk = {round(th, 7): nat for _, th, nat, _, _ in results}
    two_way = set(hk) == set(mk) == set(lk)
    counts = all(hk[k] == mk[k] == lk[k] for k in hk) if two_way else False
    print(f"\n{len(results)} cells built and checked, {failures} failed")
    print(f"angles: Hermann {len(hk)}, Mele {len(mk)}, Lopes dos Santos {len(lk)}; same set: {two_way}; "
          f"same atom counts (hence areas N sqrt(3)/2 aG^2): {counts}")
    return failures + (not two_way) + (not counts)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--bin", required=True)
    p.add_argument("--work", default="check_hermann_work")
    p.add_argument("--launcher", default="")
    p.add_argument("--only", default=None)
    p.add_argument("--tbg-sweep", type=int, default=0, metavar="NATOMS")
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--keep", action="store_true", help="--tbg-sweep: keep the run directories")
    a = p.parse_args()
    if a.tbg_sweep:
        n = tbg_sweep(a)
        print("PASS: TBG sweep" if not n else f"FAIL: TBG sweep ({n})")
        sys.exit(1 if n else 0)
    failures, built = 0, {}
    for name in CASES:
        if a.only and a.only not in name:
            continue
        print(f"\n== hermann {name}")
        try:
            problems, built[name] = check_case(a, name)
        except Exception as exc:                      # a failed run is a failed check
            problems = [f"{type(exc).__name__}: {exc}"]
        for x in problems:
            print(f"  {x}")
        print("  PASS" if not problems else "  FAIL")
        failures += bool(problems)
    for x, y, n in [(x, y, 1) for x, y in SAME] + COPIES:
        if built.get(x) is None or built.get(y) is None:
            continue
        ok = copies_of(built[x], built[y], n)
        what = "same structure" if n == 1 else f"{n} copies of the primitive structure"
        print(f"\n== hermann {y} = {what} {x}\n  {'PASS' if ok else 'FAIL'}")
        failures += not ok
    for x, y in MIRROR:
        if built.get(x) is None or built.get(y) is None:
            continue
        px, ux = built[x]
        ok = copies_of(built[y], (px * np.array([1.0, -1.0, 1.0]), ux), 1)
        print(f"\n== hermann {y} = mirror image of {x} (stacking about a hexagon centre)\n  {'PASS' if ok else 'FAIL'}")
        failures += not ok
    print()
    print("PASS: all Hermann cell checks" if not failures else f"FAIL: {failures} Hermann cell check(s)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
