#!/usr/bin/env python3
"""Known-answer checks of solver features that the examples do not exercise.

    python3 check_physics.py --bin /path/to/grabnes [--work DIR] [--launcher "mpirun -np 1"]

landau_levels   Graphene supercell (1152 atoms) with one flux quantum through the cell (137 T),
                nearest-neighbor hopping t. Analytic Dirac Landau levels: E_N = sgn(N) hbar v sqrt(2 e B |N| / hbar)
                with hbar v = 3 a_cc t / 2, so that with B = Phi_0 / A:  E_N = hbar v sqrt(4 pi |N| / A).
                Required: the N = 0 level at zero energy, twofold; E_1, E_2, E_3 within 0.5 %, 0.6 %, 0.8 % of
                the formula (the lattice corrections at this field are 0.16 %, 0.31 %, 0.47 %); the levels
                N = -3 ... 3 the same at three k-points, since the supercell is the magnetic unit cell.
hbn_gap         Monolayer hBN, nearest-neighbor hopping, 3 x 3 cell (the K points fold to Gamma). At K the
                hopping term vanishes, so the band edges are exactly the two on-site energies, each twice.
sparse_dense    Twisted bilayer (364 atoms): the 20 levels returned by the sparse solver (ARPACK, levels
                nearest to zero energy) against the dense diagonalization, in both directions.
tapw_two_way    The same bilayer with the plane-wave reduction (TAPW, 124 states per valley instead of 364):
                within 1 eV of charge neutrality the levels of the K and K' runs together must be the exact
                levels, one to one (equal counts, each within 0.05 meV), at three k-points. Control: with a
                wrong moire angle the same comparison must fail, so the test can tell a wrong basis.

input_messages  Messages of the input library: a key given twice is reported (the first value is used), and a
                key that nothing reads, here a misspelt one, is listed at the end of the run.

Exit status 1 if a check fails. Requires NumPy.
"""
import argparse
import math
import os
import shutil
import subprocess
import sys

import numpy as np


def unlimited_stack():
    """Run the solver with the largest stack the system allows: builds that keep automatic arrays on the
    stack (Intel Fortran by default) otherwise end with a segmentation fault for a few thousand atoms."""
    import resource
    hard = resource.getrlimit(resource.RLIMIT_STACK)[1]
    resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))

A_G = 2.46

BANDS = """Kubo.Calc .false.
Diag.Calc .true.
Calculate.Bands .true.
Bands.NumPoints 1
Bands.UseSameNumberOfPoints .true.
WriteDataFiles .false.
WritePos .false.
Neigh.fastNNnotsquare .true.
SuperCell 1
TB.NeighLevels 1
"""


def path(points):
    return f"&begin Bands.Path {len(points)}\n" + "\n".join(f"{a} {b} 0.0" for a, b in points) + "\n&end Bands.Path\n"


class Runner:
    def __init__(self, a):
        self.a, self.failures = a, 0

    def run(self, name, text, threads="2"):
        d = os.path.join(self.a.work, name)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        open(os.path.join(d, "Gendata.in"), "w").write("Prefix generate\n" + text)
        with open(os.path.join(d, "job.out"), "w") as out:
            rc = subprocess.run(self.a.launcher.split() + [os.path.abspath(self.a.bin), "Gendata.in"], cwd=d, stdout=out,
                                stderr=subprocess.STDOUT, preexec_fn=unlimited_stack, env=dict(os.environ, OMP_NUM_THREADS=threads, OMP_STACKSIZE="512M"), timeout=900).returncode
        if rc != 0:
            raise RuntimeError(f"the solver ended with status {rc} (see {d}/job.out)")
        # header: three lines, then "nAt nspin nk"; every k-point is its path length followed by the levels,
        # wrapped over several lines
        lines = open(os.path.join(d, "generate.bands")).read().split("\n")
        nk = int(lines[3].split()[2])
        values = np.array(" ".join(lines[4:]).split(), float)
        if nk < 1 or values.size % nk:
            raise RuntimeError(f"unexpected band file ({values.size} values for {nk} k-points) in {d}")
        return list(values.reshape(nk, -1)[:, 1:])

    def check(self, name, function):
        print(f"\n== {name}")
        try:
            problems = function()
        except Exception as exc:  # a failed run is a failed check
            problems = [str(exc)]
        for p in problems:
            print(f"  {p}")
        print("  PASS" if not problems else "  FAIL")
        self.failures += bool(problems)


def landau_levels(r):
    n = 24
    rows = r.run("landau_levels", f"TypeOfSystem Graphene\nCellSize {n}\nMagField.Integer 1\n" + BANDS
                 + path([(0.0, 0.0), (0.21, 0.37), (0.4, 0.13), (0.0, 0.0)]))
    g0 = 12.14 - 3.72 * A_G
    area = n * n * A_G * A_G * math.sqrt(3) / 2
    hv = 1.5 * A_G / math.sqrt(3) * g0
    problems = []
    ref = None
    for k, e in enumerate(rows[:3]):
        e = np.sort(e)
        pos = e[e > 1e-6]
        zero = e[abs(e) <= 1e-6]
        levels = [pos[0], pos[2], pos[4]]
        print(f"  k-point {k + 1}: N = 0 states {len(zero)}; E_1, E_2, E_3 = " + ", ".join(f"{x:.5f}" for x in levels) + " eV")
        if len(zero) != 2:
            problems.append(f"k-point {k + 1}: {len(zero)} states at zero energy instead of 2")
        for N, (x, tol) in enumerate(zip(levels, (0.005, 0.006, 0.008)), start=1):
            exact = hv * math.sqrt(4 * math.pi * N / area)
            if abs(x / exact - 1) > tol:
                problems.append(f"k-point {k + 1}: E_{N} = {x:.5f} eV, formula {exact:.5f} eV")
        if abs(pos[0] - pos[1]) > 2e-6 or abs(e[e < -1e-6][-1] + pos[0]) > 2e-6:
            problems.append(f"k-point {k + 1}: the N = 1 level is not twofold or not particle-hole symmetric")
        low = e[abs(e) < 0.8]          # Landau levels N = -3 ... 3; near the van Hove energy the bands are wide
        if ref is not None and (len(low) != len(ref) or np.max(abs(low - ref)) > 2e-6):
            problems.append(f"k-point {k + 1}: the levels below 0.8 eV depend on k")
        ref = low if ref is None else ref
    print("  formula: " + ", ".join(f"{hv * math.sqrt(4 * math.pi * N / area):.5f}" for N in (1, 2, 3)) + " eV")
    return problems


def hbn_gap(r):
    eb, en = 3.09, -1.89
    rows = r.run("hbn_gap", "TypeOfSystem BoronNitride\nCellSize 3\n" + BANDS + path([(0.0, 0.0), (0.5, 0.0), (0.0, 0.0)]))
    e = np.sort(rows[0])
    half = len(e) // 2
    edges = (e[half - 2], e[half - 1], e[half], e[half + 1])
    print("  band edges at the folded K points: " + ", ".join(f"{x:.5f}" for x in edges) + f" eV; gap {edges[2] - edges[1]:.5f} eV")
    want = (en, en, eb, eb)
    return [] if max(abs(a - b) for a, b in zip(edges, want)) < 2e-6 else [
        f"expected {want}: the on-site energies of N and B, each twice"]


def sparse_dense(r):
    common = ("TypeOfSystem TwistedBilayerBasedOnMoireCell\nTypeOfBL Koshino\nMoireCellParameters 6 5 5 6\ntwoLayers .true.\n"
              "twoLayersZ1 18.33\ntwoLayersZ2 21.67\nInterlayerDistance 3.34\nCellHeight 40.0\nNeigh.LayerNeighbors 100\n"
              "Neigh.LayerDistFactor 6.2\n" + BANDS + path([(0.21, 0.37), (0.0, 0.0), (0.5, 0.0), (0.21, 0.37)]))
    dense = r.run("sparse_dense_reference", common)
    sparse = r.run("sparse_dense", common + "sparseDiagSolver .true.\nBands.SparseNeig 20\n")
    problems = []
    for k in range(3):
        d = np.sort(dense[k])
        s = np.sort(sparse[k][sparse[k] != 0.0])
        nearest = np.sort(d[np.argsort(abs(d))[:len(s)]])
        worst = np.max(abs(s - nearest)) if len(s) == len(nearest) else float("inf")
        print(f"  k-point {k + 1}: {len(s)} sparse levels; largest difference from the {len(s)} dense levels nearest to zero: {worst:.2e} eV")
        if len(s) != 20 or worst > 5e-6:
            problems.append(f"k-point {k + 1}: the sparse levels are not the 20 dense levels nearest to zero")
    return problems


TBG = ("TypeOfSystem TwistedBilayerBasedOnMoireCell\nTypeOfBL Koshino\nMoireCellParameters 6 5 5 6\ntwoLayers .true.\n"
       "twoLayersZ1 18.33\ntwoLayersZ2 21.67\nInterlayerDistance 3.34\nCellHeight 40.0\nNeigh.LayerNeighbors 100\n"
       "Neigh.LayerDistFactor 6.2\n" + BANDS + path([(0.21, 0.37), (0.0, 0.0), (0.5, 0.0), (0.21, 0.37)]))


def tapw_two_way(r):
    exact = r.run("tapw_exact", TBG)

    def valleys(tag, angle):
        return [r.run(f"tapw_{tag}_{v}", TBG + "useTAPW .true.\nuseDenseMatrixTAPW .true.\nDiag.N_G 5\nTAPW.aG 2.46\n"
                      f"Diag.MoireAngle {angle}\nDiag.UseKprimeValley {flag}\n", threads="1")
                for v, flag in (("K", ".false."), ("Kp", ".true."))]

    def compare(runs, label):
        worst, equal_counts = 0.0, True
        for k in range(3):
            e = np.sort(exact[k])
            e0 = 0.5 * (e[len(e) // 2 - 1] + e[len(e) // 2])
            t = np.sort(np.concatenate([x[k] for x in runs]))
            t, ew = t[(t != 0.0) & (abs(t - e0) < 1.0)], e[abs(e - e0) < 1.0]
            same = len(t) == len(ew) and len(t) > 0
            d = float(np.max(abs(t - ew))) if same else float("inf")
            print(f"  {label}, k-point {k + 1}: TAPW levels {len(t)}, exact levels {len(ew)}, largest difference "
                  + (f"{d * 1e3:.4f} meV" if same else "undefined"))
            equal_counts &= same
            worst = max(worst, d)
        return equal_counts, worst

    problems = []
    ok, worst = compare(valleys("angle0", "0.0"), "moire angle 0")
    if not ok or worst > 5e-5:
        problems.append("the TAPW levels of the two valleys are not the exact levels within 1 eV of neutrality")
    ok, worst = compare(valleys("control", "26.9955"), "control (wrong angle)")
    if ok and worst < 5e-3:
        problems.append("the control with a wrong moire angle also agrees: the comparison is not sensitive")
    return problems


def input_messages(r):
    rows = r.run("input_messages", "TypeOfSystem Graphene\nCellSize 3\nTB.NeighLevel 3\nCellSize 6\n" + BANDS
                 + path([(0.0, 0.0), (0.5, 0.0), (0.0, 0.0)]))
    log = open(os.path.join(r.a.work, "input_messages", "job.out")).read()
    problems = []
    if len(rows[0]) != 18:
        problems.append(f"{len(rows[0])} levels instead of 18: the first value of the repeated key CellSize was not used")
    if 'the key "CellSize" appears 2 times' not in log:
        problems.append("the repeated key CellSize is not reported")
    listed = log.split("were not read in this run")[1].split("....")[0].split() if "were not read in this run" in log else []
    print("  keys reported as not read: " + (", ".join(x for x in listed if x not in ("input:", "(misspelt,", "or", "not", "used", "by", "the", "selected", "calculation):")) or "none"))
    if "TB.NeighLevel" not in listed:
        problems.append("the misspelt key TB.NeighLevel is not listed as not read")
    if "TB.NeighLevels" in listed or "CellSize" in listed:
        problems.append("a key that the run reads is listed as not read")
    return problems


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--bin", required=True)
    p.add_argument("--work", default="check_physics_work")
    p.add_argument("--launcher", default="")
    a = p.parse_args()
    r = Runner(a)
    r.check("landau_levels", lambda: landau_levels(r))
    r.check("hbn_gap", lambda: hbn_gap(r))
    r.check("sparse_dense", lambda: sparse_dense(r))
    r.check("tapw_two_way", lambda: tapw_two_way(r))
    r.check("input_messages", lambda: input_messages(r))
    print()
    print("PASS: all physics checks" if not r.failures else f"FAIL: {r.failures} physics check(s)")
    sys.exit(1 if r.failures else 0)


if __name__ == "__main__":
    main()
