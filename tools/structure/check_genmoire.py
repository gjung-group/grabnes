#!/usr/bin/env python3
"""Check of the cell search of usefulGeneralScripts/twistedBilayers/GenMoire.f90 against hermann.py.

Compiles GenMoire.f90 (gfortran, or --fc) and, for every commensurate twisted-bilayer angle in (0, 60) deg with at
most NATOMS atoms in its primitive cell (default: the five angles of the task's table; --all: every angle, 828 up to
12000 atoms), asks it for that angle with tight tolerances. Pass: the cell it returns has that angle and the
primitive atom count of hermann.py / Mele, its "Smallest" line too, and `GenMoire cell` reports the given cell's
primitive cell like hermann.describe.

    python3 check_genmoire.py [--all] [--natoms 12000] [--work DIR] [--fc gfortran]
"""

import argparse
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hermann  # noqa: E402

SRC = os.path.join(HERE, "..", "..", "usefulGeneralScripts", "twistedBilayers", "GenMoire.f90")
ROW = re.compile(r"^(For angle|For delta|For both|Smallest):\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(\S+)\s+(\S+)",
                 re.M)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--all", action="store_true")
    p.add_argument("--natoms", type=int, default=12000)
    p.add_argument("--work", default="check_genmoire_work")
    p.add_argument("--fc", default="gfortran")
    a = p.parse_args()
    os.makedirs(a.work, exist_ok=True)
    exe = os.path.join(a.work, "GenMoire")
    subprocess.run([a.fc, "-O2", "-o", exe, SRC], check=True)
    table = hermann.tbg_table(a.natoms)
    if not a.all:
        want = [21.7868, 38.2132, 27.7958, 46.8264, 17.8966]
        table = {k: v for k, v in table.items() if any(abs(k - w) < 1e-4 for w in want)}
    fails = 0
    for th, (nat, m) in sorted(table.items()):
        out = subprocess.run([os.path.abspath(exe), "search", f"{th:.9f}", "2.46", "1e-5", "1e-6",
                              str(nat // 4), "0", "1000000000"], cwd=a.work, capture_output=True, text=True).stdout
        rows = {r[0]: r[1:] for r in ROW.findall(out)}
        ok = True
        for key in ("For angle", "Smallest"):
            if key not in rows:
                ok = False
                continue
            c = tuple(int(x) for x in rows[key][:4])
            d = hermann.describe(c)
            if not (d["primitive"] and 4 * d["n1"] == nat and abs(abs(d["theta"]) - th) < 1e-6):
                ok = False
        # exact-geometry mode on the non-primitive double of the cell
        big = hermann.mul((2, 0), (m[0], m[1])) + hermann.mul((2, 0), (m[2], m[3]))
        out2 = subprocess.run([os.path.abspath(exe), "cell", *map(str, big)], cwd=a.work, capture_output=True,
                              text=True).stdout
        mm = re.search(r"holds\s+(\d+) primitive cells.*?Primitive cell \(same twist and lattice ratio\):\s+(-?\d+)\s+"
                       r"(-?\d+)\s+(-?\d+)\s+(-?\d+)", out2, re.S)
        if not mm or int(mm.group(1)) != 4 or not hermann.same_cell(tuple(int(x) for x in mm.groups()[1:]), m):
            ok = False
        if not ok:
            fails += 1
            print(f"FAIL {th:.7f} deg ({nat} atoms): {rows}")
    print(f"GenMoire: {len(table) - fails} of {len(table)} angles found with their primitive cell")
    print("PASS" if not fails else "FAIL")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
