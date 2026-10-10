#!/usr/bin/env python3
"""List, and with --apply remove, the subroutines and functions of the solver that nothing references.

A routine counts as unreferenced when its name appears nowhere in the sources (all of Src, the MIO library
and any extra files given with --also) except in its own definition and end statement; comments are ignored.
Removal takes the routine with the comment block directly above it and its name out of "public ::" lists.

    python3 remove_unreferenced_routines.py [--also FILE ...] [--keep NAME ...] [--apply]
"""
import argparse
import glob
import os
import re

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lanczosKuboCode", "Src")
HEAD = re.compile(r"\s*(?:recursive\s+|pure\s+|elemental\s+)*(?:(?:integer|real\s*\(\s*dp\s*\)|logical|complex\s*\(\s*dp\s*\)|double\s+precision)\s+)?"
                  r"(subroutine|function)\s+(\w+)", re.I)


def code_only(text):
    return "\n".join(l if l.lstrip().lower().startswith("!$omp") else l.split("!")[0] for l in text.split("\n")).lower()


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--also", nargs="*", default=[])
    p.add_argument("--keep", nargs="*", default=[])
    p.add_argument("--apply", action="store_true")
    a = p.parse_args()
    files = sorted(glob.glob(os.path.join(SRC, "*.[Ff]90")) + glob.glob(os.path.join(SRC, "*.inc")))
    extra = glob.glob(os.path.join(SRC, "MIO", "*.[Ff]90")) + glob.glob(os.path.join(SRC, "math", "*.[Ff]90")) + list(a.also)
    text = {f: open(f, errors="replace").read() for f in files + extra}
    everything = "\n".join(code_only(t) for t in text.values())
    keep = {k.lower() for k in a.keep}
    total = 0
    for f in files:
        if os.path.basename(f) in ("test_sparse_diag.F90", "grabnes.f90"):
            continue
        L = text[f].split("\n")
        heads = [(i, m.group(2)) for i, l in enumerate(L) if not l.lstrip().startswith("!") for m in [HEAD.match(l)] if m]
        dead = []
        for i, name in heads:
            n = name.lower()
            uses = len(re.findall(r"(?<!\w)" + re.escape(n) + r"(?!\w)", everything))
            ends = len(re.findall(r"end[ \t]+(?:subroutine|function)[ \t]+" + re.escape(n) + r"(?!\w)", everything))
            if uses - ends - 1 <= 0 and n not in keep:
                e = next((j for j in range(i + 1, len(L)) if re.match(rf"\s*end\s+(subroutine|function)\s+{name}\b", L[j], re.I)), None)
                if e is None:
                    continue
                s = i
                while s > 0 and L[s - 1].lstrip().startswith("!") and not L[s - 1].lstrip().lower().startswith("!$omp"):
                    s -= 1
                dead.append((s, e, name))
        if not dead:
            continue
        n_lines = sum(e - s + 1 for s, e, _ in dead)
        total += n_lines
        print(f"{os.path.basename(f)}: " + ", ".join(f"{name} ({e - s + 1})" for s, e, name in dead) + f"  -> {n_lines} lines")
        if a.apply:
            for s, e, name in sorted(dead, reverse=True):
                del L[s:e + 1]
                while s < len(L) and s > 0 and not L[s].strip() and not L[s - 1].strip():
                    del L[s]
            out = []
            for l in L:
                m = re.match(r"(\s*public\s*::\s*)(.*)$", l, re.I)
                if m:
                    names = [x.strip() for x in m.group(2).split(",")]
                    left = [x for x in names if x.lower() not in {d[2].lower() for d in dead}]
                    if not left:
                        continue
                    l = m.group(1) + ", ".join(left) if len(left) != len(names) else l
                out.append(l)
            open(f, "w").write("\n".join(out))
    print(f"total: {total} lines" + (" removed" if a.apply else " in unreferenced routines (use --apply to remove them)"))


if __name__ == "__main__":
    main()
