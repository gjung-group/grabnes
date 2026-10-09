#!/usr/bin/env python3
"""Insert "not set" markers (NaN) for real and complex local variables of a routine that may be used unset.

Two kinds of variables are marked:

  locals     the names given with --names (for instance the list printed by gfortran -Wmaybe-uninitialized),
             at the first executable statement of the routine;
  privates   every real or complex variable in the PRIVATE clause of an OpenMP PARALLEL DO of the routine,
             at the top of each iteration (a PRIVATE variable is undefined at the start of every iteration,
             whatever was assigned before the loop).

A model that then uses such a variable without setting it produces NaN, which the solver stops on
(HamCheckFinite), identically for every compiler. Integers, logicals and module variables are left alone.

    python3 mark_unset_variables.py FILE ROUTINE [--names a,b,c] [--no-privates] [--dry-run]

The routine must already contain "use, intrinsic :: ieee_arithmetic" with ieee_value and ieee_quiet_nan, or
the script adds it. Markers are written between "! >>> unset markers" and "! <<< unset markers" lines, so a
second run replaces them.
"""
import argparse
import re
import sys

DECL = re.compile(r"^\s*(real\s*\(\s*dp\s*\)|real\s*\*\s*8|double\s+precision|complex\s*\(\s*dp\s*\)|complex\s*\*\s*16)"
                  r"\s*(?:,[^:]*)?::\s*(.*)$", re.I)
BEGIN, END = "! >>> unset markers", "! <<< unset markers"


def split_names(text):
    """Names of a declaration list, without dimensions or initialisers."""
    depth, cur, out = 0, "", []
    for ch in text.split("!")[0]:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    out.append(cur)
    names = []
    for item in out:
        m = re.match(r"\s*([A-Za-z_]\w*)", item)
        if m and "=" not in item.split("(")[0]:
            names.append(m.group(1))
        elif m:
            names.append(m.group(1))
    return names


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("file")
    p.add_argument("routine")
    p.add_argument("--names", default="")
    p.add_argument("--no-privates", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    L = open(a.file, errors="replace").read().split("\n")
    # drop markers of an earlier run
    out, skip = [], False
    for l in L:
        if l.strip() == BEGIN:
            skip = True
        elif l.strip() == END:
            skip = False
        elif not skip:
            out.append(l)
    L = out

    start = next(i for i, l in enumerate(L) if re.match(rf"\s*subroutine\s+{a.routine}\b", l, re.I))
    end = next(i for i in range(start + 1, len(L)) if re.match(rf"\s*end\s+subroutine\s+{a.routine}\b", l := L[i], re.I))

    # declarations of the routine (continuation lines joined)
    kind, i, last_decl = {}, start + 1, start
    while i < end:
        line = L[i]
        j = i
        while line.split("!")[0].rstrip().endswith("&") and j + 1 < end:
            j += 1
            line = line.split("!")[0].rstrip()[:-1] + " " + L[j].lstrip().lstrip("&")
        m = DECL.match(line)
        if m:
            t = "complex" if m.group(1).lower().startswith("complex") else "real"
            for n in split_names(m.group(2)):
                kind[n.lower()] = t
        s = line.strip().lower()
        if m or re.match(r"(use\b|implicit\b|integer\b|logical\b|character\b|type\s*\(|real\b|complex\b|double\b|"
                         r"external\b|parameter\b|save\b|common\b|dimension\b|intrinsic\b)", s):
            last_decl = j
        elif s and not s.startswith(("!", "#", "&")) and last_decl > start and i > last_decl:
            break
        i = j + 1
    first_exec = last_decl + 1

    def marker(n, indent):
        return f"{indent}{n} = hamUnset" if kind[n] == "real" else f"{indent}{n} = cmplx(hamUnset,hamUnset,dp)"

    # OpenMP PARALLEL DO regions of the routine: (line of the "do", names)
    regions = []
    if not a.no_privates:
        i = first_exec
        while i < end:
            if re.match(r"\s*!\$OMP\s+PARALLEL\s+DO", L[i], re.I):
                d, j = L[i], i
                while L[j].rstrip().endswith("&"):
                    j += 1
                    d += " " + L[j]
                names = [x.strip().lower() for m in re.findall(r"(?<!FIRST)PRIVATE\s*\(([^)]*)\)", d, re.I)
                         for x in m.split(",")]
                k = j + 1
                while not re.match(r"\s*do\s+\w+\s*=", L[k], re.I):
                    k += 1
                regions.append((k, [n for n in dict.fromkeys(names) if n in kind]))
                i = k
            i += 1

    wanted = [n.strip().lower() for n in a.names.split(",") if n.strip()]
    locals_ = [n for n in dict.fromkeys(wanted) if n in kind]
    skipped = [n for n in wanted if n not in kind]

    print(f"{a.routine}: lines {start + 1}-{end + 1}; first executable statement at line {first_exec + 1}: {L[first_exec].strip()[:60]}")
    print(f"  locals marked: {len(locals_)}; not real/complex locals of the routine, skipped: {len(skipped)} {skipped[:12]}")
    for k, names in regions:
        print(f"  parallel loop at line {k + 1}: {len(names)} private real/complex variables marked")
    if a.dry_run:
        return

    # insert from the bottom up so that line numbers stay valid
    for k, names in sorted(regions, reverse=True):
        if not names:
            continue
        indent = re.match(r"\s*", L[k + 1] if L[k + 1].strip() else L[k]).group(0) or "   "
        L[k + 1:k + 1] = [f"{indent}{BEGIN}"] + [marker(n, indent) for n in names] + [f"{indent}{END}"]
    block = [f"   {BEGIN}", "   ! Value of a variable that has not been set: any use of it gives NaN, and the run stops",
             "   ! in HamCheckFinite instead of continuing with whatever the memory holds.",
             "   hamUnset = ieee_value(hamUnset, ieee_quiet_nan)"] + [marker(n, "   ") for n in locals_] + [f"   {END}"]
    L[first_exec:first_exec] = block
    # declaration of hamUnset and the intrinsic module
    body = "\n".join(L[start:first_exec])
    add = []
    if not re.search(r"\bhamUnset\b", body.replace("hamUnset = ", "")):
        add.append("   real(dp) :: hamUnset                    ! NaN: marks variables that have not been set")
    if add:
        L[first_exec:first_exec] = [f"   {BEGIN}"] + add + [f"   {END}"]
    if "ieee_arithmetic" not in body:
        u = next(i for i in range(start + 1, first_exec) if L[i].strip().lower().startswith("use "))
        L[u:u] = ["   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan"]
    open(a.file, "w").write("\n".join(L))


if __name__ == "__main__":
    main()
