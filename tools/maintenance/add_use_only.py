#!/usr/bin/env python3
"""Give the module imports of the solver sources an explicit `only :` list.

    python3 tools/maintenance/add_use_only.py            # report only
    python3 tools/maintenance/add_use_only.py --apply    # edit the sources

The convention (docs/development/coding-conventions.md): a routine imports what it needs with
`use module, only : names`. This tool rewrites every bare `use math`, `use constants`, `use kuboarrays`
and `use name`: the list holds the public names of the module that occur in the scope of the statement
(the routine, or the whole module for an import at module level). `use mio` stays as it is: importing the
whole library is the convention of the original code.

A statement whose scope uses none of the names is reported and left alone. The edit cannot change a result:
a name that is missing from a list is a compilation error.
"""
import argparse
import os
import re
import sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lanczosKuboCode", "Src"))
MODULES = ("math", "constants", "kuboarrays", "name")
UNIT = re.compile(r"^\s*(?:(?:recursive|pure|elemental|impure)\s+)*(?:[\w()*=, ]+?\s+)?(module|subroutine|function|program)\s+(\w+)", re.I)
END_UNIT = re.compile(r"^\s*end\s*(module|subroutine|function|program)\b", re.I)


def code_of(line):
    out, quote = [], ""
    for ch in line:
        if quote:
            if ch == quote:
                quote = ""
            out.append(" ")
            continue
        if ch in "'\"":
            quote = ch
            out.append(" ")
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def public_names(module):
    """Public names of a module of the solver, in the order of their declaration."""
    if module == "math":
        text = open(os.path.join(ROOT, "math", "math.f90")).read()
        return re.findall(r"^\s*(?:subroutine|function)\s+(\w+)", text, re.I | re.M)
    path = os.path.join(ROOT, module + ".f90")
    names = []
    for line in open(path):
        c = code_of(line)
        if "::" not in c:
            continue
        decl, rhs = c.split("::", 1)
        if module == "constants" and "public" not in decl.lower():
            continue
        if re.match(r"^\s*(use|implicit|private|public\s*$)", decl, re.I):
            continue
        depth, item, items = 0, "", []
        for ch in rhs:
            if ch in "([":
                depth += 1
            elif ch in ")]":
                depth -= 1
            if ch == "," and depth == 0:
                items.append(item)
                item = ""
            else:
                item += ch
        items.append(item)
        for it in items:
            m = re.match(r"\s*(\w+)", it)
            if m:
                names.append(m.group(1))
    return names


def process(path, names, apply):
    lines = open(path, errors="replace").read().split("\n")
    stack, scopes = [], {}                 # line index of a unit start -> line index of its end
    for i, raw in enumerate(lines):
        c = code_of(raw)
        if END_UNIT.match(c):
            if stack:
                scopes[stack.pop()] = i
            continue
        m = UNIT.match(c)
        if m and not re.match(r"^\s*module\s+procedure\b", c, re.I):
            stack.append(i)
    starts = sorted(scopes)
    changed, report = 0, []
    for i, raw in enumerate(lines):
        m = re.match(r"^(\s*)use\s+(\w+)\s*(!.*)?$", raw, re.I)
        if not m or m.group(2).lower() not in names:
            continue
        mod = m.group(2).lower()
        enclosing = [s for s in starts if s < i <= scopes[s]]
        if not enclosing:
            report.append(f"{os.path.basename(path)}:{i + 1}: use {mod}: no enclosing program unit found")
            continue
        s = max(enclosing)
        body = "\n".join(code_of(l) for l in lines[s:scopes[s] + 1]).lower()
        used = [n for n in names[mod] if re.search(r"(?<![\w%])" + re.escape(n.lower()) + r"(?!\w)", body)]
        if not used:
            report.append(f"{os.path.basename(path)}:{i + 1}: use {mod}: none of its names occurs here; left alone")
            continue
        head = f"{m.group(1)}use {m.group(2)},"
        head = head.ljust(len(m.group(1)) + 26) + "only : "
        out, cur = [], head
        for k, n in enumerate(used):
            piece = n + (", " if k < len(used) - 1 else "")
            if len(cur) + len(piece) > 118:
                out.append(cur.rstrip() + " &")
                cur = " " * len(head)
            cur += piece
        out.append(cur.rstrip())
        lines[i] = "\n".join(out)
        changed += 1
    if changed and apply:
        open(path, "w").write("\n".join(lines))
    return changed, report


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--apply", action="store_true")
    p.add_argument("--root", default=ROOT, help="directory of the sources to edit (default: the solver)")
    a = p.parse_args()
    names = {m: public_names(m) for m in MODULES}
    total, reports = 0, []
    for f in sorted(os.listdir(a.root)):
        if not f.endswith((".F90", ".f90", ".inc")) or f in ("constants.f90", "kuboarrays.f90", "name.f90"):
            continue
        n, r = process(os.path.join(a.root, f), names, a.apply)
        if n:
            print(f"{f}: {n} imports given an only list")
        total += n
        reports += r
    for r in reports:
        print(r)
    print(f"\n{total} imports changed" + ("" if a.apply else " (report only; --apply edits)"))


if __name__ == "__main__":
    sys.exit(main())
