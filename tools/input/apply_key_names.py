#!/usr/bin/env python3
"""Rename input keys in the solver sources according to a table, keeping the former names valid.

    python3 tools/input/apply_key_names.py            # report only
    python3 tools/input/apply_key_names.py --apply    # edit the sources and the alias table

The table is docs/development/input-key-renaming-proposal.md: every row `| old | new | ... |` with a
non-empty second column is applied; rows without a new name are left as they are. For each pair

  - the string passed to MIO_InputParameter / MIO_InputSearchLabel in lanczosKuboCode/Src is changed
    from the old name to the new one (the Fortran variables are not renamed);
  - a line  call InputAddAlias('new','old')  is written to Src/MIO/input_aliases.inc, so that an input
    file with the old name gives the same run.

Lines of input_aliases.inc that this tool did not write (above the marker line) are kept. The tool refuses
to do anything if two keys would become equal in the form the input library compares (case, `_` and `-`
ignored), or if a new name equals a key that already exists. It can be run again after the table has
been edited: a pair whose old name is no longer in the sources is skipped if its alias exists already.
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "lanczosKuboCode", "Src")
TABLE = os.path.join(ROOT, "docs", "development", "input-key-renaming-proposal.md")
INC = os.path.join(SRC, "MIO", "input_aliases.inc")
MARK = "! ---- pairs written by tools/input/apply_key_names.py (do not edit below this line by hand) ----"
ROW = re.compile(r"^\|\s*`([^`]+)`\s*\|\s*(?:`([^`]*)`)?\s*\|")
READ = re.compile(r"((?:MIO_InputParameter|MIO_InputSearchLabel|InputSearchLabel|InputParameter)\s*\(\s*)(['\"])([^'\"]+)\2", re.I)
ALIAS = re.compile(r"call\s+InputAddAlias\(\s*'([^']+)'\s*,\s*'([^']+)'\s*\)", re.I)


def form(key):
    return key.lower().replace("_", "").replace("-", "")


def source_files():
    out = []
    for d, _, fs in os.walk(SRC):
        out += [os.path.join(d, f) for f in fs if f.endswith((".F90", ".f90", ".inc")) and f != "input_aliases.inc"]
    return sorted(out)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--apply", action="store_true")
    a = p.parse_args()

    pairs = []
    for line in open(TABLE):
        m = ROW.match(line)
        if m and m.group(2) and m.group(1) not in ("Present key",):
            pairs.append((m.group(1), m.group(2).strip()))

    texts = {f: open(f, errors="replace").read() for f in source_files()}
    present = {}
    for f, t in texts.items():
        for m in READ.finditer(t):
            if not m.group(3).startswith("&"):
                present.setdefault(form(m.group(3)), m.group(3))

    inc = open(INC).read() if os.path.exists(INC) else ""
    head = inc.split(MARK)[0].rstrip("\n") + "\n"
    kept = ALIAS.findall(head)                      # (new, old) pairs maintained by hand
    written = ALIAS.findall(inc.split(MARK)[1]) if MARK in inc else []

    problems, todo, skipped = [], [], 0
    new_forms = {}
    for old, new in pairs:
        if "." not in new:
            problems.append(f"{old}: the new name '{new}' has no section")
        if form(new) in new_forms and new_forms[form(new)] != old:
            problems.append(f"'{new}' is proposed for both '{new_forms[form(new)]}' and '{old}'")
        new_forms[form(new)] = old
        if form(new) in present and form(new) != form(old) and (new, old) not in written:
            problems.append(f"{old}: the new name '{new}' is already a key of the solver")
        if form(old) not in present:
            if (new, old) in written:
                skipped += 1                         # applied in an earlier run
                todo.append((old, new))
                continue
            problems.append(f"{old}: not read anywhere in the sources")
            continue
        todo.append((old, new))
    if problems:
        print("\n".join(problems))
        print(f"\n{len(problems)} problem(s); nothing changed")
        return 1

    by_form = {form(o): n for o, n in todo}
    changed = 0
    for f, t in texts.items():
        def sub(m):
            nonlocal changed
            new = by_form.get(form(m.group(3)))
            if new is None or m.group(3) == new:
                return m.group(0)
            changed += 1
            return m.group(1) + m.group(2) + new + m.group(2)
        t2 = READ.sub(sub, t)
        if t2 != t and a.apply:
            open(f, "w").write(t2)

    # a key that had a former name already (hand-maintained pair) keeps it under its new name
    extra = [(by_form[form(n)], o) for n, o in kept if form(n) in by_form]
    lines = [head.rstrip("\n"), "", MARK]
    lines += [f"call InputAddAlias('{n}','{o}')" for o, n in sorted(todo, key=lambda x: x[1].lower())]
    lines += [f"call InputAddAlias('{n}','{o}')" for n, o in extra]
    if a.apply:
        open(INC, "w").write("\n".join(lines) + "\n")
    print(f"{len(todo)} pairs ({skipped} applied earlier), {changed} strings changed in the sources, "
          f"{len(extra)} hand-maintained former name(s) carried over" + ("" if a.apply else "   (report only; --apply edits)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
