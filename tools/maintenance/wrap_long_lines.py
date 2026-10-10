#!/usr/bin/env python3
"""Break the statements of the solver sources that are longer than 132 columns.

    python3 tools/maintenance/wrap_long_lines.py            # report only
    python3 tools/maintenance/wrap_long_lines.py --apply    # edit the sources

132 columns is the limit of free-form Fortran; longer lines compile only with a compiler option
(-ffree-line-length-none) or by tolerance of the compiler. A long line is broken after a comma or before a
binary operator that lies outside character strings, with `&` at the end and the continuation indented by
six columns more than the statement.

Left alone, and listed: preprocessor lines, OpenMP directives, comment lines, lines that are long only
because of a trailing comment (the comment is moved to its own line above instead), and lines that offer no
break point outside a string.

Self-check: for every file, the sequence of statements (continuations joined, blanks outside strings and
comments removed) must be the same before and after; otherwise nothing is written for that file.
"""
import argparse
import os
import re
import sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lanczosKuboCode", "Src"))
LIMIT = 132
TARGET = 120


def split_comment(line):
    """(code, comment) of a line; the comment includes the exclamation mark."""
    quote = ""
    for i, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = ""
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            return line[:i], line[i:]
    return line, ""


def break_points(code):
    """(column, parenthesis depth) of the places where the code may be cut, outside strings; the cut goes
    before the column."""
    pts, quote, depth = [], "", 0
    i, n = 0, len(code)
    while i < n:
        ch = code[i]
        if quote:
            if ch == quote:
                quote = ""
        elif ch in "'\"":
            quote = ch
        elif ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        elif ch == ",":
            pts.append((i + 1, depth))
        elif ch == "." and code[i:i + 5].lower() == ".and.":
            pts.append((i, depth))
        elif ch == "." and code[i:i + 4].lower() == ".or.":
            pts.append((i, depth))
        elif ch in "+-" and i > 0 and code[i - 1] == " " and i + 1 < n and code[i + 1] == " ":
            pts.append((i, depth))
        elif ch == "/" and code[i:i + 2] == "//":
            pts.append((i, depth))
            i += 1
        i += 1
    return pts


def wrap(code, indent):
    """Lines of a broken statement, or None if it cannot be broken. The cut is made at the shallowest
    parenthesis level available in the last half of the line, and there as far right as possible."""
    out, rest = [], code.rstrip()
    cont = indent + "      "
    while len(rest) > LIMIT - 2:
        lead = len(rest) - len(rest.lstrip())
        pts = [(p, d) for p, d in break_points(rest) if lead + 8 < p <= TARGET]
        if not pts:
            return None
        late = [(p, d) for p, d in pts if p >= TARGET // 2] or pts
        shallow = min(d for _, d in late)
        p = max(p for p, d in late if d == shallow)
        out.append(rest[:p].rstrip() + " &")
        rest = cont + rest[p:].lstrip()
    out.append(rest)
    return out


def canonical(lines):
    """Statements with continuations joined and blanks outside strings removed; comments dropped."""
    stmts, cur = [], ""
    for raw in lines:
        st = raw.strip()
        if st.startswith("#") or st.lower().startswith("!$omp"):
            stmts.append(re.sub(r"\s+", " ", st))
            continue
        code, _ = split_comment(raw.rstrip("\n"))
        code = code.strip()
        if not code:
            continue
        if code.startswith("&"):
            code = code[1:]
        if code.endswith("&"):
            cur += code[:-1]
            continue
        cur += code
        squeezed, quote = [], ""
        for ch in cur:
            if quote:
                squeezed.append(ch)
                if ch == quote:
                    quote = ""
            elif ch in "'\"":
                quote = ch
                squeezed.append(ch)
            elif not ch.isspace():
                squeezed.append(ch)
        stmts.append("".join(squeezed))
        cur = ""
    return stmts


def process(lines, report, name):
    out, changed = [], 0
    for n, raw in enumerate(lines, start=1):
        line = raw.rstrip("\n")
        if len(line) <= LIMIT:
            out.append(raw)
            continue
        st = line.lstrip()
        if st.startswith("#") or st.lower().startswith("!$omp"):
            report.append(f"{name}:{n}: left (preprocessor or OpenMP line, {len(line)} columns)")
            out.append(raw)
            continue
        if st.startswith("!"):
            report.append(f"{name}:{n}: left (comment line, {len(line)} columns)")
            out.append(raw)
            continue
        code, comment = split_comment(line)
        indent = re.match(r"\s*", code).group(0)
        new = []
        if comment and code.strip():
            new.append(indent + comment.rstrip())      # the comment goes on its own line above
            comment = ""
        if len(code.rstrip()) <= LIMIT:
            new.append(code.rstrip())
        else:
            if code.rstrip().endswith("&") or code.lstrip().startswith("&"):
                body = code.rstrip()
                tail = ""
                if body.endswith("&"):
                    body, tail = body[:-1].rstrip(), " &"
                w = wrap(body, indent)
                if w is not None:
                    w[-1] += tail
            else:
                w = wrap(code, indent)
            if w is None or any(len(x) > LIMIT for x in w):
                report.append(f"{name}:{n}: left (no break point outside a string, {len(line)} columns)")
                out.append(raw)
                continue
            new += w
        out += [x + "\n" for x in new]
        changed += 1
    return out, changed


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--apply", action="store_true")
    a = p.parse_args()
    files = []
    for d, _, fs in os.walk(ROOT):
        files += [os.path.join(d, f) for f in fs if f.endswith(".F90") or f.endswith(".f90")]
    total, report, bad = 0, [], 0
    for path in sorted(files):
        lines = open(path, errors="replace").readlines()
        name = os.path.relpath(path, ROOT)
        new, changed = process(lines, report, name)
        if not changed:
            continue
        if canonical(lines) != canonical(new):
            print(f"{name}: SELF-CHECK FAILED, file not changed")
            bad += 1
            continue
        still = sum(len(x.rstrip("\n")) > LIMIT for x in new)
        print(f"{name}: {changed} lines broken, {still} left longer than {LIMIT} columns")
        total += changed
        if a.apply:
            open(path, "w").writelines(new)
    for r in report:
        print(r)
    print(f"\n{total} lines broken" + ("" if a.apply else " (report only; --apply edits)"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
