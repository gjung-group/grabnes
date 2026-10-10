#!/usr/bin/env python3
"""Add the execution trace of the original code to the routines that lack it.

    python3 tools/maintenance/add_debug_trace.py            # report only
    python3 tools/maintenance/add_debug_trace.py --apply    # edit the sources

The convention (docs/development/coding-conventions.md): the first executable statement of a subroutine is

    #ifdef DEBUG
       call MIO_Debug('Name',0)
    #endif /* DEBUG */

and the same call with the argument 1 stands before every `return` and at the end. A build with -DDEBUG then
writes the sequence of routines to debug.log; other builds are not affected.

The trace is added to the module-level subroutines of lanczosKuboCode/Src/*.F90 that
  - have no MIO_Debug call yet,
  - are not pure, elemental or recursive,
  - are never called from inside a parallel region, directly or through a routine that is (several
    threads would write at once),
  - are not called directly from inside a loop or from a recursive routine (one line per atom, bond or
    k-point), unless named with
    --trace-looped: the loops over magnetic fields, shifts and k-points call routines that are worth tracing.
Routines whose layout the tool cannot handle safely (preprocessor blocks that mix declarations and
statements, a `return` in a one-line `if` with a continuation, ...) are left alone and listed.

A one-line `if (condition) return` becomes a block so that the trace call can precede the return; nothing
else of the code is changed. Functions are not traced.
"""
import argparse
import os
import re
import sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lanczosKuboCode", "Src"))

SUB = re.compile(r"^\s*((?:(?:recursive|pure|elemental|impure)\s+)*)subroutine\s+(\w+)", re.I)
FUN = re.compile(r"^\s*(?:(?:recursive|pure|elemental|impure)\s+)*(?:(?:integer|real|complex|logical|character|double\s+precision|type)\b[^!]*?\s+)?function\s+(\w+)\s*\(", re.I)
END_SUB = re.compile(r"^\s*end\s*subroutine\b", re.I)
END_FUN = re.compile(r"^\s*end\s*function\b", re.I)
CALL = re.compile(r"\bcall\s+(\w+)", re.I)
DO = re.compile(r"^\s*(?:\w+\s*:\s*)?do\b(?!\s*\w+\s*=\s*[^,]*$)", re.I)
END_DO = re.compile(r"^\s*end\s*do\b", re.I)
OMP_BEGIN = re.compile(r"^\s*!\$omp\s+parallel\b", re.I)
OMP_END = re.compile(r"^\s*!\$omp\s+end\s+parallel\b", re.I)
DECL = re.compile(
    r"^\s*(use\b|implicit\b|integer\b|real\b|complex\b|logical\b|character\b|double\s+precision\b|type\s*\(|class\s*\(|"
    r"parameter\b|dimension\b|external\b|intrinsic\b|save\b|common\b|data\b|include\b|procedure\b|allocatable\b|"
    r"pointer\b|target\b|intent\b|optional\b|namelist\b|equivalence\b|interface\b|end\s*interface\b|import\b|"
    r"type\b\s*(,|::|\w+\s*$)|end\s*type\b)", re.I)
RETURN_ALONE = re.compile(r"^(\s*)return\s*(!.*)?$", re.I)
RETURN_IF = re.compile(r"^(\s*)if\s*\((.*)\)\s*return\s*(!.*)?$", re.I)
RETURN_ANY = re.compile(r"\breturn\b", re.I)


def code_of(line):
    """The statement part of a line: without comment (quotes respected) and trailing blanks."""
    out, quote = [], ""
    for ch in line:
        if quote:
            if ch == quote:
                quote = ""
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            break
        out.append(ch)
    return "".join(out).rstrip()


class Routine:
    def __init__(self, name, kind, prefix, start):
        self.name, self.kind, self.prefix, self.start = name, kind, prefix.lower(), start
        self.end = None
        self.calls = []          # (callee, in_loop)
        self.internal = False    # contained in another routine


def parse(lines):
    """Module-level and internal routines of one file, with their calls."""
    routines, stack = [], []
    loop = omp = 0
    for i, raw in enumerate(lines):
        if OMP_BEGIN.match(raw):
            omp += 1
            continue
        if OMP_END.match(raw):
            omp = max(0, omp - 1)
            continue
        s = code_of(raw)
        if not s.strip() or s.lstrip().startswith("#"):
            continue
        m = SUB.match(s)
        f = None if m else FUN.match(s)
        if m or f:
            r = Routine(m.group(2) if m else f.group(1), "subroutine" if m else "function", m.group(1) if m else "", i)
            r.internal = bool(stack)
            r.saved = (loop, omp)
            loop = omp = 0
            stack.append(r)
            routines.append(r)
            continue
        if END_SUB.match(s) or END_FUN.match(s):
            if stack:
                r = stack.pop()
                r.end = i
                loop, omp = r.saved
            continue
        if not stack:
            continue
        if END_DO.match(s):
            loop = max(0, loop - 1)
        elif re.match(r"^\s*(?:\w+\s*:\s*)?do\b", s, re.I) and not re.match(r"^\s*do\w", s, re.I):
            loop += 1
        for c in CALL.findall(s):
            stack[-1].calls.append((c.lower(), omp > 0, loop > 0))
    return routines


def find_hot(all_routines):
    """Two sets of routine names. `parallel`: called from inside a parallel region, directly or through
    such a routine, or from a function (a function may be evaluated anywhere). `looped`: called directly
    from inside a loop (not propagated: the outer loops over fields or shifts run a whole calculation)."""
    parallel, looped = set(), set()
    for r in all_routines:
        for callee, in_omp, in_loop in r.calls:
            if in_omp or r.kind == "function":
                parallel.add(callee)
            if in_loop or 'recursive' in r.prefix:
                looped.add(callee)
    changed = True
    while changed:
        changed = False
        for r in all_routines:
            if r.name.lower() in parallel:
                for callee, _, _ in r.calls:
                    if callee not in parallel:
                        parallel.add(callee)
                        changed = True
    return parallel, looped


def block(name, flag, indent):
    return ["#ifdef DEBUG\n", f"{indent}call MIO_Debug('{name}',{flag})\n", "#endif /* DEBUG */\n"]


def instrument(lines, r):
    """New lines of routine r with the trace, or a string saying why it is left alone."""
    body = lines[r.start:r.end + 1]
    # statements joined over continuations: (first line index, last line index, text)
    stmts, i = [], 0
    cpp_depth = []
    depth = 0
    contains = None
    while i < len(body):
        raw = body[i]
        st = raw.lstrip()
        if st.startswith("#"):
            if re.match(r"#\s*if", st):
                depth += 1
            elif re.match(r"#\s*endif", st):
                depth -= 1
            stmts.append((i, i, "#", depth))
            i += 1
            continue
        s = code_of(raw)
        j = i
        text = s
        while text.rstrip().endswith("&") and j + 1 < len(body):
            j += 1
            nxt = code_of(body[j])
            while not nxt.strip() and j + 1 < len(body):      # comment lines inside a continuation
                j += 1
                nxt = code_of(body[j])
            text = text.rstrip()[:-1] + " " + nxt.lstrip().lstrip("&")
        stmts.append((i, j, text, depth))
        i = j + 1
    if depth != 0:
        return "unbalanced preprocessor conditionals"

    # end of the routine's own statements: an internal `contains`, or the end statement
    last = len(stmts) - 1
    for k, (a, b, text, d) in enumerate(stmts):
        if k > 0 and re.match(r"^\s*contains\s*$", text, re.I):
            contains = k
            break
    stop = contains if contains is not None else last
    if stmts[stop][3] != 0:
        return "the end of the routine lies inside a preprocessor conditional"

    # first executable statement
    first = None
    k = 1
    # the subroutine statement itself may be continued: it is stmts[0]
    while k < stop:
        a, b, text, d = stmts[k]
        if text == "#" or not text.strip() or DECL.match(text):
            k += 1
            continue
        first = k
        break
    if first is None:
        first = stop
    if stmts[first][3] != 0:
        # inside a conditional: go back to its opening line, provided nothing but declarations precede it there
        k = first
        d0 = stmts[first][3]
        opening = None
        while k > 0:
            k -= 1
            if stmts[k][2] == "#" and stmts[k][3] == d0 and re.match(r"#\s*if", body[stmts[k][0]].lstrip()):
                if d0 == 1:
                    opening = k
                break
        if opening is None:
            return "the first statement lies in a nested preprocessor conditional"
        for q in range(opening + 1, first):
            if stmts[q][2] != "#" and stmts[q][2].strip():
                return "a preprocessor conditional mixes declarations and statements"
        first = opening
    # declarations after the first executable statement would be a layout the heuristic misread
    for q in range(first + 1, stop):
        t = stmts[q][2]
        if t != "#" and re.match(r"^\s*(use\b|implicit\b)", t, re.I):
            return "a use or implicit statement follows the first statement found"

    indent_m = re.match(r"\s*", body[stmts[min(first, stop)][0]] if not body[stmts[first][0]].lstrip().startswith("#") else body[r_first_code(stmts, body, first, stop)])
    indent = indent_m.group(0) if indent_m.group(0) else "   "

    new, edits = [], {}
    # exits
    for k in range(first, stop):
        a, b, text, d = stmts[k]
        if text == "#" or not RETURN_ANY.search(text):
            continue
        if re.search(r"['\"][^'\"]*\breturn\b[^'\"]*['\"]", text) and not RETURN_ALONE.match(text) and not RETURN_IF.match(text):
            continue                                     # the word inside a string
        if d != 0 and False:
            return "a return inside a preprocessor conditional"
        m = RETURN_ALONE.match(body[a].rstrip("\n")) if a == b else None
        if m:
            edits[a] = block(r.name, 1, m.group(1)) + [body[a]]
            continue
        m = RETURN_IF.match(body[a].rstrip("\n")) if a == b else None
        if m:
            ind, cond, com = m.group(1), m.group(2), m.group(3) or ""
            if cond.count("(") != cond.count(")"):
                return "a return in an if statement the tool cannot split"
            edits[a] = ([f"{ind}if ({cond}) then" + (f"   {com}" if com else "") + "\n"] + block(r.name, 1, ind + "   ")
                        + [f"{ind}   return\n", f"{ind}end if\n"])
            continue
        return "a return statement in a form the tool does not handle: " + text.strip()[:60]

    a_first = stmts[first][0]
    a_stop = stmts[stop][0]
    for i, raw in enumerate(body):
        if i == a_first and i != a_stop:
            new += block(r.name, 0, indent) + ["\n"]
        if i == a_stop:
            if a_first == a_stop:
                new += block(r.name, 0, indent)
            # blank line before the closing block unless one is there
            if new and new[-1].strip():
                new.append("\n")
            new += block(r.name, 1, indent) + ["\n"]
        new += edits.get(i, [raw])
    return new


def r_first_code(stmts, body, first, stop):
    for q in range(first, stop + 1):
        if stmts[q][2] != "#" and stmts[q][2].strip():
            return stmts[q][0]
    return stmts[stop][0]


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--apply", action="store_true")
    p.add_argument("--skip", nargs="*", default=[], help="routine names to leave alone")
    p.add_argument("--trace-looped", nargs="*", default=[],
                   help="routines called inside a loop that are traced nevertheless (few calls per run)")
    p.add_argument("--root", default=ROOT, help="directory of the sources (default: the solver in this repository)")
    p.add_argument("--only-file", default=None, help="edit this file only (the others are read for the calls)")
    a = p.parse_args()
    globals()["ROOT"] = a.root
    keep = {s.lower() for s in a.trace_looped}
    files = sorted(f for f in os.listdir(ROOT) if f.endswith(".F90"))
    # the include files of a private build are parsed for calls as well, if present
    extra = sorted(f for f in os.listdir(ROOT) if f.endswith(".inc") or f.endswith(".f90"))
    parsed, everything = {}, []
    for f in files + extra:
        lines = open(os.path.join(ROOT, f), errors="replace").readlines()
        rs = parse(lines)
        parsed[f] = (lines, rs)
        everything += rs
    hot, looped = find_hot(everything)
    skip = {s.lower() for s in a.skip}
    total = {"added": 0, "present": 0, "parallel": 0, "loop": 0, "kind": 0, "layout": 0}
    for f in files:
        if a.only_file and f != a.only_file:
            continue
        lines, rs = parsed[f]
        out, changed = list(lines), False
        for r in sorted((r for r in rs if r.kind == "subroutine" and not r.internal and r.end is not None),
                        key=lambda r: -r.start):
            text = "".join(lines[r.start:r.end + 1])
            if re.search(r"MIO_Debug\s*\(", text, re.I):
                total["present"] += 1
                continue
            if r.prefix.strip() or r.name.lower() in skip:
                total["kind"] += 1
                print(f"{f}: {r.name}: not traced (pure, elemental, recursive, or excluded)")
                continue
            if r.name.lower() in hot:
                total["parallel"] += 1
                print(f"{f}: {r.name}: not traced (called inside a parallel region)")
                continue
            if r.name.lower() in looped and r.name.lower() not in keep:
                total["loop"] += 1
                print(f"{f}: {r.name}: not traced (called inside a loop; --trace-looped NAME to trace it)")
                continue
            res = instrument(lines, r)
            if isinstance(res, str):
                total["layout"] += 1
                print(f"{f}: {r.name}: LEFT ALONE: {res}")
                continue
            out[r.start:r.end + 1] = res
            changed = True
            total["added"] += 1
            print(f"{f}: {r.name}: trace added")
        if changed and a.apply:
            open(os.path.join(ROOT, f), "w").writelines(out)
    print("\nsummary: " + ", ".join(f"{k} {v}" for k, v in total.items()) + ("" if a.apply else "   (report only; --apply edits)"))


if __name__ == "__main__":
    sys.exit(main())
