#!/usr/bin/env python3
"""Move a top-level block of a long routine into an internal procedure of that routine.

    python3 tools/maintenance/outline_block.py FILE ROUTINE LINE NAME "what the block does" [--apply]

LINE is the first line of a block (`if (...) then`, `do ...`, `select case`) of ROUTINE that is not inside a
loop (it may be inside `if` blocks). The block becomes

    contains
       subroutine NAME()        ! internal procedure: it sees the variables of ROUTINE
          ... the block, unchanged ...
       end subroutine NAME

and its place in ROUTINE is taken by `call NAME()`. For an `if (condition) then ... end if` without an
`else` at its own level, the condition stays in ROUTINE: `if (condition) call NAME()`.

Why an internal procedure: it shares every variable of the routine (host association), so no variable has
to be passed, renamed or re-declared and the statements are moved without any change. This is the step that
can be verified exactly; turning an internal procedure into a module procedure with an argument list is a
separate step.

The tool refuses a block that
  - lies inside a loop of the routine, or inside an OpenMP parallel region that starts outside it
    (inside such a region the variables of the routine that are private to the region would be shared in the
    procedure),
  - contains `return`, `goto`, a statement label, or `cycle` / `exit` that leave the block,
  - has unbalanced preprocessor conditionals.
OpenMP directives directly before and after a `do` block (`!$OMP PARALLEL DO` ... `!$OMP END PARALLEL DO`)
are moved with it.
"""
import argparse
import re
import sys


def code(line):
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
    return "".join(out).rstrip()


OPEN = re.compile(r"^(\w+\s*:\s*)?(if\s*\(.*\)\s*then$|do(\s+\w+\s*=|\s+while|\s*$)|select\s+case)", re.I)
CLOSE = re.compile(r"^end\s*(if|do|select)\b", re.I)
OMP = re.compile(r"^\s*!\$omp\b", re.I)


def statements(lines, a, b):
    """(first, last, text, is_omp) of the statements of lines a..b-1, continuations joined."""
    out, i = [], a
    while i < b:
        raw = lines[i]
        if OMP.match(raw):
            j = i
            while lines[j].rstrip().endswith("&") and j + 1 < b:
                j += 1
            out.append((i, j, " ".join(l.strip() for l in lines[i:j + 1]), True))
            i = j + 1
            continue
        c = code(raw)
        if c.strip().startswith("#") or not c.strip():
            out.append((i, i, c.strip(), False))
            i += 1
            continue
        j = i
        while c.rstrip().endswith("&") and j + 1 < b:
            j += 1
            c = c.rstrip()[:-1] + " " + code(lines[j]).strip().lstrip("&")
        out.append((i, j, c.strip(), False))
        i = j + 1
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("file")
    p.add_argument("routine")
    p.add_argument("line", type=int)
    p.add_argument("name")
    p.add_argument("doc")
    p.add_argument("--apply", action="store_true")
    a = p.parse_args()

    lines = open(a.file).read().split("\n")
    try:
        r0 = next(i for i, l in enumerate(lines) if re.match(r"^\s*subroutine\s+" + a.routine + r"\b", l, re.I))
        r1 = next(i for i, l in enumerate(lines) if i > r0 and re.match(r"^\s*end\s+subroutine\s+" + a.routine + r"\b", l, re.I))
    except StopIteration:
        sys.exit(f"routine {a.routine} not found")
    if any(re.match(r"^\s*subroutine\s+" + a.name + r"\b", l, re.I) for l in lines):
        sys.exit(f"a routine named {a.name} exists already")
    contains = next((i for i in range(r0, r1) if re.match(r"^\s*contains\s*$", code(lines[i]), re.I)), None)
    body_end = contains if contains is not None else r1

    st = statements(lines, r0 + 1, body_end)
    depth, omp_depth, start_k, end_k = 0, 0, None, None
    kinds = []                 # kinds of the blocks that are open: "if", "do", "select"
    for k, (i, j, text, is_omp) in enumerate(st):
        if is_omp:
            t = text.lower()
            if start_k is None and re.match(r"!\$omp\s+parallel\b(?!\s+do)", t):
                omp_depth += 1
            elif start_k is None and re.match(r"!\$omp\s+end\s+parallel\b(?!\s+do)", t):
                omp_depth -= 1
            continue
        if text.startswith("#") or not text:
            continue
        text = re.sub(r"^\d+\s+", "", text)        # a statement label does not hide a block statement
        if i == a.line - 1:
            if "do" in kinds:
                sys.exit(f"line {a.line} lies inside a loop of {a.routine}")
            base_depth = depth
            if omp_depth != 0:
                sys.exit(f"line {a.line} lies inside an OpenMP parallel region that starts before it")
            if not OPEN.match(text):
                sys.exit(f"line {a.line} does not start a block: {text[:60]}")
            start_k = k
        if OPEN.match(text):
            depth += 1
            kinds.append("do" if re.match(r"^(\w+\s*:\s*)?do\b", text, re.I) else "if")
        elif CLOSE.match(text):
            depth -= 1
            if kinds:
                kinds.pop()
            if start_k is not None and depth == base_depth:
                end_k = k
                break
    if start_k is None or end_k is None:
        sys.exit(f"no complete block starts at line {a.line}")

    first, last = st[start_k][0], st[end_k][1]
    head = st[start_k][2]
    inner = st[start_k + 1:end_k]

    # an if block without else at its own level keeps its condition in the routine
    keep_condition = False
    m = re.match(r"^if\s*\((.*)\)\s*then$", head, re.I)
    if m:
        d, has_else = 0, False
        for i, j, text, is_omp in inner:
            if is_omp or text.startswith("#") or not text:
                continue
            text = re.sub(r"^\d+\s+", "", text)
            if OPEN.match(text):
                d += 1
            elif CLOSE.match(text):
                d -= 1
            elif d == 0 and re.match(r"^else\b", text, re.I):
                has_else = True
        keep_condition = not has_else
    if keep_condition:
        move_first, move_last = st[start_k][1] + 1, st[end_k][0] - 1
    else:
        move_first, move_last = first, last
        # OpenMP directives that belong to a do block
        if re.match(r"^(\w+\s*:\s*)?do\b", head, re.I):
            k = start_k - 1
            while k >= 0 and st[k][3]:
                move_first = st[k][0]
                k -= 1
            k = end_k + 1
            while k < len(st) and st[k][3] and re.match(r"!\$omp\s+end\b", st[k][2].lower()):
                move_last = st[k][1]
                k += 1

    moved = lines[move_first:move_last + 1]
    checked = statements(lines, move_first, move_last + 1)
    cpp = 0
    loop = 0
    for i, j, text, is_omp in checked:
        if text.startswith("#"):
            cpp += 1 if re.match(r"#\s*if", text) else (-1 if re.match(r"#\s*endif", text) else 0)
            if cpp < 0:
                sys.exit("unbalanced preprocessor conditionals in the block")
            continue
        if is_omp or not text:
            continue
        t = text.lower()
        if re.match(r"^(\w+\s*:\s*)?do\b", t):
            loop += 1
        elif re.match(r"^end\s*do\b", t):
            loop -= 1
        if re.search(r"(^|[\s)])return\s*$", t) or re.search(r"\bgo\s*to\b", t) or re.match(r"^\d+\s", t):
            sys.exit(f"line {i + 1}: return, goto or statement label in the block: {text[:60]}")
        if re.search(r"(^|[\s)])(cycle|exit)\b", t) and loop <= 0:
            sys.exit(f"line {i + 1}: cycle or exit that leaves the block: {text[:60]}")
    if cpp != 0:
        sys.exit("unbalanced preprocessor conditionals in the block")

    ind = re.match(r"\s*", lines[first]).group(0)
    cond_text = " ".join(l.strip() for l in lines[first:st[start_k][1] + 1])
    if keep_condition:
        cond = re.match(r"^(.*\))\s*then\s*(!.*)?$", cond_text, re.I)
        call = [f"{ind}{cond.group(1)} call {a.name}()"] if cond and "&" not in cond_text and len(cond_text) < 110 else \
            lines[first:st[start_k][1] + 1] + [f"{ind}   call {a.name}()", f"{ind}end if"]
        replaced_first, replaced_last = first, last
        # the body is moved as it is, indentation included: a continued character string would change otherwise
        body = moved
    else:
        call = [f"{ind}call {a.name}()"]
        replaced_first, replaced_last = move_first, move_last
        body = moved

    doc = ["!> " + a.doc + " (internal procedure of " + a.routine + ": it uses the variables of that routine)"]
    proc = doc + [f"subroutine {a.name}()", ""] + \
        ["#ifdef DEBUG", f"   call MIO_Debug('{a.name}',0)", "#endif /* DEBUG */", ""] + body + \
        ["", "#ifdef DEBUG", f"   call MIO_Debug('{a.name}',1)", "#endif /* DEBUG */", "", f"end subroutine {a.name}", ""]

    new = lines[:replaced_first] + call + lines[replaced_last + 1:]
    shift = len(call) - (replaced_last - replaced_first + 1)
    end_line = r1 + shift
    if contains is None:
        new[end_line:end_line] = ["contains", ""] + proc
    else:
        new[end_line:end_line] = proc
    print(f"{a.routine}: lines {move_first + 1}-{move_last + 1} ({move_last - move_first + 1} lines) -> internal procedure {a.name}"
          + ("; the condition stays in the routine" if keep_condition else "") + ("" if a.apply else "   (report only; --apply edits)"))
    if a.apply:
        open(a.file, "w").write("\n".join(new))


if __name__ == "__main__":
    main()
