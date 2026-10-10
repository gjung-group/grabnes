#!/usr/bin/env python3
"""Remove commented-out code from Fortran sources, keeping the comments that explain something.

A comment line is taken for disabled code when, after its leading "!" characters, it reads like a
Fortran statement (call, assignment, if (...), do i=..., print, write(...), allocate(...), end if, a
disabled !$OMP directive, ...), or continues such a line. Documentation comments (!> and !! blocks),
prose, and active !$OMP directives are kept. Blank lines left without neighbours are squeezed to one.

The script verifies its own result: the file with all comment-only and blank lines removed must be
identical before and after, i.e. only comments were deleted and the compiled code cannot change.

    python3 remove_commented_code.py FILE... [--dry-run] [--show N]
"""
import argparse
import re
import sys

STATEMENT = re.compile(
    r"""^(?:
        call\s+\w+\s*(?:\(|$|!) | (?:else\s*)?if\s*\( | else\b\s*$ | else\s+if\b | end\s*(?:if|do|select|where|subroutine|function)\b |
        do\s+\w+\s*=\s*[\w(] | do\s+while\s*\( | print\s*\*\s*, | print\s*'\( | write\s*\( | read\s*\( | open\s*\( | close\s*\( |
        allocate\s*\( | deallocate\s*\( | return\s*$ | cycle\s*$ | exit\s*$ | stop\b | goto\s+\d | continue\s*$ |
        use\s+\w+\s*(?:,|$) | implicit\s+none | (?:real|integer|logical|complex|character)\s*(?:\(|,|::|\*) |
        \$OMP\b | select\s+case\s*\( | case\s*\( | where\s*\( | format\s*\( | \d+\s+(?:continue|format)\b |
        [A-Za-z_]\w*(?:%\w+)*\s*(?:\([^!]*\))?\s*=\s*[^=\s]
    )""", re.I | re.X)
PROSE = re.compile(r"\b(?:the|this|that|is|are|was|were|we|not|for|with|from|which|because|should|must|when|then use|only|"
                   r"note|todo|see|where the|if the|if we|if a|if it)\b", re.I)


def is_code_comment(line, previous_was_code):
    s = line.strip()
    if not s.startswith("!") or s.startswith(("!>", "!$OMP", "!$omp")):
        return False
    if re.match(r"!!(?!\$)", s) and not re.match(r"!!\s*(call|if\s*\(|do\s+\w+\s*=|print|write\s*\()", s, re.I):
        return False                                   # "!!" starts a documentation line unless it is plainly code
    body = s.lstrip("!").strip()
    if not body:
        return False
    if previous_was_code and (body.startswith("&") or previous_was_code == "&"):
        return True
    if not STATEMENT.match(body):
        return False
    # an assignment-looking sentence ("H = -hopp is the convention", "a1 = first lattice vector") is prose
    m = re.match(r"[A-Za-z_]\w*(?:%\w+)*\s*(?:\([^!]*?\))?\s*=(.*)", body)
    if m:
        rhs = re.sub(r"'[^']*'|\"[^\"]*\"", "", m.group(1).split("!")[0])
        if len(PROSE.findall(body.split("!")[0])) >= 2 or re.search(r"\b[A-Za-z]{2,}\s+[A-Za-z]{2,}\b", rhs):
            return False
    if re.match(r"(?:else\s*)?if\s*\(", body, re.I) and not re.search(r"\)\s*(then\b|\w+\s*=|call\b|cycle|exit|return|stop|&)", body, re.I):
        return False                                   # "if (...)" inside a sentence
    return True


def strip_comments(lines):
    return [l.rstrip() for l in lines if l.strip() and not (l.lstrip().startswith("!") and not l.lstrip().upper().startswith("!$OMP"))]


def process(path, dry, show):
    lines = open(path, errors="replace").read().split("\n")
    out, removed, prev = [], [], False
    for l in lines:
        if is_code_comment(l, prev):
            removed.append(l)
            prev = "&" if l.split("!")[-1].rstrip().endswith("&") else True
        else:
            out.append(l)
            prev = False
    squeezed = []
    for l in out:
        if not l.strip() and squeezed and not squeezed[-1].strip():
            continue
        squeezed.append(l)
    assert strip_comments(lines) == strip_comments(squeezed), f"{path}: code lines changed"
    kept = sum(1 for l in squeezed if l.lstrip().startswith("!") and not l.lstrip().upper().startswith("!$OMP"))
    print(f"{path}: {len(lines)} -> {len(squeezed)} lines; {len(removed)} commented-out code lines removed, {kept} comment lines kept")
    for l in removed[:show]:
        print("   - " + l.strip()[:110])
    if not dry:
        open(path, "w").write("\n".join(squeezed))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("files", nargs="+")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--show", type=int, default=0)
    a = p.parse_args()
    for f in a.files:
        process(f, a.dry_run, a.show)
