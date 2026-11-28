#!/usr/bin/env python3
"""Generate the Fortran API reference of GRABNES from the code repository.

Usage:
    uv run python tools/gen_api_reference.py [--src PATH] [--out PATH]

``--src`` is ``lanczosKuboCode/Src`` of a checkout of the code repository;
``--out`` is the overview page (default ``reference/units.md``). One page per
program unit is written next to it, under ``reference/api/``, and one page per
public procedure under ``reference/api/<unit>/`` (the directory is emptied
first, so units removed upstream disappear). The pages are MyST Markdown for
MyST-MD (front matter ``edit_url: null``: generated pages have no edit link),
and their navigation is written into ``myst.yml`` between the ``BEGIN/END
generated`` markers of ``project.toc`` (each unit with its procedures).

What is shown, and where it comes from:

* Program units, procedures, arguments, results, generic interfaces, derived
  types and module variables: parsed by Ford. Ford runs on a normalized copy in
  a temporary directory where ordinary ``!`` comments are blanked (Ford would
  otherwise mistake them for documentation); ``!!`` lines are kept. The cpp
  template fragments ``MIO/mem_temp.f90`` and ``MIO/MPI/mpi_type.f90`` are
  include files with ``[_TYPE]`` placeholders and are skipped.
* Summaries, details and argument descriptions: the Doxygen tags ``@brief``,
  ``@details``, ``@note``, ``@par``, ``@see``, ``@param`` and ``@return`` (or an
  untagged ``!>`` block) of the comment block right above a
  ``module``/``subroutine``/``function`` statement, copied verbatim. Without a
  Doxygen block, the plain ``!`` comment above the statement is used
  (``parse_comment``): the MIO ``Subroutine: Name`` headers give the description
  and, from their ``INPUT``/``OUTPUT`` lists, the argument descriptions. Ford does
  not read these comments, so they are read from the sources directly. Module
  variables, derived types and their components get the comment directly above
  their declaration.
* Calls and callers: Ford's ``correlate()``, with the names renamed in ``use``
  statements (``use sys, MIO_Kill => SysKill``) mapped back. correlate() drops
  private procedures and the ``!!`` documentation and adds re-exported entities to
  the public lists, so these are read before it runs.
* Input keys: the ``!! - `Key`: description`` item that follows every line
  reading an input key (collected by Ford into the routine's ``doc_list``).
* Grouping: the object list ``OBJS`` of ``Src/Makefile`` (solver), the
  ``MIO/`` and ``math/`` directories (libraries); programs inside the library
  directories are build-time code generators; anything else is not compiled.
* Source links: the definition line in the original file, at the exact commit.

The output depends only on the sources (the date shown is the commit date), so
regenerating an unchanged commit gives an identical page.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

CODE_REPO = "https://github.com/gjung-group/grabnes"
SRC_PREFIX = "lanczosKuboCode/Src"
TEMPLATE_FRAGMENTS = {"mem_temp.f90", "mpi_type.f90"}

# ---------------------------------------------------------------------------
# Source normalization (for Ford)


def strip_comments(text: str) -> str:
    """Blank ordinary ``!`` comments; keep ``!!`` lines (Ford documentation).

    Line numbers are preserved. A ``!`` inside a quoted string does not start
    a comment.
    """
    out = []
    for line in text.splitlines():
        in_str: str | None = None
        cut = -1
        for i, ch in enumerate(line):
            if in_str:
                if ch == in_str:
                    in_str = None
            elif ch in "\"'":
                in_str = ch
            elif ch == "!":
                cut = i
                break
        if cut < 0:
            out.append(line)
        elif line.lstrip().startswith("!!"):
            out.append(line)
        elif not line[:cut].strip():
            out.append("")  # whole-line comment
        else:
            out.append(line[:cut].rstrip())  # inline comment
    return "\n".join(out) + "\n"


def fortran_files(src: Path) -> list[Path]:
    return [
        f for f in sorted(src.rglob("*"))
        if f.is_file() and f.suffix.lower() == ".f90" and f.name not in TEMPLATE_FRAGMENTS
    ]


def normalize_tree(src: Path, dst: Path) -> None:
    for f in fortran_files(src):
        target = dst / f.relative_to(src)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(strip_comments(f.read_text(errors="replace")))


# ---------------------------------------------------------------------------
# Doxygen comments and definition lines, read from the original sources

DEF_RE = re.compile(
    r"^\s*(?!end\b)(?:[\w(),=*:\s]*?\s)?(subroutine|function|module|program)\s+(\w+)\s*(?:\(|$|&|!|result\b)",
    re.IGNORECASE,
)
TAG_RE = re.compile(r"^@(\w+)(?:\[[^\]]*\])?\s*(.*)$")


@dataclass
class Doxy:
    """Documentation of one entity, from a Doxygen block or a plain comment."""
    brief: str = ""
    details: list[str] = field(default_factory=list)  # paragraphs
    notes: list[str] = field(default_factory=list)
    params: dict[str, str] = field(default_factory=dict)  # lower-case name -> text
    returns: str = ""
    pars: list[tuple[str, list[str]]] = field(default_factory=list)  # @par title, items
    see: list[str] = field(default_factory=list)  # @see names
    # INPUT/OUTPUT entries of a plain comment header: (declaration, name, text)
    entries: list[tuple[str, str, str]] = field(default_factory=list)
    banner: str = ""  # name in a "Subroutine: Name" banner of a plain comment

    def __bool__(self) -> bool:
        return bool(self.brief or self.details or self.notes or self.params or self.returns
                    or self.pars or self.see or self.entries)

    def summary(self) -> str:
        """First sentence of the description (for summary tables)."""
        text = self.brief or (self.details[0] if self.details else "")
        m = re.match(r"(.+?[.!?])(\s|$)", text)
        return m.group(1) if m else text


def parse_doxygen(block: list[str]) -> Doxy:
    """Parse ``!>``/``!!`` comment lines (top to bottom) with Doxygen tags.

    A block opened by ``!>`` without tags is a description (its first paragraph
    is the summary). ``!!`` lines alone are Ford documentation of the preceding
    statement and are not read here.
    """
    d = Doxy()
    marks = [ln.strip()[:2] for ln in block if ln.strip().startswith(("!>", "!!"))]
    if not any("@" in ln for ln in block) and "!>" not in marks:
        return d
    target: tuple[str, str] | None = None  # (tag, param name)
    buf: list[str] = []

    def flush():
        nonlocal buf
        lines = [s for s in buf if s]
        text = " ".join(lines).strip()
        buf = []
        if not target or not text:
            return
        tag, name = target
        if tag == "brief":
            d.brief = f"{d.brief} {text}".strip()
        elif tag == "details":
            d.details.append(text)
        elif tag == "note":
            d.notes.append(text)
        elif tag == "param" and name:
            d.params[name.lower()] = f"{d.params.get(name.lower(), '')} {text}".strip()
        elif tag in ("return", "returns", "retval"):
            d.returns = f"{d.returns} {text}".strip()
        elif tag == "see":
            d.see += [n for n in re.split(r"[\s,]+", text) if n]
        elif tag == "par":
            items: list[str] = []
            for s in lines:
                if s.startswith("- ") or not items:
                    items.append(s.removeprefix("- ").strip())
                else:
                    items[-1] = f"{items[-1]} {s}"
            d.pars.append((name, items))

    prev_mark = ""
    for raw in block:
        s = raw.strip()
        if not (s.startswith("!>") or s.startswith("!!")):
            continue
        if re.fullmatch(r"!+", s):
            continue  # a line of '!' (box border)
        mark, s = s[:2], s[2:].strip()
        new_block = mark == "!>" and prev_mark == "!!"  # Doxygen: '!>' after '!!' opens a new block
        prev_mark = mark
        m = TAG_RE.match(s)
        if new_block and not m and s:
            flush()
            target = ("details", "") if d.brief else ("brief", "")
            buf = [s]
            continue
        if m:
            flush()
            tag, rest = m.group(1).lower(), m.group(2)
            if tag == "param":
                name, _, rest = rest.partition(" ")
                target = ("param", name.strip())
            elif tag == "par":
                target = ("par", rest.strip())
                rest = ""
            elif tag in ("brief", "details", "note", "return", "returns", "retval", "see"):
                target = (tag, "")
            else:
                target = None  # @file, @author, ...: not shown
            buf = [rest.strip()]
        elif not s:
            flush()  # paragraph break
            if target and target[0] == "brief":
                target = ("details", "")
        else:
            if target is None and not d.brief and not d.details:
                target = ("brief", "")
            buf.append(s)
    flush()
    return d


BANNER_RE = re.compile(r"^\**\s*(subroutine|function|module|program)\s*:\s*(\w+)[\s*]*$", re.IGNORECASE)
LABEL_RE = re.compile(r"^description\s*:\s*(.*)$", re.IGNORECASE)
SECTION_RE = re.compile(r"^(input|output)\s*-+\s*$", re.IGNORECASE)
ENTRY_RE = re.compile(r"^(\S.*?)\s+:\s+(.*)$")


def parse_comment(block: list[str]) -> Doxy:
    """Parse a plain ``!`` comment (top to bottom) written above a definition.

    Two layouts occur in the sources besides free text: the header of the MIO
    library (a ``Subroutine: Name`` banner between lines of ``*``, the
    description, then ``INPUT ---``/``OUTPUT ---`` lists of ``declaration name :
    text``) and a ``SUBROUTINE: Name`` / ``DESCRIPTION:`` header. The banner,
    the labels and the border lines are dropped; every other word is kept as
    written. Lines are joined into paragraphs at blank comment lines (and after
    a line ending with ``:``); continuation lines of an entry are appended to it.
    """
    d = Doxy()
    if not block or any(ln.strip().startswith("!$") for ln in block):
        return d  # OpenMP directives: code, not a comment
    boxed = re.fullmatch(r"!{3,}", block[0].strip()) is not None
    paras: list[list[str]] = [[]]
    section = ""
    for raw in block:
        c = re.sub(r"^!+", "", raw.strip())
        if boxed:
            c = re.sub(r"\s*!\s*$", "", c)  # right border of a box of '!'
        c = c.strip()
        if re.fullmatch(r"[\s*=!\-]*", c):
            if not section:
                paras.append([])
            continue
        m = BANNER_RE.match(c)
        if m:
            d.banner = m.group(2)
            continue
        m = SECTION_RE.match(c)
        if m:
            section = m.group(1).lower()
            continue
        m = LABEL_RE.match(c)
        if m:
            c = m.group(1).strip()
            if not c:
                continue
        if section:
            m = ENTRY_RE.match(c)
            if m:
                decl = m.group(1).strip()
                d.entries.append((decl, decl.split()[-1], m.group(2).strip()))
            elif d.entries:
                decl, name, text = d.entries[-1]
                d.entries[-1] = (decl, name, f"{text} {c}")
            continue
        if re.fullmatch(r"(subroutine|function)\s*:?", c, re.IGNORECASE):
            continue
        paras[-1].append(c)
        if c.endswith(":"):
            paras.append([])
    texts = [" ".join(p) for p in paras if p]
    if texts:
        d.brief, d.details = texts[0], texts[1:]
    return d


class SourceFile:
    """Original text of one Fortran file: definition lines and Doxygen blocks."""

    def __init__(self, path: Path):
        self.lines = path.read_text(errors="replace").splitlines()
        self.defs: dict[str, list[int]] = defaultdict(list)  # lower name -> 1-based lines
        for i, line in enumerate(self.lines, 1):
            if line.lstrip().startswith(("!", "#")):
                continue
            m = DEF_RE.match(line)
            if m and not re.match(r"^\s*module\s+procedure\b", line, re.IGNORECASE):
                self.defs[m.group(2).lower()].append(i)

    def line_of(self, name: str, after: int = 0) -> int | None:
        for ln in self.defs.get(name.lower(), []):
            if ln > after:
                return ln
        return None

    def comment_above(self, line: int | None) -> list[str]:
        """The comment lines directly above a 1-based line, top to bottom."""
        if not line:
            return []
        i = line - 2
        block: list[str] = []
        while i >= 0 and self.lines[i].strip().startswith("!"):
            block.append(self.lines[i])
            i -= 1
        return block[::-1]

    def doxygen_above(self, line: int | None) -> Doxy:
        return parse_doxygen(self.comment_above(line))

    def doc_above(self, line: int | None, name: str = "") -> Doxy:
        """Documentation of the definition at ``line``: its Doxygen block or,
        without one, the plain comment above it (not used when its banner names
        another entity, as for a comment that describes a generic interface)."""
        block = self.comment_above(line)
        d = parse_doxygen(block)
        if d:
            return d
        plain = [ln for ln in block if not ln.strip().startswith(("!>", "!!")) or re.fullmatch(r"!{3,}", ln.strip())]
        d = parse_comment(plain)
        if d.banner and name and d.banner.lower() != name.lower():
            return Doxy(banner=d.banner)
        return d

    def banner_docs(self) -> dict[str, tuple[int, Doxy]]:
        """Plain comment headers by the name in their banner (lower case)."""
        out: dict[str, tuple[int, Doxy]] = {}
        for i, line in enumerate(self.lines):
            c = re.sub(r"^!+", "", line.strip()) if line.strip().startswith("!") else ""
            m = BANNER_RE.match(c.strip())
            if not m:
                continue
            j = i
            while j < len(self.lines) and self.lines[j].strip().startswith("!"):
                j += 1
            d = parse_comment(self.comment_above(j + 1))
            out.setdefault(m.group(2).lower(), (j + 1, d))
        return out

    def use_statements(self) -> list[str]:
        """``use`` statements, continuation lines joined, comments removed."""
        out, cur = [], None
        for line in self.lines:
            code_part = line.split("!")[0].rstrip()
            if cur is not None:
                cur += " " + code_part.strip().lstrip("&")
            elif re.match(r"^\s*use\b", code_part, re.IGNORECASE):
                cur = code_part.strip()
            else:
                continue
            if cur.endswith("&"):
                cur = cur[:-1]
            else:
                out.append(cur)
                cur = None
        return out

    def decl_line(self, name: str, start: int = 0, kind: str = "variable") -> int | None:
        """1-based line declaring a module variable (``:: name``) or a derived
        type (``type ... name``), searched from ``start`` up to ``contains``."""
        if kind == "type":
            pat = re.compile(rf"^\s*type\b(?!\s*\().*?(::)?\s*{re.escape(name)}\s*(!.*)?$", re.IGNORECASE)
        else:
            pat = re.compile(rf"::.*?(?<![\w%]){re.escape(name)}\b(?!\s*%)", re.IGNORECASE)
        for i in range(start, len(self.lines)):
            s = self.lines[i]
            if re.match(r"^\s*contains\b", s, re.IGNORECASE):
                return None
            if s.lstrip().startswith(("!", "#")) or re.match(
                    r"^\s*(public|private|protected|save|use|namelist|volatile|target)\b", s, re.IGNORECASE):
                continue
            if pat.search(s.split("!")[0] if kind != "type" else s):
                return i + 1
        return None

    def header_doxygen(self) -> Doxy:
        """Doxygen block at the top of the file (e.g. ``!> @file`` ... ``@brief``)."""
        block = []
        for line in self.lines:
            if line.strip().startswith("!"):
                block.append(line)
            elif line.strip():
                break
            elif block:
                break
        return parse_doxygen(block)


# ---------------------------------------------------------------------------
# Input-key items (``!! - `Key`: description``)

KEY_CACHE: dict[int, list[tuple[str, str]]] = {}  # id(procedure) -> items, read before correlate()
KEY_ITEM = re.compile(r"^-\s*`([^`]+)`:\s*(.*)$")


def key_items(p) -> list[tuple[str, str]]:
    """Input-key descriptions in a procedure's ``!!`` comments (Ford ``doc_list``).

    Upstream writes one Markdown item after every line that reads a key; its
    continuation lines are indented. Other ``!!`` lines are ignored.
    """
    if id(p) in KEY_CACHE:
        return KEY_CACHE[id(p)]
    items: list[list[str]] = []
    current: list[str] | None = None
    for raw in getattr(p, "doc_list", None) or []:
        m = KEY_ITEM.match(raw.strip())
        if m:
            current = [m.group(1), m.group(2).strip()]
            items.append(current)
        elif current is not None and raw.startswith("  ") and raw.strip() and not raw.strip().startswith("!"):
            # (a line starting with '!' is commented-out code written as '!!   ! ...')
            current[1] = f"{current[1]} {raw.strip()}".strip()
        else:
            current = None
    return [(k, v) for k, v in items]


# ---------------------------------------------------------------------------
# Markdown helpers


def md_text(text: str) -> str:
    """Escape what MyST would read as markup in plain text (``*``/``_`` emphasis,
    ``$`` math, ``[``/``]`` links, ``<`` HTML), outside inline code spans."""
    parts = re.split(r"(`[^`]*`)", text)
    return "".join(
        part if part.startswith("`") else re.sub(r"([*_$\[\]<])", r"\\\1", part)
        for part in parts
    )


def cell(text: str) -> str:
    """Plain text for a table cell."""
    return md_text(text).replace("|", "\\|") if text else ""


def code(text: str) -> str:
    """Inline code that survives backticks and table pipes."""
    text = str(text).replace("|", "\\|")
    return f"`` {text} ``" if "`" in text else f"`{text}`"


# Generated pages have no "Edit this page" link (MyST-MD front matter).
FRONTMATTER = "---\nedit_url: null\n---\n\n"

# The navigation of the generated pages is written into myst.yml between these
# lines (MyST-MD reads the table of contents from myst.yml only).
TOC_BEGIN = "            # BEGIN generated by tools/gen_api_reference.py: do not edit by hand\n"
TOC_END = "            # END generated by tools/gen_api_reference.py\n"


def write_toc(myst_yml: Path, entries: list[tuple[str, list[str]]]) -> bool:
    """Replace the generated part of the table of contents in myst.yml:
    every unit page, with its procedure pages as children."""
    text = myst_yml.read_text()
    if TOC_BEGIN not in text or TOC_END not in text:
        return False
    head, rest = text.split(TOC_BEGIN, 1)
    _, tail = rest.split(TOC_END, 1)
    lines = []
    for unit, procs in entries:
        lines.append(f"            - file: {unit}\n")
        if procs:
            lines.append("              children:\n")
            lines += [f"                - file: {p}\n" for p in procs]
    myst_yml.write_text(head + TOC_BEGIN + "".join(lines) + TOC_END + tail)
    return True


def table(header: list[str], rows: list[list[str]]) -> list[str]:
    """A ``list-table`` with a header row."""
    out = ["```{list-table}", ":header-rows: 1", ""]
    for r in [header] + rows:
        out.append(f"* - {r[0] or ' '}")
        out += [f"  - {c or ' '}" for c in r[1:]]
    return out + ["```", ""]


def ident(name: str) -> str:
    return re.sub(r"[^a-z0-9_]+", "-", name.lower())


# ---------------------------------------------------------------------------
# Ford entities -> text


def proto_name(v) -> str:
    p = getattr(v, "proto", None)
    if not p:
        return ""
    first = p[0] if isinstance(p, (list, tuple)) else p
    return first if isinstance(first, str) else getattr(first, "name", str(first))


def type_spec(v) -> str:
    """Fortran type specification, e.g. ``real(dp)``, ``character(len=*)``."""
    t = (getattr(v, "vartype", "") or "").lower()
    if t in ("type", "class", "procedure"):
        return f"{t}({proto_name(v) or '*'})"
    if t == "character":
        ln = getattr(v, "strlen", None)
        kind = getattr(v, "kind", None)
        inner = ", ".join(x for x in [f"len={ln}" if ln else "", f"kind={kind}" if kind else ""] if x)
        return f"character({inner})" if inner else "character"
    kind = getattr(v, "kind", None)
    return f"{t}({kind})" if kind else t


def var_attribs(v, with_intent: bool = False) -> list[str]:
    attrs = []
    dim = getattr(v, "dimension", "") or ""
    if dim:
        attrs.append(f"dimension{dim}" if dim.startswith("(") else f"dimension({dim})")
    if getattr(v, "parameter", False):
        attrs.append("parameter")
    if with_intent and getattr(v, "intent", ""):
        attrs.append(f"intent({v.intent})")
    if getattr(v, "optional", False):
        attrs.append("optional")
    for a in getattr(v, "attribs", None) or []:
        a = str(a).strip()
        if a and a.lower() not in {x.lower() for x in attrs} and not a.lower().startswith(("intent", "dimension")):
            attrs.append(a)
    return attrs


def declaration(v, with_intent: bool = False) -> str:
    """Declaration as in the source, e.g. ``real(dp), dimension(3), intent(in)``."""
    return code(", ".join([type_spec(v)] + var_attribs(v, with_intent)))


def is_function(p) -> bool:
    return type(p).__name__ == "FortranFunction"


def signature(p) -> str:
    args = ", ".join(a.name for a in (getattr(p, "args", None) or []))
    prefix = " ".join(a for a in (getattr(p, "attribs", None) or []) if a.lower() in
                      ("pure", "elemental", "recursive", "impure", "module"))
    prefix = f"{prefix} " if prefix else ""
    if is_function(p):
        r = getattr(p, "retvar", None)
        rname = getattr(r, "name", p.name) if r is not None else p.name
        rtype = type_spec(r) if r is not None and getattr(r, "vartype", None) else ""
        result = f" result({rname})" if rname.lower() != p.name.lower() else ""
        typed = f"{rtype} " if rtype and not result else ""
        return f"{prefix}{typed}function {p.name}({args}){result}"
    return f"{prefix}subroutine {p.name}({args})"


def own(unit, entity) -> bool:
    """True if ``entity`` is defined in ``unit`` (Ford's public lists also hold
    what a module re-exports from the modules it uses)."""
    return getattr(entity, "parent", None) is unit


PROC_CACHE: dict[int, list] = {}  # id(entity) -> its procedures, read before correlate()


def procs_of(entity) -> list:
    if id(entity) in PROC_CACHE:
        return PROC_CACHE[id(entity)]
    return list(getattr(entity, "subroutines", None) or []) + list(getattr(entity, "functions", None) or [])


def all_procs(entity, parents=()):
    """(procedure, parent names) for every procedure, nested ones included."""
    for p in procs_of(entity):
        yield p, parents
        yield from all_procs(p, parents + (p.name,))


def use_names(d) -> list[str]:
    if isinstance(d, str):
        head = d.split(",")[0].split("(")[0].strip()
        return [] if not head or head.lower().startswith("only") else [head]
    if isinstance(d, (list, tuple, set)):
        return [n for x in d for n in use_names(x)]
    return [getattr(d, "name", str(d)).strip()]


# ---------------------------------------------------------------------------
# Model


@dataclass
class Unit:
    kind: str  # "program" | "module"
    obj: object
    path: str  # relative to Src, posix
    group: str
    label: str = ""
    page: str = ""
    uses: list[str] = field(default_factory=list)
    used_by: list[str] = field(default_factory=list)
    pub: dict = field(default_factory=dict)  # public lists of the unit, as parsed

    @property
    def name(self) -> str:
        return self.obj.name


GROUPS = [
    ("solver", "Solver", "The `grabnes` program and the modules compiled into it (the object list `OBJS` of `Src/Makefile`)."),
    ("mio", "MIO library", "Memory, input and output library (`MIO/`, built as `libMIO.a`): input parsing, messages, files, timers, OpenMP and MPI support."),
    ("math", "Math library", "Numerical routines (`math/`, built as `libmath.a`)."),
    ("generator", "Build-time code generators", "Programs compiled and run by the makefiles of the libraries to generate Fortran sources before the build."),
    ("unbuilt", "Sources not in the solver build", "Files in `Src/` that are not in the object list of `Src/Makefile` and are therefore not compiled into `grabnes`."),
]


def makefile_objects(src: Path) -> set[str]:
    mk = src / "Makefile"
    if not mk.exists():
        return set()
    text = mk.read_text(errors="replace").replace("\\\n", " ")
    m = re.search(r"^OBJS\s*:?=\s*(.*)$", text, re.MULTILINE)
    return {Path(o).stem.lower() for o in m.group(1).split() if o.endswith(".o")} if m else set()


def classify(path: str, kind: str, objs: set[str]) -> str:
    top = path.split("/")[0]
    stem = Path(path).stem.lower()
    if top in ("MIO", "math"):
        return "generator" if kind == "program" else top.lower()
    if stem in objs or (kind == "program" and stem == "grabnes"):
        return "solver"
    return "unbuilt"


# ---------------------------------------------------------------------------
# Rendering


class Renderer:
    def __init__(self, units: list[Unit], sources: dict[str, SourceFile], commit: str, date: str, ford_version: str):
        self.units = units
        self.sources = sources
        self.commit = commit
        self.date = date
        self.ford_version = ford_version
        self.by_name: dict[str, Unit] = {}
        for u in units:  # first unit of a name wins (the built one sorts first)
            self.by_name.setdefault(u.name.lower(), u)

    # -- links ------------------------------------------------------------
    def src_url(self, path: str, line: int | None = None) -> str:
        if not self.commit:
            return ""
        url = f"{CODE_REPO}/blob/{self.commit}/{SRC_PREFIX}/{path}"
        return f"{url}#L{line}" if line else url

    def src_link(self, path: str, line: int | None = None, text: str | None = None) -> str:
        label = text or (f"{SRC_PREFIX}/{path}" + (f", line {line}" if line else ""))
        url = self.src_url(path, line)
        return f"[{code(label) if text is None else label}]({url})" if url else code(label)

    def unit_ref(self, name: str) -> str:
        u = self.by_name.get(name.lower())
        return f"{{ref}}`{name} <{u.label}>`" if u else code(name)

    # -- overview ---------------------------------------------------------
    def overview(self) -> str:
        stats = self.stats()
        short = self.commit[:7] if self.commit else "unknown"
        commit_link = f"[`{short}`]({CODE_REPO}/tree/{self.commit})" if self.commit else "`unknown`"
        out = [
            "(api-units)=",
            "# Fortran unit reference",
            "",
            "Reference documentation of every program unit in `lanczosKuboCode/Src/` of the",
            f"[GRABNES code repository]({CODE_REPO}). Each unit has a page with its",
            "generic interfaces, derived types, module variables and the units it uses",
            "and is used by; each public procedure has its own page with its signature,",
            "parameters, return value, the input keys it reads, and the procedures it",
            "calls and is called by.",
            "",
            ":::{admonition} Provenance",
            ":class: note",
            "",
            f"Generated by `tools/gen_api_reference.py` with Ford {self.ford_version} from commit",
            f"{commit_link} of the code repository ({self.date}). Do not edit by hand; the page is",
            "regenerated whenever the synchronization with the code repository brings new sources.",
            "Descriptions are quoted verbatim from the sources: the Doxygen block (`@brief`, `@param`,",
            "...) above a definition or, without one, the plain comment directly above it (its",
            "`Subroutine: Name` banner, labels and border lines dropped; its `INPUT`/`OUTPUT` lists",
            "give the parameter descriptions), and the Ford (`!!`) comments for the input keys.",
            "An entry without a description has none in the source. Calls and callers come from",
            "Ford's call analysis, with the names renamed in `use` statements resolved.",
            ":::",
            "",
            "::::{grid} 2 2 4 4",
            "",
        ]
        for value, what in stats:
            out += [f":::{{card}} {value}", "", what, ":::", ""]
        out += ["::::", ""]

        for gid, title, blurb in GROUPS:
            members = [u for u in self.units if u.group == gid]
            if not members:
                continue
            out += [f"(api-group-{gid})=", f"## {title}", "", blurb, ""]
            rows = []
            for u in members:
                pubs = [p for p in self.public_procs(u)]
                nkeys = sum(len(key_items(p)) for p, _ in self.all_procs(u.obj))
                rows.append([
                    self.unit_ref(u.name) if self.by_name.get(u.name.lower()) is u else f"{{ref}}`{u.name} <{u.label}>`",
                    u.kind,
                    self.src_link(u.path, text=code(u.path)),
                    str(len(pubs)) if pubs else "—",
                    str(nkeys) if nkeys else "—",
                ])
            out += table(["Unit", "Kind", "Source", "Public procedures", "Input-key reads"], rows)

        out += [
            "```{note}",
            "The cpp template fragments `MIO/mem_temp.f90` and `MIO/MPI/mpi_type.f90` are include",
            "files with `[_TYPE]` placeholders, expanded by the code generators above; they are",
            "not program units and are not listed.",
            "```",
            "",
        ]
        return FRONTMATTER + "\n".join(out)

    def stats(self) -> list[tuple[str, str]]:
        n_units = len(self.units)
        n_pub = sum(len(self.public_procs(u)) for u in self.units)
        n_if = sum(len(self.interfaces(u)) for u in self.units)
        n_keys = sum(len(key_items(p)) for u in self.units for p, _ in self.all_procs(u.obj))
        n_doc = sum(1 for u in self.units for p in self.public_procs(u) if self.proc_doc(u, p))
        return [
            (str(n_units), "program units"),
            (str(n_pub), f"public procedures ({n_doc} with a description in the source)"),
            (str(n_if), "generic interfaces"),
            (str(n_keys), "documented input-key reads"),
        ]

    # -- unit queries -----------------------------------------------------
    def unit_doxy(self, u: Unit) -> Doxy:
        sf = self.sources[u.path]
        d = sf.doc_above(sf.line_of(u.name), u.name)
        return d if d else sf.header_doxygen()

    def public_procs(self, u: Unit) -> list:
        if u.kind == "program":
            return []
        return sorted(
            (p for p in (u.pub["pub_procs"]).values()
             if type(p).__name__ in ("FortranSubroutine", "FortranFunction") and own(u.obj, p)),
            key=lambda p: p.name.lower(),
        )

    def interfaces(self, u: Unit) -> list:
        if u.kind == "program":
            return []
        return sorted(
            (p for p in (u.pub["pub_procs"]).values()
             if type(p).__name__ == "FortranInterface" and own(u.obj, p)),
            key=lambda p: p.name.lower(),
        )

    all_procs = staticmethod(lambda entity, parents=(): all_procs(entity, parents))

    def proc_doc(self, u: Unit, p) -> Doxy:
        sf = self.sources[u.path]
        return sf.doc_above(sf.line_of(p.name), p.name)

    def proc_label(self, u: Unit, p) -> str:
        return f"{u.label}-{ident(p.name)}"

    def proc_page(self, u: Unit, p) -> str:
        return f"{u.page}/{ident(p.name)}"

    def proc_ref(self, u: Unit, p) -> str:
        return f"{{ref}}`{p.name} <{self.proc_label(u, p)}>`"

    def iface_label(self, u: Unit, i) -> str:
        return f"{u.label}-interface-{ident(i.name)}"

    # -- call graph -------------------------------------------------------
    def build_graph(self) -> None:
        """Procedures called by every public procedure (with its internal
        procedures) and by every program, from Ford's call analysis; the
        reverse gives "called by". Only procedures of the sources are kept."""
        self.pub_index: dict[str, list[tuple[Unit, object]]] = defaultdict(list)
        self.iface_index: dict[str, list[tuple[Unit, object]]] = defaultdict(list)
        self.priv_index: dict[str, list[Unit]] = defaultdict(list)
        self.public: set[int] = set()
        for u in self.units:
            for p in self.public_procs(u):
                self.pub_index[p.name.lower()].append((u, p))
                self.public.add(id(p))
            for i in self.interfaces(u):
                self.iface_index[i.name.lower()].append((u, i))
            for p, parents in all_procs(u.obj):
                if not parents:
                    self.priv_index[p.name.lower()].append(u)
        # local names given by a rename in a use statement (``use sys, MIO_Kill => SysKill``)
        self.aliases: dict[str, str] = {}
        self.alias_names: dict[str, str] = {}  # lower -> as written
        for sf in self.sources.values():
            for stmt in sf.use_statements():
                for alias, target in re.findall(r"(\w+)\s*=>\s*(\w+)", stmt):
                    self.aliases.setdefault(alias.lower(), target.lower())
                    self.alias_names.setdefault(alias.lower(), alias)
        self.calls: dict[tuple[str, str], list] = {}
        self.called_by: dict[tuple[str, str], list] = defaultdict(list)
        for u in self.units:
            callers = [(p, [p] + [q for q, _ in all_procs(p)]) for p in procs_of(u.obj)]
            if u.kind == "program":
                callers.append((None, [u.obj]))
            for p, bodies in callers:
                targets, seen = [], set()
                for body in bodies:
                    for c in getattr(body, "calls", None) or []:
                        t = self.resolve(c, u)
                        if t and t[2] not in seen and not (p is not None and t[2] == p.name.lower() and t[1] is u):
                            seen.add(t[2])
                            targets.append(t)
                key = (u.label, p.name.lower() if p is not None else "")
                self.calls[key] = sorted(targets, key=lambda t: t[2])
                for kind, tu, tname, tobj, _ in targets:
                    if kind == "public":
                        entry = (u, p)
                        if entry not in self.called_by[(tu.label, tname)]:
                            self.called_by[(tu.label, tname)].append(entry)

    @staticmethod
    def owner_module(obj) -> str:
        o = getattr(obj, "parent", None)
        while o is not None and type(o).__name__ not in ("FortranModule", "FortranProgram"):
            o = getattr(o, "parent", None)
        return o.name.lower() if o is not None else ""

    def resolve(self, c, caller: Unit):
        """(kind, unit, lower name, object) of a called procedure, or None."""
        name = (c if isinstance(c, str) else getattr(c, "name", "")).lower()
        if not name:
            return None
        shown = c if isinstance(c, str) else getattr(c, "name", "")
        if isinstance(c, str) and name not in self.pub_index and name not in self.iface_index:
            shown = self.alias_names.get(name, shown)
            name = self.aliases.get(name, name)
        owner = "" if isinstance(c, str) else self.owner_module(c)

        def pick(cands):
            for want in (owner, caller.name.lower()):
                for cu, obj in cands:
                    if want and cu.name.lower() == want:
                        return cu, obj
            return cands[0] if len(cands) == 1 else None

        if name in self.pub_index:
            hit = pick(self.pub_index[name])
            if hit:
                return ("public", hit[0], name, hit[1], shown)
        if name in self.iface_index:
            hit = pick(self.iface_index[name])
            if hit:
                return ("interface", hit[0], name, hit[1], shown)
        if not isinstance(c, str) and name in self.priv_index:
            us = self.priv_index[name]
            cu = next((x for x in us if x.name.lower() in (owner, caller.name.lower())), us[0] if len(us) == 1 else None)
            if cu:
                return ("private", cu, name, c, shown)
        return None

    def target_ref(self, t) -> str:
        kind, tu, name, obj, shown = t
        # the name used in the call (a renamed import, e.g. MIO_Kill for SysKill)
        alias = shown if shown and shown.lower() != obj.name.lower() else ""
        if kind == "public":
            return self.proc_ref(tu, obj) if not alias else \
                f"{{ref}}`{alias} <{self.proc_label(tu, obj)}>` (= {code(obj.name)})"
        if kind == "interface":
            label = self.iface_label(tu, obj)
            return f"{{ref}}`{obj.name} <{label}>`" if not alias else \
                f"{{ref}}`{alias} <{label}>` (= {code(obj.name)})"
        return self.private_ref(tu, obj)

    def private_ref(self, u: Unit, p) -> str:
        """A procedure without a page: its name, linked to its definition."""
        line = self.sources[u.path].line_of(p.name)
        name = self.src_link(u.path, line, text=md_text(p.name)) if line else md_text(p.name)
        return f"{name} ({{ref}}`{u.name} <{u.label}>`)"

    def caller_ref(self, entry) -> str:
        cu, p = entry
        if p is None:
            return f"{{ref}}`{cu.name} <{cu.label}>` (program)"
        if id(p) in self.public:
            return self.proc_ref(cu, p)
        return self.private_ref(cu, p)

    # -- shared blocks ----------------------------------------------------
    @staticmethod
    def description(d: Doxy) -> list[str]:
        out = []
        if d.brief:
            out += [md_text(d.brief), ""]
        for para in d.details:
            out += [md_text(para), ""]
        for note in d.notes:
            out += ["```{note}", md_text(note), "```", ""]
        for title, items in d.pars:
            out += [f"**{md_text(title) or 'Note'}**", ""]
            out += [f"- {md_text(i)}" for i in items] + [""]
        return out

    @staticmethod
    def deflist(rows: list[tuple[str, str]]) -> list[str]:
        out = []
        for term, text in rows:
            out += [term, f": {text or '—'}", ""]
        return out

    def var_doc(self, sf: SourceFile, line: int | None) -> str:
        d = sf.doc_above(line)
        return " ".join([d.brief] + d.details) if d and not d.entries else ""

    # -- unit page --------------------------------------------------------
    def unit_page(self, u: Unit) -> str:
        sf = self.sources[u.path]
        uline = sf.line_of(u.name)
        doxy = self.unit_doxy(u)
        pubs = self.public_procs(u)
        ifaces = self.interfaces(u)
        types = sorted((t for t in (u.pub["pub_types"]).values() if own(u.obj, t)),
                       key=lambda t: t.name.lower()) if u.kind == "module" else []
        variables = sorted(
            (v for v in (u.pub["pub_vars"]).values()
             if not type(v).__name__.startswith("FortranNamelist") and own(u.obj, v)),
            key=lambda v: v.name.lower(),
        ) if u.kind == "module" else []

        dup = sum(1 for v in self.units if v.name.lower() == u.name.lower()) > 1
        title = f"{code(u.name)} ({code(Path(u.path).name)})" if dup else code(u.name)
        out = [f"({u.label})=", f"# {title}", ""]
        out += self.description(doxy)

        group_title = next(t for g, t, _ in GROUPS if g == u.group)
        nkeys = sum(len(key_items(p)) for p, _ in self.all_procs(u.obj))
        facts = [
            ["Kind", u.kind.capitalize()],
            ["Source", self.src_link(u.path, uline, text=code(f"{SRC_PREFIX}/{u.path}"))],
            ["Group", f"{{ref}}`{group_title} <api-group-{u.group}>`"],
            ["Uses", ", ".join(self.unit_ref(n) for n in u.uses) or "—"],
            ["Used by", ", ".join(self.unit_ref(n) for n in u.used_by) or "—"],
        ]
        if u.kind == "module":
            facts.append(["Public interface", ", ".join(x for x in [
                f"{len(pubs)} procedure{'s' if len(pubs) != 1 else ''}" if pubs else "",
                f"{len(ifaces)} generic interface{'s' if len(ifaces) != 1 else ''}" if ifaces else "",
                f"{len(types)} derived type{'s' if len(types) != 1 else ''}" if types else "",
                f"{len(variables)} variable{'s' if len(variables) != 1 else ''}" if variables else "",
            ] if x) or "none"])
        if nkeys:
            facts.append(["Input keys", f"{nkeys} reads"])
        out += [f"- **{k}:** {v}" for k, v in facts] + [""]

        if u.kind == "program":
            calls = self.calls.get((u.label, ""), [])
            if calls:
                out += ["## Calls", "", ", ".join(self.target_ref(t) for t in calls), ""]

        if pubs:
            out += ["## Procedures", ""]
            rows = [[self.proc_ref(u, p), "function" if is_function(p) else "subroutine",
                     cell(self.proc_doc(u, p).summary()) or "—"] for p in pubs]
            out += table(["Procedure", "Kind", "Summary"], rows)

        if ifaces:
            banners = sf.banner_docs()
            out += ["## Generic interfaces", ""]
            pub_names = {p.name.lower(): p for p in pubs}
            for i in ifaces:
                specifics = [getattr(m, "name", str(m)) for m in (getattr(i, "modprocs", None) or [])]
                specifics += [r.name for r in (getattr(i, "routines", None) or []) if getattr(r, "name", None)
                              and r.name not in specifics]
                links = [self.proc_ref(u, pub_names[s.lower()]) if s.lower() in pub_names else code(s) for s in specifics]
                iline = sf.line_of(i.name) or next((n for n, s in enumerate(sf.lines, 1) if re.match(
                    rf"^\s*interface\s+{re.escape(i.name)}\b", s, re.IGNORECASE)), None)
                out += [f"({self.iface_label(u, i)})=", f"### {code(i.name)}", ""]
                d = banners.get(i.name.lower(), (0, Doxy()))[1]
                out += self.description(d)
                if d.entries:
                    out += ["**Arguments (from the source comment)**", ""]
                    out += self.deflist([(code(decl), md_text(t)) for decl, _, t in d.entries])
                out += [f"Specific procedures: {', '.join(links) or '—'}.", ""]
                if iline:
                    out += [f"Defined in {self.src_link(u.path, iline)}.", ""]

        if types:
            out += ["## Derived types", ""]
            for t in types:
                tline = sf.decl_line(t.name, uline or 0, "type") or sf.line_of(t.name)
                out += [f"({u.label}-type-{ident(t.name)})=", f"### {code(t.name)}", ""]
                tdoc = self.var_doc(sf, tline)
                if tdoc:
                    out += [md_text(tdoc), ""]
                comps = list(getattr(t, "variables", None) or [])
                if comps:
                    out += table(["Component", "Declaration", "Default", "Description"], [
                        [code(c.name), declaration(c),
                         code(c.initial) if getattr(c, "initial", None) else "—",
                         cell(self.var_doc(sf, sf.decl_line(c.name, tline or 0))) or "—"] for c in comps])
                out += [f"Defined in {self.src_link(u.path, tline)}." if tline else "", ""]

        if variables:
            out += ["## Module variables", ""]
            out += table(["Name", "Declaration", "Initial value", "Description"], [
                [code(v.name), declaration(v),
                 code(v.initial) if getattr(v, "initial", None) else "—",
                 cell(self.var_doc(sf, sf.decl_line(v.name, uline or 0))) or "—"] for v in variables])

        anchors = {p.name.lower() for p in pubs}
        others = [(p, parents) for p, parents in self.all_procs(u.obj)
                  if key_items(p) and not (not parents and p.name.lower() in anchors)
                  and not (parents and parents[0].lower() in anchors)]
        if others:
            out += [
                "## Input keys read by private procedures", "",
                "Procedures that are not part of the public interface but read input keys.", "",
            ]
            for p, parents in others:
                out += self.key_block(u, sf, p, parents)
        return FRONTMATTER + "\n".join(out).rstrip() + "\n"

    # -- procedure page ---------------------------------------------------
    def procedure_page(self, u: Unit, p) -> str:
        sf = self.sources[u.path]
        line = sf.line_of(p.name)
        d = self.proc_doc(u, p)
        kind = "function" if is_function(p) else "subroutine"
        out = [
            f"({self.proc_label(u, p)})=",
            f"# {code(p.name)}",
            "",
            f"```{{index}} single: {p.name} ({kind} in {u.name})", "```",
            "",
            "```fortran", signature(p), "```", "",
            f"- **Module:** {{ref}}`{u.name} <{u.label}>`",
            f"- **Kind:** {kind.capitalize()}",
        ]
        if line:
            out.append(f"- **Source:** {self.src_link(u.path, line)}")
        out.append("")
        out += self.description(d)

        args = list(getattr(p, "args", None) or [])
        r = getattr(p, "retvar", None) if is_function(p) else None
        entries = {n.lower(): (decl, t) for decl, n, t in d.entries}
        used: set[str] = set()
        if args:
            out += ["## Parameters", ""]
            rows = []
            for a in args:
                text = d.params.get(a.name.lower(), "")
                if not text and a.name.lower() in entries:
                    text = entries[a.name.lower()][1]
                    used.add(a.name.lower())
                rows.append((f"{code(a.name)} : {declaration(a, with_intent=True)}", md_text(text)))
            out += self.deflist(rows)
        if r is not None:
            text = d.returns
            if not text and r.name.lower() in entries:
                text = entries[r.name.lower()][1]
                used.add(r.name.lower())
            attrs = ", ".join([type_spec(r)] + var_attribs(r))
            out += ["## Returns", ""]
            out += self.deflist([(f"{code(r.name)} : {code(attrs)}", md_text(text))])
        extra = [(decl, n, t) for decl, n, t in d.entries if n.lower() not in used]
        if extra:
            out += ["## Also described in the source comment", ""]
            out += self.deflist([(code(decl), md_text(t)) for decl, _, t in extra])

        own = key_items(p)
        nested = [(q, parents) for q, parents in self.all_procs(p, (p.name,)) if key_items(q)]
        if own or nested:
            n = len(own) + sum(len(key_items(q)) for q, _ in nested)
            out += [f"## Input keys read ({n})", ""]
            if own:
                out += table(["Key", "Description"], [[code(k), cell(v)] for k, v in own])
            for q, parents in nested:
                out += self.key_block(u, sf, q, parents, after=line or 0)

        see = []
        for name in d.see:
            hits = self.pub_index.get(name.lower(), [])
            hit = next(((cu, o) for cu, o in hits if cu is u), hits[0] if hits else None)
            see.append(self.proc_ref(*hit) if hit else code(name))
        if see:
            out += ["## See also", "", ", ".join(see), ""]
        calls = self.calls.get((u.label, p.name.lower()), [])
        if calls:
            out += ["## Calls", "", ", ".join(self.target_ref(t) for t in calls), ""]
        callers = sorted(self.called_by.get((u.label, p.name.lower()), []),
                         key=lambda e: ((e[1].name if e[1] is not None else e[0].name).lower(), e[0].name.lower()))
        if callers:
            out += ["## Called by", "", ", ".join(self.caller_ref(e) for e in callers), ""]
        return FRONTMATTER + "\n".join(out).rstrip() + "\n"

    def key_block(self, u: Unit, sf: SourceFile, p, parents, after: int = 0) -> list[str]:
        line = sf.line_of(p.name, after)
        where = f"internal procedure of {code(parents[-1])}" if parents else (getattr(p, "permission", "") or "private")
        src = f", {self.src_link(u.path, line, text='source')}" if line else ""
        out = [f"**{code(p.name)}** ({where}{src})", ""]
        out += table(["Key", "Description"], [[code(k), cell(v)] for k, v in key_items(p)])
        return out


# ---------------------------------------------------------------------------


def git(src: Path, *args: str) -> str:
    try:
        return subprocess.run(["git", "-C", str(src), *args], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default="../lanczosKuboCode/Src")
    ap.add_argument("--out", default="reference/units.md")
    a = ap.parse_args()

    src = Path(a.src).resolve()
    if not src.is_dir():
        sys.exit(f"source directory not found: {src}")
    import ford
    from ford import ProjectSettings
    from ford.fortran_project import Project

    commit = git(src, "rev-parse", "HEAD")
    date = git(src, "log", "-1", "--format=%cs") or "unknown date"
    sources = {f.relative_to(src).as_posix(): SourceFile(f) for f in fortran_files(src)}
    objs = makefile_objects(src)

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        normalize_tree(src, tmp)
        project = Project(ProjectSettings(src_dir=[tmp], extensions=["f90"], preprocessor="",
                                          display=["public", "protected", "private"]))

        units: list[Unit] = []
        for kind, objs_ in (("program", project.programs), ("module", project.modules)):
            for o in objs_:
                path = Path(o.source_file.path).relative_to(tmp).as_posix()
                units.append(Unit(kind, o, path, classify(path, kind, objs)))

        # Ford's correlate() resolves the calls between procedures, but it also
        # consumes the ``!!`` documentation and rewrites the public lists with
        # what the modules re-export: keep both as parsed.
        for u in units:
            u.pub = {k: dict(getattr(u.obj, k, None) or {}) for k in ("pub_procs", "pub_types", "pub_vars")}
            PROC_CACHE[id(u.obj)] = procs_of(u.obj)
            for p, _ in all_procs(u.obj):
                PROC_CACHE[id(p)] = procs_of(p)
                KEY_CACHE[id(p)] = key_items(p)
        project.correlate()

        order = {g: i for i, (g, _, _) in enumerate(GROUPS)}
        units.sort(key=lambda u: (order[u.group], u.kind != "program", u.name.lower(), u.path))
        out = Path(a.out)
        api_dir = out.parent / "api"
        seen: set[str] = set()
        for u in units:
            label = f"api-{u.name.lower()}"
            if label in seen:
                label = f"api-{u.name.lower()}-{Path(u.path).stem.lower()}"
            seen.add(label)
            u.label = label
            u.page = f"api/{label.removeprefix('api-')}"
            # use statements of the unit and of every procedure in it
            entities = [u.obj] + [p for p, _ in all_procs(u.obj)]
            u.uses = sorted({n for e in entities for d in (getattr(e, "uses", None) or []) for n in use_names(d)
                             if n and n.lower() != u.name.lower()}, key=str.lower)
        names = {u.name.lower() for u in units}
        for u in units:
            for n in u.uses:
                if n.lower() in names:
                    for v in units:
                        if v.name.lower() == n.lower() and u.name not in v.used_by:
                            v.used_by.append(u.name)
                            break
        for u in units:
            u.used_by.sort(key=str.lower)

        r = Renderer(units, sources, commit, date, ford.__version__)
        if api_dir.exists():
            shutil.rmtree(api_dir)
        api_dir.mkdir(parents=True)
        r.build_graph()
        n_procs = 0
        toc: list[tuple[str, list[str]]] = []
        for u in units:
            (out.parent / f"{u.page}.md").write_text(r.unit_page(u))
            toc.append((f"{u.page}.md", []))
            for p in r.public_procs(u):
                page = out.parent / f"{r.proc_page(u, p)}.md"
                page.parent.mkdir(parents=True, exist_ok=True)
                page.write_text(r.procedure_page(u, p))
                toc[-1][1].append(f"{r.proc_page(u, p)}.md")
                n_procs += 1
        out.write_text(r.overview())

    # myst.yml lists the pages relative to the repository root
    root = Path(__file__).resolve().parent.parent
    base = out.parent.resolve()
    if base.is_relative_to(root) and write_toc(root / "myst.yml", [
            ((base / unit).relative_to(root).as_posix(), [(base / q).relative_to(root).as_posix() for q in procs])
            for unit, procs in toc]):
        print("updated the API reference navigation in myst.yml")
    else:
        print("myst.yml not updated (output outside the repository, or no generated-toc markers)")

    print(f"wrote {out}, {len(units)} unit pages and {n_procs} procedure pages in {api_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
