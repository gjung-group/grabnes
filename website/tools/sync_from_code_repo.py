#!/usr/bin/env python3
"""Keep the generated and derived parts of this site in step with the GRABNES code.

The source is the HEAD commit of the code repository that contains this site
(``--src PATH`` points to another checkout); the working tree is not read.

Usage:
    uv run python tools/sync_from_code_repo.py            # apply
    uv run python tools/sync_from_code_repo.py --check    # report only, exit 1 on drift

The site lives in the ``website/`` folder of the code repository and uses the
repository's files directly (docs/development, docs/user-guide, examples/,
lanczosKuboCode/docs/): nothing is mirrored. Generated here:
    lanczosKuboCode/Src                         -> reference/units.md, reference/api/ (gen_api_reference.py)

Hand-written pages that are derived from upstream sources cannot be synced
automatically. Their upstream git object hashes are recorded in
tools/sync-state.json and any change is reported as MANUAL:
    docs/development/README.md        -> development/index.md
    docs/theory/tapw-procedure.tex    -> theory/tapw-procedure.md
    examples/README.md, examples/*/README.md -> examples/*.md
    tools/input/{gendata,parser,gendata_schema}.py -> user-guide/yaml-json-input.md
    README.md, lanczosKuboCode/README.md -> installation.md, status.md, citing.md
    CITATION.cff                      -> citing.md
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from git import Repo

ROOT = Path(__file__).resolve().parent.parent
STATE = ROOT / "tools" / "sync-state.json"
SRC_UNITS = "lanczosKuboCode/Src"
WATCH = [
    "README.md",
    "CITATION.cff",
    "lanczosKuboCode/README.md",
    "docs/development/README.md",
    "docs/theory/tapw-procedure.tex",
    "examples/README.md",
    "tools/input/gendata.py",
    "tools/input/parser.py",
    "tools/input/gendata_schema.py",
]

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=str(ROOT.parent), help="checkout of the code repo (default: the one containing this site)")
    ap.add_argument("--check", action="store_true", help="report drift without writing; exit 1 if any")
    a = ap.parse_args()

    repo = Repo(a.src)
    commit = repo.head.commit.hexsha
    print(f"source: {repo.working_dir} @ {commit[:7]}")
    tree = repo.head.commit.tree

    # hand-written pages and the generated reference: compare recorded object hashes
    state = json.loads(STATE.read_text()) if STATE.exists() else {"trees": {}}
    watch = WATCH + sorted(
        b.path for b in tree["examples"].traverse() if b.type == "blob" and b.path.count("/") == 2
        and b.path.endswith("/README.md")
    )
    hashes = {p: tree[p].hexsha for p in watch}
    hashes[SRC_UNITS] = tree[SRC_UNITS].hexsha
    manual = [p for p in watch if state["trees"].get(p) != hashes[p]]
    units_stale = state["trees"].get(SRC_UNITS) != hashes[SRC_UNITS] or not (ROOT / "reference/units.md").exists()

    if units_stale:
        print("  regenerate  reference/units.md")
    for p in manual:
        print(f"  MANUAL      upstream {p} changed (or no baseline): review the page derived from it")

    drift = bool(units_stale or manual)
    if not drift:
        print("up to date")
        return 0
    if a.check:
        print("drift found (check mode, nothing written)")
        return 1

    if units_stale:
        # gen_api_reference.py reads sources from disk, so give it the checkout
        subprocess.run(
            [sys.executable, str(ROOT / "tools/gen_api_reference.py"), "--src", str(Path(repo.working_dir) / SRC_UNITS)],
            check=True, cwd=ROOT,
        )
    STATE.write_text(json.dumps({"commit": commit, "trees": hashes}, indent=2, sort_keys=True) + "\n")
    print(f"synced to {commit[:7]}; review MANUAL items and run the -W build")
    return 0


if __name__ == "__main__":
    sys.exit(main())
