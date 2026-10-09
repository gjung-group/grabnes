#!/usr/bin/env python3
"""Survey of the model switches of GRABNES: does each one run, and what Hamiltonian does it build?

Every case is a small input (tens to hundreds of atoms) derived from a base input by a few
key changes. The solver only assembles the Hamiltonian (a ten-step recursion follows, so that
the run ends normally) and writes its tables (WriteDataFiles). For every case the script reports

  status     ok        the run ended normally and the tables are consistent
             INERT     as ok, but the Hamiltonian is identical to that of the base input:
                       the switch had no effect in this configuration
             NONHERM   the two directions of some bond differ, or a reverse entry is missing
             NAN       the tables contain NaN or infinite values
             REFUSED   the solver stopped with its own error message (shown)
             CRASH     the run ended abnormally (runtime check, signal, time limit; shown)
  atoms, entries, largest asymmetry, and a fingerprint of the Hamiltonian

and, with --check, compares status and fingerprint with the stored reference
(model_survey_reference.json), so that a change of any Hamiltonian is detected.

    python3 model_survey.py --bin /path/to/grabnes [--work DIR] [--only PATTERN]
                            [--check | --update-reference] [--report FILE.md]

The fingerprint is a hash of the sorted list of (atom, neighbour, lattice translation, hopping)
and of the on-site energies, rounded to 1e-9 g0. It does not depend on the order of the
neighbour list. Use a checked build (tests/regression/config/gfortran.debug.make.sys) to turn
out-of-bounds accesses into CRASH lines instead of silent wrong numbers.

Requires NumPy and tools/hamiltonian/verify_tables.py (for the table reader).
"""

import argparse
import fnmatch
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools", "hamiltonian"))
from verify_tables import load  # noqa: E402

sys.path.insert(0, HERE)
from model_cases import BASES, CASES, STRUCTURES  # noqa: E402

REFERENCE = os.path.join(HERE, "model_survey_reference.json")

# Appended to every input: assemble the Hamiltonian, write the tables, stop after a short recursion.
RUN_MODE = {
    "Prefix": "generate", "Kubo.Calc": ".true.", "Diag.Calc": ".false.", "Calculate.OnlyDOS": ".true.",
    "Calculate.Bands": ".false.", "RecursionNumber": "10", "NumberofEnergyPoints": "11", "Epsilon": "0.05",
    "EnergyMin": "0-9.0", "EnergyMax": "9.0", "setSeed": ".true.", "seedValue": "123456",
    "WriteDataFiles": ".true.", "WritePos": ".false.",
}


def build_input(base, changes):
    """Base text with the keys of `changes` replaced (first occurrence) or appended; value None removes a key."""
    lines = BASES[base].strip("\n").split("\n")
    todo = dict(RUN_MODE)
    todo.update(changes)
    out, inblock = [], False
    for line in lines:
        s = line.strip()
        if s.startswith("&begin"):
            inblock = True
        if s.startswith("&end"):
            inblock = False
            out.append(line)
            continue
        key = s.split()[0] if s and not s.startswith("#") and not inblock else None
        hit = next((k for k in todo if key is not None and k.lower() == key.lower()), None)
        if hit is None:
            out.append(line)
        else:
            value = todo.pop(hit)
            if value is not None:
                out.append(f"{hit} {value}")
    # a key starting with "&" carries a whole input block as its value
    out += [v if k.startswith("&") else f"{k} {v}" for k, v in todo.items() if v is not None]
    return "\n".join(out) + "\n"


def classify(run_dir, rc, timed_out):
    log = open(os.path.join(run_dir, "job.out"), errors="replace").read() if os.path.exists(
        os.path.join(run_dir, "job.out")) else ""
    tail = [l.strip() for l in log.split("\n") if l.strip()]
    finished = "Program finished" in log
    info = dict(status="ok", message="", atoms=0, entries=0, asym=0.0, unpaired=0, fingerprint="")
    if timed_out:
        info.update(status="CRASH", message="time limit")
        return info
    if not finished:
        runtime = [l for l in tail if re.search(r"Fortran runtime error|SIGSEGV|SIGFPE|SIGABRT|severe|"
                                                r"Segmentation|Backtrace|Index '", l)]
        killed = [l for l in tail if re.search(r"\bERROR\b|Error:|killed|MIO_Kill|not supported|must be", l)]
        if runtime:
            info.update(status="CRASH", message=" | ".join(runtime[:2])[:300])
        elif killed:
            start = next(i for i, l in enumerate(tail) if l == killed[0])
            info.update(status="REFUSED", message=" ".join(tail[start:start + 3])[:300])
        else:
            info.update(status="CRASH", message=("rc=%d; last line: " % rc + (tail[-1] if tail else "no output"))[:300])
        return info
    try:
        t = load(run_dir, "generate")
    except Exception as exc:  # the run ended but the tables are unusable
        smag = os.path.join(run_dir, "generate.s.mag")
        if os.path.exists(smag) and re.search(r"NaN|Infinity", open(smag).read()):
            info.update(status="NAN", message="non-finite values in the hopping table")
        else:
            info.update(status="CRASH", message=f"tables not readable: {exc}"[:300])
        return info
    n = len(t["h0"])
    info.update(atoms=n, entries=len(t["hop"]))
    if not (np.all(np.isfinite(t["hop"].real)) and np.all(np.isfinite(t["hop"].imag)) and np.all(np.isfinite(t["h0"]))):
        info.update(status="NAN", message="non-finite values in the hopping or on-site table")
        return info
    # reverse entry of (i -> m, R) is (m -> i, -R) with the conjugate value
    key = {}
    for e in range(len(t["hop"])):
        key.setdefault((int(t["owner"][e]), int(t["neighbor"][e]), tuple(int(x) for x in t["image"][e])), []).append(e)
    asym, unpaired, duplicates = 0.0, 0, sum(len(v) - 1 for v in key.values())
    for (i, m, r), es in key.items():
        partner = key.get((m, i, tuple(-x for x in r)))
        if partner is None:
            unpaired += 1
            continue
        asym = max(asym, abs(t["hop"][es[0]] - np.conj(t["hop"][partner[0]])))
    info.update(asym=float(asym), unpaired=unpaired)
    rows = sorted((i, m) + r + (round(float(t["hop"][es[0]].real), 9) + 0.0, round(float(t["hop"][es[0]].imag), 9) + 0.0)
                  for (i, m, r), es in key.items())
    h = hashlib.sha1()
    h.update(repr(rows).encode())
    h.update(repr([round(float(x), 9) + 0.0 for x in t["h0"]]).encode())
    info["fingerprint"] = h.hexdigest()[:16]
    notes = []
    if duplicates:
        notes.append(f"{duplicates} duplicate entries")
    if unpaired:
        notes.append(f"{unpaired} entries without reverse entry")
    if asym > 1e-9:
        notes.append(f"largest asymmetry {asym:.3e} g0")
    if unpaired or asym > 1e-9:
        info["status"] = "NONHERM"
    warn = [l for l in tail if "WARNING" in l]
    if warn:
        notes.append("log: " + warn[0][:120])
    info["message"] = "; ".join(notes)
    return info


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--bin", required=True)
    p.add_argument("--work", default=os.path.join(HERE, "work", "model_survey"))
    p.add_argument("--only", default="*", help="shell pattern on the case name")
    p.add_argument("--launcher", default="mpirun -np 1")
    p.add_argument("--timeout", type=int, default=120)
    p.add_argument("--check", action="store_true", help="compare with the stored reference; exit 1 on a difference")
    p.add_argument("--update-reference", action="store_true")
    p.add_argument("--report", help="write the table as Markdown")
    a = p.parse_args()

    env = dict(os.environ, OMP_NUM_THREADS="2")
    results, base_print = {}, {}
    names = [c for c in CASES if fnmatch.fnmatch(c, a.only) or c in BASES]
    for name in names:
        base, changes, *rest = CASES[name]
        files = rest[0] if rest else None
        d = os.path.join(a.work, name)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        open(os.path.join(d, "Gendata.in"), "w").write(build_input(base, changes))
        if files:
            STRUCTURES[files](d)
        timed_out, rc = False, 0
        with open(os.path.join(d, "job.out"), "w") as out:
            try:
                rc = subprocess.run(a.launcher.split() + [os.path.abspath(a.bin), "Gendata.in"], cwd=d, stdout=out,
                                    stderr=subprocess.STDOUT, env=env, timeout=a.timeout).returncode
            except subprocess.TimeoutExpired:
                timed_out = True
        r = classify(d, rc, timed_out)
        r["base"] = base
        if name == base:
            base_print[base] = r["fingerprint"]
        elif r["status"] == "ok" and r["fingerprint"] and r["fingerprint"] == base_print.get(base):
            r["status"] = "INERT"
        results[name] = r
        print(f"{name:34s} {r['status']:8s} atoms {r['atoms']:5d} entries {r['entries']:7d} {r['fingerprint']:16s} "
              f"{r['message']}", flush=True)

    count = {}
    for r in results.values():
        count[r["status"]] = count.get(r["status"], 0) + 1
    print("\nsummary: " + ", ".join(f"{k} {v}" for k, v in sorted(count.items())))

    if a.report:
        with open(a.report, "w") as f:
            f.write("| case | base | status | atoms | entries | fingerprint | remarks |\n|---|---|---|---|---|---|---|\n")
            for name, r in results.items():
                f.write(f"| `{name}` | `{r['base']}` | {r['status']} | {r['atoms']} | {r['entries']} | "
                        f"`{r['fingerprint']}` | {r['message'].replace('|', '/')} |\n")
    compact = {k: dict(status=v["status"], fingerprint=v["fingerprint"]) for k, v in results.items()}
    if a.update_reference:
        ref = json.load(open(REFERENCE)) if os.path.exists(REFERENCE) else {}
        ref.update(compact)
        json.dump(ref, open(REFERENCE, "w"), indent=1, sort_keys=True)
        print(f"reference updated: {REFERENCE}")
    if a.check:
        ref = json.load(open(REFERENCE))
        bad = [k for k, v in compact.items() if ref.get(k) != v]
        for k in bad:
            print(f"DIFFERS from the reference: {k}: {ref.get(k)} -> {compact[k]}")
        print("PASS: every case equals its reference" if not bad else f"FAIL: {len(bad)} case(s) differ")
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
