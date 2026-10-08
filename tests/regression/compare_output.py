#!/usr/bin/env python3
"""Compare a GRABNES output file with its reference data.

Only the Python standard library is used so that the regression tests do not
depend on NumPy. The comparison is numerical: both files must contain the same
number of lines and the same number of values on every line, and every value
must agree within ``atol + rtol * |reference|``. Byte identity is reported for
information; it is only required with ``--exact``.

Exit status: 0 on agreement, 1 on disagreement, 2 on unreadable input.
"""

import argparse
import math
import sys


def read_table(path):
    rows = []
    try:
        with open(path, "r") as handle:
            for number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                try:
                    rows.append([float(token) for token in line.split()])
                except ValueError:
                    raise SystemExit2(f"{path}:{number}: not a numeric line: {line.rstrip()!r}")
    except OSError as exc:
        raise SystemExit2(f"cannot read {path}: {exc}")
    if not rows:
        raise SystemExit2(f"{path}: no data")
    return rows


class SystemExit2(Exception):
    """Unusable input (exit status 2)."""


def check_bands(rows, label):
    """Structural checks of a GRABNES ``*.bands`` file."""
    problems = []
    if len(rows) < 5 or len(rows[3]) != 3:
        return [f"{label}: missing the 4-line band-file header"]
    bands, spins, points = (int(value) for value in rows[3])
    width = 1 + bands * spins
    values = [value for row in rows[4:] for value in row]
    if len(values) != width * points:
        return [
            f"{label}: header announces {points} k-points x {bands} bands x {spins} spin(s) "
            f"but the file holds {len(values)} values"
        ]
    previous_x = -math.inf
    for point in range(points):
        record = values[point * width : (point + 1) * width]
        if record[0] < previous_x:
            problems.append(f"{label}: k-path coordinate decreases at point {point + 1}")
        previous_x = record[0]
        for spin in range(spins):
            energies = record[1 + spin * bands : 1 + (spin + 1) * bands]
            if any(b < a for a, b in zip(energies, energies[1:])):
                problems.append(f"{label}: eigenvalues not in ascending order at point {point + 1}")
                break
    return problems


def check_dos(rows, label):
    """Structural checks of a GRABNES ``*.diag.DOS`` file."""
    problems = []
    if any(len(row) != len(rows[0]) or len(row) < 2 for row in rows):
        return [f"{label}: rows do not all have the same (>= 2) number of columns"]
    energies = [row[0] for row in rows]
    if any(b <= a for a, b in zip(energies, energies[1:])):
        problems.append(f"{label}: energy grid is not strictly increasing")
    if any(value < 0.0 for row in rows for value in row[1:]):
        problems.append(f"{label}: negative density of states")
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("result")
    parser.add_argument("reference")
    parser.add_argument("--kind", choices=("bands", "dos", "table"), default="table")
    parser.add_argument("--atol", type=float, default=0.0)
    parser.add_argument("--rtol", type=float, default=0.0)
    parser.add_argument("--exact", action="store_true", help="also require byte identity")
    args = parser.parse_args()

    try:
        result = read_table(args.result)
        reference = read_table(args.reference)
    except SystemExit2 as exc:
        print(f"  ERROR: {exc}")
        return 2

    problems = []
    if len(result) != len(reference):
        problems.append(f"line count differs: {len(result)} vs reference {len(reference)}")
    else:
        for number, (row, ref_row) in enumerate(zip(result, reference), start=1):
            if len(row) != len(ref_row):
                problems.append(
                    f"data line {number}: {len(row)} values vs reference {len(ref_row)}"
                )
                break

    max_abs = max_rel = 0.0
    worst = None
    count = failures = 0
    if not problems:
        for number, (row, ref_row) in enumerate(zip(result, reference), start=1):
            for column, (value, ref) in enumerate(zip(row, ref_row), start=1):
                count += 1
                if not math.isfinite(value):
                    problems.append(f"data line {number}, column {column}: non-finite value")
                    continue
                diff = abs(value - ref)
                if diff > max_abs:
                    max_abs, worst = diff, (number, column, value, ref)
                if ref != 0.0:
                    max_rel = max(max_rel, diff / abs(ref))
                if diff > args.atol + args.rtol * abs(ref):
                    failures += 1
        if failures:
            number, column, value, ref = worst
            problems.append(
                f"{failures} of {count} values outside tolerance; largest deviation at data "
                f"line {number}, column {column}: {value!r} vs reference {ref!r}"
            )
        checker = {"bands": check_bands, "dos": check_dos}.get(args.kind)
        if checker:
            problems.extend(checker(result, "result"))
            problems.extend(checker(reference, "reference"))

    with open(args.result, "rb") as a, open(args.reference, "rb") as b:
        identical = a.read() == b.read()
    if args.exact and not identical:
        problems.append("files are not byte-identical (--exact)")

    print(
        f"  values: {count}   max |diff|: {max_abs:.3e}   max rel diff: {max_rel:.3e}   "
        f"byte-identical: {'yes' if identical else 'no'}   "
        f"(atol {args.atol:g}, rtol {args.rtol:g})"
    )
    for problem in problems:
        print(f"  MISMATCH: {problem}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
