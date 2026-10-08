#!/usr/bin/env python3
"""Check a stochastic Kubo (recursion) DOS against an exact reference.

A random-phase state gives the DOS with a statistical error that decreases as
1/sqrt(number of atoms) and depends on the compiler's random-number generator,
so the result cannot be compared digit by digit. This check requires

- the same energy grid as the reference;
- a non-negative DOS;
- a maximum and a root-mean-square deviation below the given thresholds.

Only the Python standard library is used. Exit status: 0 pass, 1 fail,
2 unreadable input.
"""

import argparse
import math
import sys


def read(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.strip():
                rows.append([float(token) for token in line.split()[:2]])
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("result")
    parser.add_argument("reference")
    parser.add_argument("--max-dev", type=float, required=True)
    parser.add_argument("--rms-dev", type=float, required=True)
    args = parser.parse_args()

    try:
        result, reference = read(args.result), read(args.reference)
    except (OSError, ValueError) as exc:
        print(f"  ERROR: {exc}")
        return 2

    problems = []
    if len(result) != len(reference):
        problems.append(f"{len(result)} rows vs reference {len(reference)}")
    else:
        if max(abs(a[0] - b[0]) for a, b in zip(result, reference)) > 1e-9:
            problems.append("energy grid differs from the reference")
        if any(not math.isfinite(a[1]) or a[1] < 0.0 for a in result):
            problems.append("negative or non-finite DOS")
        diffs = [a[1] - b[1] for a, b in zip(result, reference)]
        max_dev = max(abs(d) for d in diffs)
        rms_dev = math.sqrt(sum(d * d for d in diffs) / len(diffs))
        peak = max(b[1] for b in reference)
        print(
            f"  values: {len(diffs)}   max |diff|: {max_dev:.4f} (limit {args.max_dev:g})   "
            f"rms diff: {rms_dev:.4f} (limit {args.rms_dev:g})   reference peak: {peak:.3f}"
        )
        if max_dev > args.max_dev:
            problems.append("maximum deviation above the limit")
        if rms_dev > args.rms_dev:
            problems.append("rms deviation above the limit")
    for problem in problems:
        print(f"  MISMATCH: {problem}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
