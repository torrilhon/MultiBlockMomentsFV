#!/usr/bin/env python3
# compare_solution.py -- numeric regression comparison of two Tecplot ASCII
# solution files written by the solver.
#
#   python3 tools/compare_solution.py <new.dat> <reference.dat> [abstol] [reltol]
#
# Compares the two files field-by-field with a tolerance, ignoring harmless
# differences (signed zeros, whitespace, Tecplot STRANDID metadata) that are
# not physically meaningful. Exit code 0 = match (PASS), 1 = mismatch (FAIL).
#
# Defaults: abstol=1e-10, reltol=1e-8  (suitable for a deterministic build).

import sys


def load(path):
    """Return (zone_headers, data_rows) where data_rows is a flat list of float tuples."""
    headers, rows = [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith(("TITLE", "VARIABLES")):
                continue
            if s.startswith("ZONE"):
                # keep zone identity but drop volatile STRANDID metadata
                headers.append(",".join(
                    tok for tok in s.split(",") if "STRANDID" not in tok))
                continue
            rows.append(tuple(float(x) for x in s.split()))
    return headers, rows


def main():
    if len(sys.argv) < 3:
        print("usage: compare_solution.py <new.dat> <reference.dat> [abstol] [reltol]")
        return 2
    new, ref = sys.argv[1], sys.argv[2]
    abstol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10
    reltol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-8

    hn, dn = load(new)
    hr, dr = load(ref)

    if hn != hr:
        print(f"FAIL: zone structure differs ({len(hn)} vs {len(hr)} zones)")
        for a, b in zip(hn, hr):
            if a != b:
                print(f"   new: {a}\n   ref: {b}")
                break
        return 1
    if len(dn) != len(dr):
        print(f"FAIL: row count differs ({len(dn)} vs {len(dr)})")
        return 1

    max_abs = 0.0
    worst = None
    for i, (a, b) in enumerate(zip(dn, dr)):
        if len(a) != len(b):
            print(f"FAIL: column count differs at data row {i}")
            return 1
        for j, (x, y) in enumerate(zip(a, b)):
            d = abs(x - y)
            if d > max_abs:
                max_abs = d
                worst = (i, j, x, y)
            if d > abstol + reltol * abs(y):
                print(f"FAIL: value mismatch at row {i}, col {j}: "
                      f"{x:.12g} vs {y:.12g} (|diff| = {d:.3e})")
                return 1

    print(f"PASS: {len(dn)} rows match within tol "
          f"(abstol={abstol:g}, reltol={reltol:g}); max |diff| = {max_abs:.3e}")
    if worst and max_abs > 0:
        i, j, x, y = worst
        print(f"      largest diff at row {i}, col {j}: {x:.12g} vs {y:.12g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
