#!/usr/bin/env python3
# compare_solution.py -- numeric regression comparison of two binary (.hdr/.dat)
# solution results written by the solver.
#
#   python3 tools/compare_solution.py <new> <reference> [abstol] [reltol]
#
# Each argument may be the .hdr, the .dat, or the shared basename; the companion
# file is located automatically. Compares the two results field-by-field with a
# tolerance, ignoring harmless differences (volatile STRANDID / provenance
# metadata) that are not physically meaningful.
# Exit code 0 = match (PASS), 1 = mismatch (FAIL).
#
# Defaults: abstol=1e-10, reltol=1e-8  (suitable for a deterministic build).

import sys
import os
import re
import struct


def companion(path):
    """Return (hdr, dat) paths from a .hdr, .dat, or bare basename."""
    base = path[:-4] if path.endswith((".hdr", ".dat")) else path
    return base + ".hdr", base + ".dat"


def parse_zones(hdr):
    """Parse the [ZONES] records, returning a list of zone dicts."""
    zones = []
    for line in open(hdr):
        s = line.strip()
        if not s.startswith("ZONE"):
            continue
        g = lambda k: int(re.search(rf"\b{k}\s*=\s*(-?\d+)", s).group(1))
        zones.append(dict(block=g("block"), I=g("I"), J=g("J"),
                          nvars=g("NVARS"), offset=g("OFFSET"),
                          nbytes=g("NBYTES")))
    return zones


def load(path):
    """Return (zone_keys, values) where zone_keys is a structure fingerprint
    list (block,I,J,NVARS per zone, STRANDID-free) and values is a flat list of
    floats over all zones in file order."""
    hdr, dat = companion(path)
    zones = parse_zones(hdr)
    keys = [(z["block"], z["I"], z["J"], z["nvars"]) for z in zones]
    values = []
    with open(dat, "rb") as f:
        for z in zones:
            f.seek(z["offset"])
            n = z["nvars"] * z["I"] * z["J"]
            values.extend(struct.unpack(f"<{n}d", f.read(z["nbytes"])))
    return keys, values


def main():
    if len(sys.argv) < 3:
        print("usage: compare_solution.py <new> <reference> [abstol] [reltol]")
        return 2
    new, ref = sys.argv[1], sys.argv[2]
    abstol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10
    reltol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-8

    kn, dn = load(new)
    kr, dr = load(ref)

    if kn != kr:
        print(f"FAIL: zone structure differs ({len(kn)} vs {len(kr)} zones)")
        for a, b in zip(kn, kr):
            if a != b:
                print(f"   new: {a}\n   ref: {b}")
                break
        return 1
    if len(dn) != len(dr):
        print(f"FAIL: value count differs ({len(dn)} vs {len(dr)})")
        return 1

    max_abs = 0.0
    worst = None
    for i, (x, y) in enumerate(zip(dn, dr)):
        d = abs(x - y)
        if d > max_abs:
            max_abs = d
            worst = (i, x, y)
        if d > abstol + reltol * abs(y):
            print(f"FAIL: value mismatch at index {i}: "
                  f"{x:.12g} vs {y:.12g} (|diff| = {d:.3e})")
            return 1

    print(f"PASS: {len(dn)} values match within tol "
          f"(abstol={abstol:g}, reltol={reltol:g}); max |diff| = {max_abs:.3e}")
    if worst and max_abs > 0:
        i, x, y = worst
        print(f"      largest diff at index {i}: {x:.12g} vs {y:.12g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
