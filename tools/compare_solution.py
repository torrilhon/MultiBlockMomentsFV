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
# The .dat payload precision is read from the .hdr [BINARY] PRECISION field
# (Float32 default, Float64 with -double). Each file is decoded at its own
# precision and widened to Python float for the comparison, so a Float32 run can
# be checked against a Float64 reference (or vice versa).
#
# Default tolerances adapt to the LOOSER of the two precisions: a Float32
# payload only carries ~7 significant digits, so comparing it bit-for-bit
# against a double reference is meaningless. Defaults:
#   both Float64           -> abstol=1e-10, reltol=1e-8
#   either side Float32    -> abstol=1e-6,  reltol=1e-5

import sys
import os
import re
import struct


def companion(path):
    """Return (hdr, dat) paths from a .hdr, .dat, or bare basename."""
    base = path[:-4] if path.endswith((".hdr", ".dat")) else path
    return base + ".hdr", base + ".dat"


def parse_precision(hdr):
    """Return the payload precision ('Float32' or 'Float64') from the [BINARY]
    PRECISION field. Defaults to Float64 for older headers without the field."""
    for line in open(hdr):
        m = re.match(r"\s*PRECISION\s*=\s*(\S+)", line)
        if m:
            return m.group(1)
    return "Float64"


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
    """Return (zone_keys, values, precision) where zone_keys is a structure
    fingerprint list (block,I,J,NVARS per zone, STRANDID-free) and values is a
    flat list of floats (widened from the file precision) over all zones in file
    order."""
    hdr, dat = companion(path)
    prec = parse_precision(hdr)
    code, esize = ("f", 4) if prec == "Float32" else ("d", 8)
    zones = parse_zones(hdr)
    keys = [(z["block"], z["I"], z["J"], z["nvars"]) for z in zones]
    values = []
    with open(dat, "rb") as f:
        for z in zones:
            f.seek(z["offset"])
            n = z["nvars"] * z["I"] * z["J"]
            values.extend(struct.unpack(f"<{n}{code}", f.read(n * esize)))
    return keys, values, prec


def main():
    if len(sys.argv) < 3:
        print("usage: compare_solution.py <new> <reference> [abstol] [reltol]")
        return 2
    new, ref = sys.argv[1], sys.argv[2]

    kn, dn, pn = load(new)
    kr, dr, pr = load(ref)

    # Default tolerances follow the looser precision present on either side.
    single = "Float32" in (pn, pr)
    abstol = float(sys.argv[3]) if len(sys.argv) > 3 else (1e-6 if single else 1e-10)
    reltol = float(sys.argv[4]) if len(sys.argv) > 4 else (1e-5 if single else 1e-8)

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
          f"(abstol={abstol:g}, reltol={reltol:g}, precision new={pn}/ref={pr}); "
          f"max |diff| = {max_abs:.3e}")
    if worst and max_abs > 0:
        i, x, y = worst
        print(f"      largest diff at index {i}: {x:.12g} vs {y:.12g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
