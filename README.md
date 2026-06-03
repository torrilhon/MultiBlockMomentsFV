# MultiBlockMomentsFV

A classic finite-volume solver for the 2D **13-moment (Grad/Regularized-13)** equations
of rarefied gas dynamics on Cartesian multi-block meshes.

The scheme is a second-order MUSCL reconstruction with an HLL-type numerical flux,
SSP-RK2 (Heun) time stepping, and a direct (relaxation) integration of the stiff
production terms. The numerical method is described in the accompanying paper
(`literature/Torrilhon_MMS.pdf`).

## Repository layout

```
MultiBlockMomentsFV/
├── Makefile          build / run / check / clean
├── src/              solver sources (.c/.h)
├── cases/            input decks ("block files")
│   ├── square4.txt        4-block square domain (standard test, t_end=0.3)
│   ├── squareSmall.txt    coarse version of square4 (quick smoke test)
│   └── Lshape.txt         3-block L-shaped domain (t_end=0.6)
├── reference/        golden results for validation (see note below)
├── tools/            post-processing (plot_solution.jl)
├── literature/       Torrilhon_MMS.pdf (method description)
├── build/            compiled objects + binary  (git-ignored)
└── runs/             simulation output           (git-ignored)
```

## Build

Requires a C compiler with the math library (`gcc` recommended).

```sh
make            # compiles src/*.c -> build/moments
make clean      # removes build/
```

## Run

Each case writes into its own `runs/<case>/` directory:

```sh
make run CASE=cases/squareSmall.txt          # quick test (~100 steps)
make run CASE=cases/square4.txt              # standard reference run
make run CASE=cases/Lshape.txt               # L-shaped domain
```

Pass extra solver flags via `ARGS`:

```sh
make run CASE=cases/square4.txt ARGS="-Kn 0.01 -CFL 0.9"
```

Or run the binary directly:

```sh
./build/moments -bfile cases/squareSmall.txt -out runs/squareSmall
```

### Command-line flags

| Flag           | Meaning                                              | Default |
|----------------|------------------------------------------------------|---------|
| `-bfile <file>`| block/input file defining grid + output times        | `blockfile.txt` |
| `-out <dir>`   | directory for output files (created if missing)      | current dir |
| `-Kn <value>`  | Knudsen number                                       | `0.001` |
| `-R13 <value>` | regularization switch: `<1` = Grad13, `>=1` = Reg13  | `0.0`   |
| `-CFL <value>` | target CFL number for adaptive time stepping         | `0.95`  |
| `-sfile <file>`| restart: load initial state from a previous output   | (none)  |

## Block file format

A block file has a header for the output times followed by one line per block:

```
Output times, maximal 20 times, increasing or dt (negative) and tend:
-0.05  0.3
Block data (Maximal NB blocks!):
number/xmin,xmax/ymin,ymax/nstx,nsty/leftbound,rightbound,bottombound,topbound/
1/-1.0,0.5/-1.0,0.2/15,12/0,2,-1,4/
...
```

- **Output times line:** either a list of increasing times, or `-dt  tend`
  (a negative first value means "every `dt` until `tend`").
- **Per-block line:** block number, x-range, y-range, cell counts `nstx,nsty`,
  and the four neighbor codes for left/right/bottom/top edges.
- **Boundary codes:** `0` = extrapolation (open), `-1` = solid wall (adiabatic),
  `>0` = synchronize with that block number (internal multi-block interface).

The maximum number of blocks is `NB` in `src/defs.h` (currently 10); a block
number outside `[1, NB]` is rejected with an error.

## Output

Results are written as **Tecplot ASCII** `.dat` files named
`Grad13<Type>_nst<gridsize>_Kn<value>_tAll.dat`, containing all output time levels
as separate `ZONE`s. Variables per point:

```
x, y, rho, vx, vy, T, p, Pxx, Pxy, Pyy, qx, qy
```

View with Tecplot/ParaView, or plot in Julia:

```sh
julia tools/plot_solution.jl runs/squareSmall/Grad13Riemann2D_*_tAll.dat
```

## Validation

```sh
make check CASE=cases/squareSmall.txt
```

This runs the case and diffs its output against `reference/<case>.ref.dat`.

> **Note:** the reference data currently stored in `reference/` is known to be
> corrupted and needs to be regenerated. Until then `make check` reports the
> diff but does not fail the build. Once a correct golden file is in place,
> remove the `|| true` fallback in the Makefile to make `check` a hard gate.

---
(c) Manuel Torrilhon, ETH Zurich 2010 – RWTH Aachen 2026
