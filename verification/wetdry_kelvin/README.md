# wetdry_kelvin: a coastal Kelvin wave over drying flats and into a fjord

A Kelvin-wave pulse travels north along a steep coast, turns into a fjord,
runs along the fjord walls and comes back out. On its way it floods and
drains intertidal flats along the coast and inside the fjord. It tests
wetting and drying in a short 3-D run with rotation: faces open and close
at every stage, and volume must stay conserved to round-off.

## Set-up

- 100 x 100 x 15 grid, dx = dy = 1 km (4 tiles of 50 x 50). The 15 levels
  are 1-8 m thick (54 m in total), r* (`select_rStar=2`).
- Coast on the east side. From west to east: 50 m deep ocean up to the
  coast, a steep slope (50 m to 2 m over 2 km), 3 km of intertidal flats
  (bed +0.2, +0.6, +1.0 m above sea level), then land.
- Fjord at 15 < y < 26 km, 11 km wide, cutting 10 km into the land (head at
  x = 85 km): a 30 m deep channel with steep walls and flats (bed +0.1 and
  +0.3 m) along both walls and at its head.
- The west edge is a land column. The north and south edges are open ocean:
  the domain is periodic in y, so a wave that leaves the north edge comes
  back in at the south edge (little of it within 60 steps).
- f = 3e-3 1/s (about 20 times the Earth value), so the Rossby radius
  sqrt(g H)/f = 7.4 km is smaller than the fjord width and the wave turns
  into the fjord along its south wall.
- Initial state: a Kelvin-wave pulse (amplitude 2.5 m, centered at
  y = 10 km, e-folding length 4 km along the coast, decaying offshore over
  the Rossby radius) with the northward velocity in geostrophic balance.
  All flats start dry (film of 5 cm). Salinity is uniform (35).
- The model rest level is 4 m above sea level (`DATUM` in `gendata.py`).
- Time step 90 s; 60 steps (90 min). The wave moves about 2 km per step.
- Implicit bottom drag (`selectImplicitDrag=2`, `#define
  ALLOW_SOLVE4_PS_AND_DRAG`). Without it, the run stops at step 6 with
  34 m/s in a thin film at the wetting front.
- Packages: `gfd diagnostics wetdry`.

## What happens

| Steps | Event |
|---|---|
| 1-10 | the pulse moves north and floods the coast flats (up to 13 wet cells) |
| 5-26 | the wave turns into the fjord, runs east along the south wall, across the head and back along the north wall; up to 50 of the 66 fjord flat cells are wet |
| 10-18 | the coast flats drain |
| 26-60 | the fjord flats drain (53 of the 54 cells that were wet are dry again by step 60) |

The `%WETDRY_MON` volume budget error stays below 2e-15 of the total volume.

## Inputs

All binary inputs (`bathy.bin`, `eta.bin`, `V.bin`, 32-bit) are in the
repository. To regenerate them, run `python3 gendata.py` in `input/`.
`data.diagnostics` writes snapshots of `ETAN`, `WDdepth`, `WDdryC`,
`WDmaskW` and the surface velocity every 10 steps. A file labeled n holds
the state after step n+1. Set the frequencies to `-90.` for every step.

## Build and run

With testreport, from MITgcm's `verification/` directory:

```sh
./testreport -t wetdry_kelvin            # add -of <optfile>, -mpi as usual
```

By hand:

```sh
cd verification/wetdry_kelvin
mkdir -p build run
cd build
../../../tools/genmake2 -mods ../code -of <optfile>
make depend && make
cd ../run
ln -s ../input/* .
../build/mitgcmuv > output.txt
```

## Checks

- `PROGRAM MAIN: Execution ended Normally` after 60 steps.
- `%WETDRY_MON: ... budgetErr/vol0` below 2e-15 at every monitor output.
- `nFlips` > 0 in most monitor intervals (faces open and close), and `nDry`
  falls from 333 to about 265 and rises again.
- Salinity stays 35 everywhere (`%MON dynstat_salt_min/max`).
