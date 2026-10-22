#!/usr/bin/env python3
"""Inputs for wetdry_kelvin: a coastal Kelvin wave that floods and drains
intertidal flats along a coast and inside a fjord, on an f-plane at 60 N.

Grid 100 x 100 x 15, dx = dy = 24 km (2400 x 2400 km). The geometry below
is set in grid cells. North-south coast on the east side (f > 0, so a
Kelvin wave travels north with the coast on its right). From west to east:
50 m deep ocean up to the coast, a steep slope (50 m to 2 m over 2 cells),
3 cells of intertidal flats (bed +0.2, +0.6, +1.0 m), then land. A fjord
(a wide gulf at this grid spacing), 11 cells wide, cuts 10 cells east into
the land at cells 15-26 in y (head at cell 85 in x): a 30 m deep channel
with steep walls and intertidal flats (bed +0.1 and +0.3 m) along both
walls and at its head. The west edge is a land column. The north and south
edges are open ocean: the domain is periodic in y, so a wave that leaves
the north edge comes back in at the south edge (little of it within 60
steps).

f = 2 Omega sin(60 deg) = 1.263e-4 1/s (f-plane, no beta). The Rossby
radius sqrt(g H)/f = 175 km (7.3 cells) is smaller than the fjord width,
so the wave turns into the fjord along its south wall. c = sqrt(g H) =
22 m/s, so with deltaT = 2160 s it moves about 2 cells per step. In 60
steps (36 h) the coast flats flood and drain (steps 1-38) and the fjord
flats flood (from step 5, up to 51 cells at step 29) and drain again (48
of 52 dry by step 60). The bottom drag (in data) is low, 1e-4, so that
the thin water on the 24 km wide flats drains within 60 steps.

Initial state: a Kelvin-wave pulse of amplitude A = 2.5 m near the
southern edge (centred at cell 10 in y, e-folding length 4 cells along the
coast, decaying offshore as exp(-(x_foot - x)/R) from the foot of the
slope), with the northward velocity in geostrophic balance,
v = (g/f) d(eta)/dx. The pulse is put only on columns whose bed is below
sea level, so all flats start dry (film of 5 cm). Salinity is uniform
(sRef = 35, set in data).

Model rest level DATUM = 4 m above sea level (r_low = bed - DATUM): the
model surface starts at eta = -DATUM.

Usage:  python3 gendata.py [OUTDIR]     (default: this directory)
"""
import os
import sys

import numpy as np

NX, NY, NR = 100, 100, 15
DX = 24.0e3                     # grid spacing (m); geometry below in cells
G = 9.81
LAT = 60.0                      # f-plane latitude (deg N)
F0 = 2 * 7.292115e-5 * np.sin(np.deg2rad(LAT))   # 1.263e-4 1/s
H0 = 50.0                       # offshore depth (m)
DATUM = 4.0                     # model rest level above sea level (m)
FILM = 0.05                     # = wdMinDepth
A, Y0, LY = 2.5, 10 * DX, 4 * DX  # pulse amplitude (m), centre, e-folding (m)
DELR = np.array([1, 1, 1, 1.5, 1.5, 2, 2, 3, 3, 4, 5, 6, 7, 8, 8.0])
assert DELR.size == NR and abs(DELR.sum() - (H0 + DATUM)) < 1e-9

x = (np.arange(NX) + 0.5) * DX
y = (np.arange(NY) + 0.5) * DX
X, Y = np.meshgrid(x, y)        # shape (NY, NX): [j, i]

# --- coast: deep ocean, slope, flats, land (bed z, positive up, m) ---
X_FOOT, X_TOP, X_COAST = 70 * DX, 72 * DX, 75 * DX
zb = np.full((NY, NX), -H0)
slope = (X >= X_FOOT) & (X < X_TOP)
zb[slope] = np.where(X[slope] < X_FOOT + DX, -25.0, -2.0)
flat = (X >= X_TOP) & (X < X_COAST)
zb[flat] = 0.2 + 0.4 * np.floor((X[flat] - X_TOP) / DX)
ocean = X < X_COAST

# --- fjord: 30 m channel, flats along the walls and at the head ---
FJ_S, FJ_N, FJ_HEAD = 15 * DX, 26 * DX, 85 * DX
fj = (Y > FJ_S) & (Y < FJ_N) & (X >= X_FOOT) & (X < FJ_HEAD)
dwall = np.minimum(np.minimum(Y - FJ_S, FJ_N - Y), FJ_HEAD - X)  # to a wall
zfj = np.where(dwall < 1 * DX, 0.3, np.where(dwall < 2 * DX, 0.1, -30.0))
zb[fj] = np.where(X[fj] < X_COAST, np.minimum(zb[fj], zfj[fj]), zfj[fj])
ocean |= fj

# --- west edge: one land column; the north and south edges stay ocean
# (periodic in y) ---
ocean[:, 0] = False
r_low = np.where(ocean, zb - DATUM, 0.0)

# --- Kelvin pulse on wet columns only (bed below sea level) ---
wet0 = ocean & (zb < -FILM)
R = np.sqrt(G * H0) / F0
xprof = np.exp(-np.clip(X_FOOT - X, 0, None) / R)
eta1 = A * xprof * np.exp(-((Y - Y0) / LY) ** 2)
eta1 = np.where(wet0 & (X < X_TOP) & ~fj, eta1, 0.0)
eta = np.where(ocean, np.maximum(-DATUM + eta1, r_low + FILM), 0.0)
# v at the southern face of cell j: v = (g/f) d(eta)/dx, in deep water
# only (zero on the slope and flats, where the pulse is flat across shore)
detadx = np.zeros_like(eta1)
detadx[:, 1:-1] = (eta1[:, 2:] - eta1[:, :-2]) / (2 * DX)
vc = np.where(X < X_FOOT, G / F0 * detadx, 0.0)
v = np.zeros_like(vc)
v[1:, :] = 0.5 * (vc[1:, :] + vc[:-1, :])
v[1:, :] *= (ocean[1:, :] & ocean[:-1, :])
v3 = np.repeat(v[None], NR, 0)          # depth-uniform (barotropic)

if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.dirname(os.path.abspath(__file__))

    def w(name, a):
        a.astype('>f4').tofile(os.path.join(out, name))   # readBinaryPrec=32

    w('bathy.bin', r_low)
    w('eta.bin', eta)
    w('V.bin', v3)
    D = eta - r_low
    print(f'R = {R/1e3:.1f} km, c = {np.sqrt(G*H0):.1f} m/s;',
          'ocean cells', int(ocean.sum()), ', dry at start',
          int((ocean & (D < FILM + 0.01)).sum()),
          ', intertidal (bed above sea level)',
          int((ocean & (zb > -1.5) & (zb < 2.5)).sum()),
          f', max v {v.max():.2f} m/s')
