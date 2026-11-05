#!/usr/bin/env python3
"""Plot an external-flow (staircase body) result: velocity, pressure, Cd/Cl history.

Usage: python3 scripts/plot_airfoil.py results_external/naca0012_Re500_aoa5 [more case dirs]
Writes <case>/flow.png and <case>/forces.png.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm

BODY = ListedColormap(['#5f5e5a'])  # neutral gray for solid cells


def load_fields(case: Path):
    a = np.genfromtxt(case/'fields.csv', delimiter=',', names=True)
    x, y = np.unique(a['x']), np.unique(a['y'])
    shape = (len(y), len(x))  # file is written row by row (j outer, i inner)
    grid = {k: a[k].reshape(shape) for k in ('u', 'v', 'p', 'solid')}
    return x, y, grid


def title(case: Path) -> str:
    s = json.loads((case/'summary.json').read_text())
    name = f"NACA {s['naca']}, AoA {s['aoa_deg']:g}°" if s['body'] == 'naca' else s['body']
    cl = round(s['cl'], 3) + 0.0  # avoid printing -0.000
    return f"{name}, Re {s['Re']:g}  (Cd = {s['cd']:.3f}, Cl = {cl:.3f})"


def plot_flow(case: Path) -> None:
    L = 'D' if json.loads((case/'summary.json').read_text())['body'] == 'cylinder' else 'c'
    x, y, f = load_fields(case)
    solid = f['solid'] == 1
    speed = np.ma.masked_where(solid, np.hypot(f['u'], f['v']))
    # The solver pins p=0 in the inlet corner cell (free stream), so Cp = p/(0.5*rho*U^2) = 2p.
    cp = np.ma.masked_where(solid, 2*f['p'])
    body = np.ma.masked_where(~solid, np.ones_like(f['u']))
    # Zoom on the body: one length upstream, three downstream.
    js, is_ = np.nonzero(solid)
    xlim = (x[is_.min()] - 1.0, x[is_.max()] + 3.0)
    ylim = (y[js.min()] - 1.0, y[js.max()] + 1.0)
    fig, axes = plt.subplots(2, 1, figsize=(9, 7.5), constrained_layout=True)
    m = axes[0].pcolormesh(x, y, speed, cmap='Blues', shading='nearest')
    fig.colorbar(m, ax=axes[0], label='|U| / U∞')
    axes[0].streamplot(x, y, np.where(solid, 0, f['u']), np.where(solid, 0, f['v']),
                       color='#2c2c2a', linewidth=0.6, density=1.6, arrowsize=0.6)
    axes[0].set_title('Velocity magnitude and streamlines')
    lim = np.abs(cp).max()
    m = axes[1].pcolormesh(x, y, cp, cmap='RdBu_r', norm=TwoSlopeNorm(0, -lim, lim), shading='nearest')
    fig.colorbar(m, ax=axes[1], label='Cp')
    axes[1].set_title('Pressure coefficient')
    for ax in axes:
        ax.pcolormesh(x, y, body, cmap=BODY, shading='nearest')
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_aspect('equal')
        ax.set_xlabel(f'x / {L}'); ax.set_ylabel(f'y / {L}')
    fig.suptitle(title(case))
    fig.savefig(case/'flow.png', dpi=140)
    plt.close(fig)


def plot_forces(case: Path) -> None:
    h = np.genfromtxt(case/'forces.csv', delimiter=',', names=True)
    h = h[min(10, len(h)//10):]  # skip the start-up transient, which would dominate the y-range
    fig, axes = plt.subplots(2, 1, figsize=(7, 5), sharex=True, constrained_layout=True)
    for ax, key, label in ((axes[0], 'cd', 'Cd'), (axes[1], 'cl', 'Cl')):
        ax.plot(h['step'], h[key], color='#185fa5', lw=2)
        ax.set_ylabel(label)
        ax.grid(color='0.9'); ax.spines[['top', 'right']].set_visible(False)
    axes[1].set_xlabel('SIMPLE iteration' if h['time'][-1] == 0 else 'time step')
    fig.suptitle('Force coefficient history')
    fig.savefig(case/'forces.png', dpi=140)
    plt.close(fig)


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    for arg in sys.argv[1:]:
        case = Path(arg)
        plot_flow(case)
        plot_forces(case)
        print(f'{case}: flow.png, forces.png')
    return 0


if __name__ == '__main__':
    sys.exit(main())
