#!/usr/bin/env python3
"""Compare NACA 0012 force coefficients across runs: Cl(AoA), Cd(AoA) and the drag polar.

Usage: python3 scripts/plot_polar.py [runs] [out_dir]
Writes <out_dir>/naca0012_polar.png and <out_dir>/naca0012_polar.csv (out_dir defaults to the runs folder).

For each run the force history (forces.csv) is reduced as follows:
  * converged run          -> final value, no range bar
  * not converged (SIMPLE stalled, or PIMPLE still oscillating at max_steps)
                           -> mean over the last third of the run, with a bar for the min..max range
The bar shows how much the coefficient was still moving (vortex shedding at high AoA),
not a statistical uncertainty.
"""
from __future__ import annotations
import csv
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# (label, Re, algorithm, colour, marker) - colours are categorical slots 1-3 of the reference palette.
SERIES = [
    ('Re 1000, PIMPLE', 1000, 'pimple', '#2a78d6', 'o'),
    ('Re 500, PIMPLE', 500, 'pimple', '#eb6834', 's'),
    ('Re 500, SIMPLE', 500, 'simple', '#1baf7a', '^'),
]
# Published Re = 1000 values (text/tables only), plotted as black open markers.
REF_FILE = Path(__file__).resolve().parent.parent/'data'/'naca0012_re1000_reference.csv'
REF_MARK = {'Kurtulus 2022': 'D', 'Mittal & Tezduyar 1994': 'P', 'Kurtulus 2015': 'X', 'Di Ilio et al. 2018': 'v'}
INK, MUTED, GRID = '#1f1e1c', '#6b6a64', '#e4e3de'


def reduce(case: Path) -> dict | None:
    sfile = case/'summary.json'
    if not sfile.exists():  # run never finished (e.g. the aborted *_fine cases)
        return None
    s = json.loads(sfile.read_text())
    if s.get('body') != 'naca' or s.get('nx') != 350:  # keep one grid so runs are comparable
        return None
    rows = [r for r in csv.DictReader(open(case/'forces.csv')) if r.get('cl')]
    window = rows if s['converged'] else rows[-max(1, len(rows)//3):]
    out = dict(case=case.name, Re=s['Re'], algorithm=s['algorithm'], aoa=s['aoa_deg'],
               steps=s['steps'], converged=s['converged'])
    for k in ('cd', 'cl'):
        v = [float(r[k]) for r in window]
        if s['converged']:
            out[k] = out[f'{k}_min'] = out[f'{k}_max'] = s[k]
        else:
            out[k], out[f'{k}_min'], out[f'{k}_max'] = sum(v)/len(v), min(v), max(v)
    return out


def style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, color=MUTED)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, length=0)


def main() -> None:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else 'runs')
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else root
    out.mkdir(parents=True, exist_ok=True)
    data = [d for d in (reduce(c) for c in sorted(root.glob('naca0012_*'))) if d]
    with open(out/'naca0012_polar.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(data[0]))
        w.writeheader()
        for d in sorted(data, key=lambda d: (d['Re'], d['algorithm'], d['aoa'])):
            w.writerow(d)

    fig, (acl, acd, apol) = plt.subplots(1, 3, figsize=(16, 5.2))
    for label, re, alg, col, mk in SERIES:
        pts = sorted((d for d in data if d['Re'] == re and d['algorithm'] == alg), key=lambda d: d['aoa'])
        if not pts:
            continue
        a = [d['aoa'] for d in pts]
        for ax, k in ((acl, 'cl'), (acd, 'cd')):
            y = [d[k] for d in pts]
            err = [[d[k]-d[f'{k}_min'] for d in pts], [d[f'{k}_max']-d[k] for d in pts]]
            ax.plot(a, y, '-', color=col, lw=2, zorder=2)
            ax.errorbar(a, y, yerr=err, fmt=mk, color=col, ms=8, mec='white', mew=1.5,
                        elinewidth=1.5, capsize=4, label=label, zorder=3)
        apol.plot([d['cd'] for d in pts], [d['cl'] for d in pts], '-', color=col, lw=2, zorder=2)
        apol.plot([d['cd'] for d in pts], [d['cl'] for d in pts], mk, color=col, ms=8,
                  mec='white', mew=1.5, label=label, zorder=3)

    if REF_FILE.exists():
        refs = list(csv.DictReader(l for l in open(REF_FILE) if not l.startswith('#')))
        for src in dict.fromkeys(r['source'] for r in refs):
            rs = [r for r in refs if r['source'] == src]
            kw = dict(ls='none', marker=REF_MARK.get(src, 'D'), ms=9, mfc='none', mec=INK, mew=1.6,
                      label=f'{src} (Re 1000, ref.)', zorder=4)
            acl.plot([float(r['aoa_deg']) for r in rs], [float(r['cl']) for r in rs], **kw)
            acd.plot([float(r['aoa_deg']) for r in rs], [float(r['cd']) for r in rs], **kw)
            apol.plot([float(r['cd']) for r in rs], [float(r['cl']) for r in rs], **kw)
    style(acl, 'Angle of attack α (deg)', 'Lift coefficient $C_l$')
    style(acd, 'Angle of attack α (deg)', 'Drag coefficient $C_d$')
    style(apol, 'Drag coefficient $C_d$', 'Lift coefficient $C_l$')
    acl.axhline(0, color=MUTED, lw=0.8)
    acl.set_title('Lift vs angle of attack', loc='left', color=INK)
    acd.set_title('Drag vs angle of attack', loc='left', color=INK)
    apol.set_title('Drag polar', loc='left', color=INK)
    acd.legend(frameon=False, loc='upper center', fontsize=9)
    fig.suptitle('NACA 0012, laminar, 350×250 staircase grid. Bars: min–max over the last third of unconverged runs',
                 x=0.01, ha='left', fontsize=10, color=MUTED)
    fig.tight_layout()
    fig.savefig(out/'naca0012_polar.png', dpi=150)
    print(f'wrote {out/"naca0012_polar.png"} and {out/"naca0012_polar.csv"} ({len(data)} runs)')


if __name__ == '__main__':
    main()
