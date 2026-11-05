#!/usr/bin/env python3
"""Build every figure and table for the external-flow part of the report/presentation.

Usage (from code/):  python3 scripts/postprocess_report.py [runs] [../results]
  arg 1: folder with the solver output (default runs/; ../results/data also works)
  arg 2: results folder (default ../results); figures -> figures/, tables -> tables/

Figures (PNG, 200 dpi) and tables (Markdown + CSV):
  fig_polar.png              Cl(α), Cd(α), drag polar vs published Re = 1000 values
  fig_force_history.png      Cl(t), Cd(t) for every Re = 1000 angle (small multiples)
  fig_strouhal.png           shedding Strouhal number vs α, with Kurtulus (2015) values
  fig_vorticity_gallery.png  final vorticity field for every Re = 1000 angle
  fig_flow_aoa*.png          velocity/streamlines + Cp for selected angles
  fig_algorithms.png         SIMPLE vs PIMPLE at Re = 500 (only if those runs are present)
  fig_grid_study.png         coarse vs medium grid at α = 10° vs reference (if medium has finished)
  tables.md / *.csv          every number quoted in the figures
"""
from __future__ import annotations
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from make_vortex_video import Series  # noqa: E402

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else 'runs')
RESULTS = Path(sys.argv[2] if len(sys.argv) > 2 else '../results')
OUT, TAB = RESULTS/'figures', RESULTS/'tables'
OUT.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
RE1000 = [0, 4, 8, 10, 11, 12, 15, 20]
INK, MUTED, GRID = '#1f1e1c', '#6b6a64', '#e4e3de'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
REF = {'cl10': 0.418, 'cd10': 0.165, 'cd0': 0.119}  # Kurtulus 2022 Table 1; Di Ilio et al. Table 1
# Kurtulus (2015) Sec. 3.3, values stated in the text: St = f c / U at the given angle.
KURTULUS_ST = [(8, 0.876), (9, 0.89), (23, 0.45)]
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'axes.titleweight': 'bold',
                     'axes.labelcolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED})


def style(ax, xlabel=None, ylabel=None, title=None):
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    if title: ax.set_title(title, loc='left', color=INK)
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'): ax.spines[s].set_color(GRID)
    ax.tick_params(length=0)


def forces(case: Path):
    # skip a half-written last line (a run that was stopped mid-write)
    rows = [r for r in csv.DictReader(open(case/'forces.csv')) if r.get('cd_viscous') not in (None, '')]
    return {k: np.array([float(r[k]) for r in rows]) for k in ('step', 'time', 'cd', 'cl', 'cd_pressure', 'cd_viscous')}


def summary(case: Path) -> dict:
    return json.loads((case/'summary.json').read_text())


def window_stats(f, t_from, t_to=np.inf):
    m = (f['time'] >= t_from) & (f['time'] <= t_to)
    out = {}
    for k in ('cl', 'cd', 'cd_pressure', 'cd_viscous'):
        out[k], out[k+'_min'], out[k+'_max'] = f[k][m].mean(), f[k][m].min(), f[k][m].max()
    return out


def strouhal(f, t_from=40.0, min_amp=2e-3):
    """Shedding frequency from upward zero crossings of Cl - mean (needs >= 2 full cycles)."""
    m = f['time'] >= t_from
    t, y = f['time'][m], f['cl'][m]-f['cl'][m].mean()
    if (y.max()-y.min())/2 < min_amp:
        return None
    i = np.nonzero((y[:-1] < 0) & (y[1:] >= 0))[0]
    tc = t[i]-y[i]*(t[i+1]-t[i])/(y[i+1]-y[i])
    return None if len(tc) < 3 else 1.0/np.diff(tc).mean()


def write_table(name, header, rows, md):
    with open(TAB/f'{name}.csv', 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    md.append('| '+' | '.join(header)+' |')
    md.append('|'+'---|'*len(header))
    for r in rows:
        md.append('| '+' | '.join(f'{v:.4f}' if isinstance(v, float) else str(v) for v in r)+' |')
    md.append('')


def main():
    md = ['# External-flow results: tables for the report', '',
          'All values non-dimensional (chord c, U∞). Unsteady/unconverged runs: mean over t ≥ 40 (Re 1000) or the last third of the run (Re 500).', '']

    # ---- Re = 1000 sweep: statistics, Strouhal ------------------------------------------------
    sweep = []
    for a in RE1000:
        case = ROOT/f'naca0012_Re1000_aoa{a}_pimple'
        f, s = forces(case), summary(case)
        st = window_stats(f, 40.0 if f['time'][-1] > 45 else f['time'][-1]*2/3)
        sweep.append(dict(aoa=a, f=f, s=s, st=st, St=strouhal(f)))
    md += ['## Re = 1000 sweep (350×250, PIMPLE)', '']
    write_table('re1000_sweep', ['AoA', 'Cl mean', 'Cl min', 'Cl max', 'Cd mean', 'Cd pressure', 'Cd viscous', 'Cl/Cd', 'St', 'steps', 'converged', 'wall [s]'],
                [[d['aoa'], d['st']['cl'], d['st']['cl_min'], d['st']['cl_max'], d['st']['cd'], d['st']['cd_pressure'], d['st']['cd_viscous'],
                  d['st']['cl']/d['st']['cd'], d['St'] if d['St'] else '-', d['s']['steps'], d['s']['converged'], round(d['s']['wall_seconds'])]
                 for d in sweep], md)
    d10 = next(d for d in sweep if d['aoa'] == 10)
    d0 = next(d for d in sweep if d['aoa'] == 0)
    best = max(sweep, key=lambda d: d['st']['cl']/d['st']['cd'])
    best_ld = best['st']['cl']/best['st']['cd']
    md += ['## Validation against published Re = 1000 values', '']
    write_table('validation', ['Quantity', 'Reference', 'Source', 'Present', 'Error %'],
                [['Cd, AoA 0', REF['cd0'], 'Di Ilio et al. (2018) Tab. 1', d0['st']['cd'], 100*(d0['st']['cd']/REF['cd0']-1)],
                 ['Cl, AoA 10', REF['cl10'], 'Kurtulus (2022) Tab. 1', d10['st']['cl'], 100*(d10['st']['cl']/REF['cl10']-1)],
                 ['Cd, AoA 10', REF['cd10'], 'Kurtulus (2022) Tab. 1', d10['st']['cd'], 100*(d10['st']['cd']/REF['cd10']-1)],
                 [f"max Cl/Cd (ref. at 11 deg, present at {best['aoa']} deg)", 2.55, 'Kurtulus (2015) Sec. 3.1', best_ld, 100*(best_ld/2.55-1)]],
                md)

    # ---- Fig: polar (reuse plot_polar.py so the figure stays in one place) -----------------------
    subprocess.run([sys.executable, str(HERE/'plot_polar.py'), str(ROOT), str(TAB)], check=True, stdout=subprocess.DEVNULL)
    shutil.move(TAB/'naca0012_polar.png', OUT/'fig_polar.png')

    # ---- Fig: force histories --------------------------------------------------------------------
    fig, axes = plt.subplots(2, 4, figsize=(16, 6.5), sharex=True)
    for ax, d in zip(axes.flat, sweep):
        f = d['f']; m = f['time'] >= 2
        ax.plot(f['time'][m], f['cl'][m], color=BLUE, lw=1.6, label='$C_l$')
        ax.plot(f['time'][m], f['cd'][m], color=ORANGE, lw=1.6, label='$C_d$')
        ax.axvspan(40, 60, color=GRID, alpha=0.5, lw=0)
        style(ax, title=f"α = {d['aoa']}°")
    for ax in axes[1]: ax.set_xlabel('t U∞/c')
    axes[0, 0].legend(frameon=False)
    fig.suptitle('Force-coefficient histories, NACA 0012, Re = 1000 (shaded: averaging window)', x=0.01, ha='left', color=INK, weight='bold')
    fig.tight_layout(); fig.savefig(OUT/'fig_force_history.png', dpi=200); plt.close(fig)

    # ---- Fig: Strouhal --------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.5))
    pts = [(d['aoa'], d['St']) for d in sweep if d['St']]
    ax.plot(*zip(*pts), '-o', color=BLUE, lw=2, ms=8, mec='white', mew=1.5, label='Present (350×250, PIMPLE)')
    ax.plot(*zip(*KURTULUS_ST), 'D', mfc='none', mec=INK, ms=9, mew=1.6, ls='none', label='Kurtulus (2015), values stated in text')
    style(ax, 'Angle of attack α (deg)', 'Strouhal number  St = f c / U∞', 'Vortex-shedding frequency')
    ax.set_ylim(0, 1); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(OUT/'fig_strouhal.png', dpi=200); plt.close(fig)

    # ---- Fig: vorticity gallery (last snapshot of each angle) -------------------------------------
    cmap = matplotlib.colormaps['berlin'] if 'berlin' in matplotlib.colormaps else matplotlib.colormaps['RdBu_r']
    fig, axes = plt.subplots(4, 2, figsize=(13, 9.0), gridspec_kw=dict(hspace=0.05, wspace=0.02, left=0.005, right=0.995, top=0.955, bottom=0.1))
    for ax, d in zip(axes.T.flat, sweep):
        s = Series(ROOT/f"naca0012_Re1000_aoa{d['aoa']}_pimple")
        k = len(s.frames)-1
        w = np.ma.masked_where(s.field(k, 'solid') > 0.5, s.field(k, 'vorticity'))
        cmap.set_bad('#d9d7cf')
        im = ax.imshow(w, origin='lower', extent=s.extent, cmap=cmap, vmin=-6, vmax=6, interpolation='bicubic')
        ax.set_xlim(-0.7, 4.3); ax.set_ylim(-1.1, 0.9)
        ax.set_xticks([]); ax.set_yticks([])
        ax.text(0.015, 0.92, f"α = {d['aoa']}°   t = {s.frames[k][1]:.1f}", transform=ax.transAxes, color='white', fontsize=12, weight='bold', va='top')
    cb = fig.colorbar(im, cax=fig.add_axes([0.3, 0.05, 0.4, 0.015]), orientation='horizontal')
    cb.set_label('vorticity  ω c / U∞')
    fig.suptitle('Instantaneous vorticity at the end of each run, NACA 0012, Re = 1000', x=0.01, y=0.995, ha='left', color=INK, weight='bold')
    fig.savefig(OUT/'fig_vorticity_gallery.png', dpi=200, bbox_inches='tight'); plt.close(fig)

    # ---- Fig: velocity/Cp plots (existing per-case plotter) ----------------------------------------
    for a in (4, 10, 15, 20):
        case = ROOT/f'naca0012_Re1000_aoa{a}_pimple'
        subprocess.run([sys.executable, str(HERE/'plot_airfoil.py'), str(case)], check=True, stdout=subprocess.DEVNULL)
        shutil.copy(case/'flow.png', OUT/f'fig_flow_aoa{a}.png')

    # ---- Fig + table: SIMPLE vs PIMPLE (Re 500), only if those runs exist -------------------------
    if all((ROOT/f'naca0012_Re500_aoa{a}{t}'/'summary.json').exists() for a in (0, 5, 10, 15) for t in ('', '_pimple')):
        algorithm_comparison(md)
    grid_study(md)
    (TAB/'tables.md').write_text('\n'.join(md))
    print('figures:', ', '.join(sorted(p.name for p in OUT.iterdir())))
    print('tables: ', ', '.join(sorted(p.name for p in TAB.iterdir())))


def algorithm_comparison(md):
    rows, fig = [], plt.figure(figsize=(14, 4.6))
    a1, a2, a3 = fig.subplots(1, 3)
    for alg, col, tag in (('simple', AQUA, ''), ('pimple', ORANGE, '_pimple')):
        for a in (0, 5, 10, 15):
            case = ROOT/f'naca0012_Re500_aoa{a}{tag}'
            s, f = summary(case), forces(case)
            n = 1 if s['converged'] else max(1, len(f['cl'])//3)  # converged: final value; else mean of last third
            rows.append([alg.upper(), a, s['steps'], s['converged'], round(s['wall_seconds']), s['pressure_iterations'], s['momentum_iterations'],
                         float(np.mean(f['cl'][-n:])), float(np.mean(f['cd'][-n:]))])
            if a == 5:
                r = np.genfromtxt(case/'residuals.csv', delimiter=',', names=True)
                a1.semilogy(r['step'], np.maximum(r['u_momentum_inf'], r['v_momentum_inf']), color=col, lw=1.6, label=alg.upper())
                a2.semilogy(r['step'], r['divergence_inf'], color=col, lw=1.6, label=alg.upper())
    style(a1, 'iteration / time step', 'max momentum residual (∞-norm)', 'Convergence, α = 5°')
    style(a2, 'iteration / time step', 'continuity residual (∞-norm)', 'Mass conservation, α = 5°')
    angles = [0, 5, 10, 15]
    x = np.arange(len(angles))
    for k, (alg, col) in enumerate((('SIMPLE', AQUA), ('PIMPLE', ORANGE))):
        wall = [r[4]/60 for r in rows if r[0] == alg]
        a3.bar(x+(k-0.5)*0.38, wall, 0.36, color=col, label=alg)
    a3.set_xticks(x, [f'{a}°' for a in angles])
    style(a3, 'Angle of attack', 'wall time [min]', 'Cost to finish the run')
    for ax in (a1, a2, a3): ax.legend(frameon=False)
    fig.suptitle('Pressure–velocity coupling: SIMPLE vs PIMPLE, NACA 0012, Re = 500, 350×250', x=0.01, ha='left', color=INK, weight='bold')
    fig.tight_layout(); fig.savefig(OUT/'fig_algorithms.png', dpi=200); plt.close(fig)
    md += ['## SIMPLE vs PIMPLE at Re = 500 (350×250)', '', 'SIMPLE: iterations; PIMPLE: time steps (dt = 0.05). Cl/Cd: final value if converged, else mean of the last third.',
           'Wall times are indicative only: the runs shared the CPU with a varying number of other jobs.', '']
    write_table('algorithms', ['Algorithm', 'AoA', 'steps', 'converged', 'wall [s]', 'pressure CG its', 'momentum its', 'Cl', 'Cd'], rows, md)



def grid_study(md):
    # coarse = the sweep case at 10 deg, medium = runs/grid_study/aoa10_h013_peclet.
    # The medium run may have been stopped early (no summary.json), so the grid size comes from
    # input.cfg and both grids are averaged over the same window: t = 40 up to the end of the shorter run.
    med = ROOT/'grid_study'/'aoa10_h013_peclet'
    coarse = ROOT/'naca0012_Re1000_aoa10_pimple'
    if not (med/'forces.csv').exists():
        md += ['## Grid study', '', '_No medium-grid data yet (runs/grid_study/aoa10_h013_peclet)._', '']
        return
    levels = [('Coarse 350×250 (h = 0.02)', coarse, BLUE, 0.02), ('Medium 525×375 (h = 0.0133)', med, ORANGE, 0.02/1.5)]
    t_end = min(forces(c)['time'][-1] for _, c, _, _ in levels)
    if t_end < 45:
        md += ['## Grid study', '', f'_Medium grid only reached t = {t_end:.1f}; need t >= 45 for a usable average._', '']
        return
    fig, (b1, b2, b3) = plt.subplots(1, 3, figsize=(15, 4.6))
    grow = []
    for lab, case, col, h in levels:
        f = forces(case)
        cfg = dict(l.split(' = ', 1) for l in (case/'input.cfg').read_text().splitlines() if ' = ' in l)
        st = window_stats(f, 40, t_end)
        done = (case/'summary.json').exists()
        grow.append([lab, int(cfg['nx']), int(cfg['ny']), h, st['cl'], st['cl_min'], st['cl_max'], st['cd'], st['cd_pressure'], st['cd_viscous'],
                     100*(st['cl']/REF['cl10']-1), 100*(st['cd']/REF['cd10']-1), f"{f['time'][-1]:.1f}" + ('' if done else ' (stopped)')])
        m = f['time'] >= 2
        b1.plot(f['time'][m], f['cl'][m], color=col, lw=1.6, label=lab)
        b2.plot(f['time'][m], f['cd'][m], color=col, lw=1.6, label=lab)
    for ax, key in ((b1, 'cl10'), (b2, 'cd10')):
        ax.axhline(REF[key], color=INK, ls='--', lw=1.2, label='Kurtulus (2022)')
        ax.axvspan(40, t_end, color=GRID, alpha=0.5, lw=0)
    style(b1, 't U∞/c', '$C_l$', 'Lift history, α = 10°'); style(b2, 't U∞/c', '$C_d$', 'Drag history, α = 10°')
    b1.legend(frameon=False, fontsize=9)
    hs = [g[3] for g in grow]
    for k, mk, col, lab in ((4, 'o', BLUE, '$C_l$'), (7, 's', ORANGE, '$C_d$')):
        b3.errorbar(hs, [g[k] for g in grow], yerr=[[g[k]-g[k+1] for g in grow], [g[k+2]-g[k] for g in grow]] if k == 4 else None,
                    fmt='-'+mk, color=col, lw=2, ms=8, capsize=4, label=lab)
    b3.axhline(REF['cl10'], color=BLUE, ls='--', lw=1); b3.axhline(REF['cd10'], color=ORANGE, ls='--', lw=1)
    b3.text(0.0235, REF['cl10']+0.008, 'ref. Cl 0.418', color=BLUE, fontsize=9)
    b3.text(0.0235, REF['cd10']+0.008, 'ref. Cd 0.165', color=ORANGE, fontsize=9)
    b3.set_xlim(0.024, 0.008)
    style(b3, 'grid spacing h / c  (refinement →)', f'mean over t = 40–{t_end:.1f}', 'Grid convergence (dashed: reference)')
    b3.legend(frameon=False, loc='center right')
    fig.suptitle('Grid refinement, NACA 0012, Re = 1000, α = 10°', x=0.01, ha='left', color=INK, weight='bold')
    fig.tight_layout(); fig.savefig(OUT/'fig_grid_study.png', dpi=200); plt.close(fig)
    md += [f'## Grid study, Re = 1000, α = 10° (both grids averaged over t = 40–{t_end:.1f})', '',
           'Reference: Kurtulus (2022) Cl = 0.418, Cd = 0.165. Cl min/max show the shedding oscillation.', '']
    write_table('grid_study', ['Grid', 'nx', 'ny', 'h', 'Cl', 'Cl min', 'Cl max', 'Cd', 'Cd pressure', 'Cd viscous', 'Cl error %', 'Cd error %', 'run ends at t'], grow, md)


if __name__ == '__main__':
    main()
