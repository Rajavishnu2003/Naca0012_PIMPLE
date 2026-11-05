#!/usr/bin/env python3
"""Render vortex-shedding animations (MP4) from the solver's VTK snapshot series.

Single case (vorticity field + live lift-coefficient trace):
    python3 scripts/make_vortex_video.py results_external/video/naca0012_Re1000_aoa20 -o aoa20.mp4
Montage of several cases side by side (one panel per case):
    python3 scripts/make_vortex_video.py CASE1 CASE2 ... --montage -o montage.mp4

Each case directory needs flow.vtk.series + vtk/*.vtk (written with vtk_every > 0), summary.json
or input.cfg for the angle of attack, and forces.csv for the Cl trace.
ffmpeg: taken from --ffmpeg, else from the imageio-ffmpeg package if installed, else 'ffmpeg' on PATH.
"""
from __future__ import annotations
import argparse
import csv
import json
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import ListedColormap

INK_BG = '#0e0f13'
BG, INK, MUTED, BODY = INK_BG, '#ecebe6', '#9a988f', '#d9d7cf'
CL_LINE = '#3987e5'


class Series:
    """Frames of one case, read lazily from binary legacy VTK (big-endian float32)."""

    def __init__(self, case: Path):
        self.case = case
        meta = json.loads((case/'flow.vtk.series').read_text())
        self.frames = [(case/f["name"], float(f['time'])) for f in meta['files']]
        head = self.frames[0][0].read_bytes()[:600].decode('latin-1')
        self.nx, self.ny = map(int, re.search(r'DIMENSIONS (\d+) (\d+)', head).groups())
        x0, y0 = map(float, re.search(r'ORIGIN (\S+) (\S+)', head).groups())
        dx, dy = map(float, re.search(r'SPACING (\S+) (\S+)', head).groups())
        self.extent = (x0-dx/2, x0+(self.nx-0.5)*dx, y0-dy/2, y0+(self.ny-0.5)*dy)
        cfg = (case/'input.cfg').read_text()
        self.aoa = float(re.search(r'^aoa_deg = (\S+)', cfg, re.M).group(1))
        self.re = float(re.search(r'^reynolds = (\S+)', cfg, re.M).group(1))
        rows = [r for r in csv.DictReader(open(case/'forces.csv')) if r.get('cl')]
        self.t = np.array([float(r['time']) for r in rows])
        self.cl = np.array([float(r['cl']) for r in rows])
        self.cd = np.array([float(r['cd']) for r in rows])

    def field(self, k: int, name: str) -> np.ndarray:
        raw = self.frames[k][0].read_bytes()
        tag = f'SCALARS {name} float 1\nLOOKUP_TABLE default\n'.encode()
        i = raw.index(tag)+len(tag)
        n = self.nx*self.ny
        return np.frombuffer(raw, '>f4', n, i).reshape(self.ny, self.nx).astype(float)

    def coeffs_at(self, time: float) -> tuple[float, float]:
        j = int(np.abs(self.t-time).argmin())
        return self.cl[j], self.cd[j]


def body_layer(s: Series) -> np.ndarray:
    solid = s.field(0, 'solid') > 0.5
    rgba = np.zeros(solid.shape+(4,))
    rgba[solid] = matplotlib.colors.to_rgba(BODY)
    return rgba


def zoom(s: Series, aspect: float):
    solid = s.field(0, 'solid') > 0.5
    js, is_ = np.nonzero(solid)
    x0, _, y0, _ = s.extent
    dx = (s.extent[1]-s.extent[0])/s.nx
    dy = (s.extent[3]-s.extent[2])/s.ny
    xb = x0+np.array([is_.min(), is_.max()+1])*dx
    yc = y0+0.5*(js.min()+js.max()+1)*dy
    xl = (max(s.extent[0], xb[0]-0.7), min(s.extent[1], xb[1]+3.6))
    half = (xl[1]-xl[0])/aspect/2  # window matches the panel's width/height ratio
    return xl, (max(s.extent[2], yc-half), min(s.extent[3], yc+half))


def cmap():
    # Diverging map whose zero is the background, so irrotational flow disappears and only vortices glow.
    try:
        return matplotlib.colormaps['berlin']
    except KeyError:
        return ListedColormap(plt.cm.RdBu_r(np.linspace(0, 1, 256)))


def set_background():
    global BG
    BG = matplotlib.colors.to_hex(cmap()(0.5))


def style_axes(ax):
    ax.set_facecolor(BG)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def render_single(s: Series, out: Path, fps: int, wlim: float, ffmpeg_path: str, dpi: int = 120):
    plt.rcParams['animation.ffmpeg_path'] = ffmpeg_path
    fig = plt.figure(figsize=(16, 9), facecolor=BG)
    ax = fig.add_axes([0.02, 0.27, 0.96, 0.66])
    ac = fig.add_axes([0.07, 0.06, 0.88, 0.16])
    style_axes(ax)
    xl, yl = zoom(s, (0.96*16)/(0.66*9))
    im = ax.imshow(s.field(0, 'vorticity'), origin='lower', extent=s.extent, cmap=cmap(),
                   vmin=-wlim, vmax=wlim, interpolation='bicubic')
    ax.imshow(body_layer(s), origin='lower', extent=s.extent, interpolation='nearest')
    ax.set_xlim(xl); ax.set_ylim(yl)
    fig.text(0.02, 0.955, f'NACA 0012  ·  Re = {s.re:g}  ·  α = {s.aoa:g}°', color=INK, fontsize=22, weight='bold')
    cax = fig.add_axes([0.80, 0.958, 0.17, 0.014])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[-wlim, 0, wlim])
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=MUTED, labelsize=11, length=0)
    fig.text(0.785, 0.958, 'vorticity ωc/U∞', color=MUTED, fontsize=13, ha='right', va='bottom')
    tlabel = ax.text(0.012, 0.04, '', transform=ax.transAxes, color=INK, fontsize=16, family='monospace')
    t0, t1 = s.frames[0][1], s.frames[-1][1]
    m = (s.t >= t0-1e-9) & (s.t <= t1+1e-9)
    ac.set_facecolor(BG)
    ac.plot(s.t[m], s.cl[m], color=MUTED, lw=1.2, alpha=0.5)
    trace, = ac.plot([], [], color=CL_LINE, lw=2.2)
    dot, = ac.plot([], [], 'o', color=INK, ms=7)
    ac.set_xlim(t0, t1)
    pad = max(1e-3, 0.15*(s.cl[m].max()-s.cl[m].min()))
    ac.set_ylim(s.cl[m].min()-pad, s.cl[m].max()+pad)
    ac.set_ylabel('$C_l$', color=MUTED, fontsize=14)
    ac.set_xlabel('t U∞ / c', color=MUTED, fontsize=12)
    ac.tick_params(colors=MUTED, labelsize=11)
    for k, sp in ac.spines.items():
        sp.set_color('#3a3b40'); sp.set_visible(k in ('left', 'bottom'))
    writer = FFMpegWriter(fps=fps, bitrate=12000, codec='libx264',
                          extra_args=['-pix_fmt', 'yuv420p', '-preset', 'slow', '-crf', '16'])
    with writer.saving(fig, str(out), dpi):
        for k, (_, time) in enumerate(s.frames):
            im.set_data(s.field(k, 'vorticity'))
            cl, cd = s.coeffs_at(time)
            tlabel.set_text(f't = {time:6.2f}   Cl = {cl:6.3f}   Cd = {cd:6.3f}')
            sel = m & (s.t <= time+1e-9)
            trace.set_data(s.t[sel], s.cl[sel])
            dot.set_data([time], [cl])
            writer.grab_frame(facecolor=BG)
    plt.close(fig)


def render_montage(cases: list[Series], out: Path, fps: int, wlim: float, ffmpeg_path: str, dpi: int = 120):
    plt.rcParams['animation.ffmpeg_path'] = ffmpeg_path
    n = len(cases)
    cols = 2 if n <= 4 else 3
    rows = -(-n//cols)
    fig = plt.figure(figsize=(16, 9), facecolor=BG)
    top, gap = 0.90, 0.012
    w, h = (1-gap*(cols+1))/cols, (top-gap*(rows+1))/rows
    nframes = min(len(s.frames) for s in cases)
    ims, labels = [], []
    for k, s in enumerate(cases):
        r, c = divmod(k, cols)
        ax = fig.add_axes([gap+c*(w+gap), top-(r+1)*(h+gap), w, h])
        style_axes(ax)
        xl, yl = zoom(s, (w*16)/(h*9))
        ims.append(ax.imshow(s.field(0, 'vorticity'), origin='lower', extent=s.extent, cmap=cmap(),
                             vmin=-wlim, vmax=wlim, interpolation='bicubic'))
        ax.imshow(body_layer(s), origin='lower', extent=s.extent, interpolation='nearest')
        ax.set_xlim(xl); ax.set_ylim(yl)
        ax.text(0.02, 0.93, f'α = {s.aoa:g}°', transform=ax.transAxes, color=INK, fontsize=18, weight='bold', va='top')
        labels.append(ax.text(0.02, 0.05, '', transform=ax.transAxes, color=MUTED, fontsize=12, family='monospace'))
    fig.text(0.012, 0.945, f'NACA 0012  ·  Re = {cases[0].re:g}  ·  vorticity ωc/U∞', color=INK, fontsize=20, weight='bold')
    cax = fig.add_axes([0.62, 0.95, 0.16, 0.014])
    cb = fig.colorbar(ims[0], cax=cax, orientation='horizontal', ticks=[-wlim, 0, wlim])
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=MUTED, labelsize=11, length=0)
    tglobal = fig.text(0.988, 0.945, '', color=MUTED, fontsize=16, ha='right', family='monospace')
    writer = FFMpegWriter(fps=fps, bitrate=12000, codec='libx264',
                          extra_args=['-pix_fmt', 'yuv420p', '-preset', 'slow', '-crf', '16'])
    with writer.saving(fig, str(out), dpi):
        for k in range(nframes):
            for s, im, lab in zip(cases, ims, labels):
                time = s.frames[k][1]
                im.set_data(s.field(k, 'vorticity'))
                cl, cd = s.coeffs_at(time)
                lab.set_text(f'Cl {cl:5.3f}  Cd {cd:5.3f}')
            tglobal.set_text(f't = {cases[0].frames[k][1]:6.2f}')
            writer.grab_frame(facecolor=BG)
    plt.close(fig)


def default_ffmpeg() -> str:
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        return 'ffmpeg'


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cases', nargs='+', type=Path)
    p.add_argument('-o', '--out', type=Path, required=True)
    p.add_argument('--montage', action='store_true', help='one panel per case in a single video')
    p.add_argument('--fps', type=int, default=25)
    p.add_argument('--wlim', type=float, default=6.0, help='colour range: vorticity in [-wlim, wlim]')
    p.add_argument('--ffmpeg', default=None)
    a = p.parse_args()
    ff = a.ffmpeg or default_ffmpeg()
    set_background()
    series = [Series(c) for c in a.cases]
    if a.montage:
        render_montage(series, a.out, a.fps, a.wlim, ff)
    else:
        render_single(series[0], a.out, a.fps, a.wlim, ff)
    print(f'wrote {a.out}')


if __name__ == '__main__':
    main()
