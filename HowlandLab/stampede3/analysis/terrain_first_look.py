#!/usr/bin/env python3
"""terrain_first_look.py -- first look at the ERF terrain test cases/cbl_woa_dx100 (Witch-of-Agnesi ridge in x, hmax 200 m,
L 1 km, terrain-fitted mesh, submitted 22 Oct 2026): an x-z cross-section through the domain centre of u, w and theta from
the last plotfile, the plane-averaged profiles, and a few numbers (speed-up over the crest, wake behind it).

The fitted mesh is terrain-following: cell k sits at z = h(x) + (k + 1/2) * dz * (1 - h(x)/H) (ERF's basic terrain-
following mapping with domain height H).  That mapping is used for the plot; if the plotfile carries a node-height
variable it is listed, use it instead when precision matters.

Usage: python3 terrain_first_look.py [case_dir]      (default cases/cbl_woa_dx100; any ERF case works, flat = h = 0)
"""
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from erf_planes import Plotfile

H = "/work2/11859/hsna/stampede3/HowlandLab"
case = sys.argv[1] if len(sys.argv) > 1 else f"{H}/ERF_runs/cases/cbl_woa_dx100"
name = os.path.basename(case)


def inputs(case):
    d = {}
    for l in open(f"{case}/inputs"):
        l = l.split("#")[0].strip()
        if "=" in l: k, v = l.split("=", 1); d[k.strip()] = v.strip().strip('"')
    return d


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    inp = inputs(case); hmax = float(inp.get("prob.hmax", 0)); L = float(inp.get("prob.L", 1)); ext = [float(x) for x in inp["geometry.prob_extent"].split()]
    pls = sorted(glob.glob(f"{case}/plt*"), key=lambda p: Plotfile(p).time)
    if not pls: sys.exit(f"{name}: no plotfile yet")
    p = Plotfile(pls[-1]); print(f"{name}: plotfile {os.path.basename(pls[-1])} t = {p.time:.0f} s, vars = {p.vars}, n = {p.n}, dx = {p.dx}")
    nx, ny, nz = p.n; dz = p.dx[2]; x = (np.arange(nx) + 0.5) * p.dx[0]; xc = 0.5 * ext[0]
    h = hmax / (1 + ((x - xc) / L) ** 2) if hmax else np.zeros(nx)
    j = ny // 2
    sec = {v: np.array([p.plane(v, k)[j, :] for k in range(nz)]) for v in ("x_velocity", "z_velocity", "theta")}   # (nz, nx)
    Z = h[None, :] + (np.arange(nz)[:, None] + 0.5) * dz * (1 - h[None, :] / ext[2]); X = np.repeat(x[None, :], nz, 0)
    umean = sec["x_velocity"].mean(1)                                   # crude plane mean along x at this j (for the speed-up)
    kh = int(np.argmin(abs(np.arange(nz) * dz - 100)))                  # ~100 m above the local surface
    ic = int(np.argmin(abs(x - xc))); iu = int(np.argmin(abs(x - (xc - 3 * L))))
    print(f"  at ~100 m above ground: u over the crest {sec['x_velocity'][kh, ic]:.2f} m/s vs 3 L upstream {sec['x_velocity'][kh, iu]:.2f} m/s "
          f"-> speed-up {sec['x_velocity'][kh, ic] / sec['x_velocity'][kh, iu] - 1:+.0%} (Jackson-Hunt order 2 hmax/L = {2 * hmax / L:.2f} for a shallow hill)")
    neg = (sec["x_velocity"][:3, :] < 0).any(0); print(f"  reversed flow (u < 0) in the lowest 3 cells at {neg.sum()} of {nx} columns" + (f", x = {x[neg].min():.0f}-{x[neg].max():.0f} m" if neg.any() else ""))
    fig, ax = plt.subplots(3, 1, figsize=(10, 8.5), sharex=True)
    for a, v, lab, cmap, lim in zip(ax, ("x_velocity", "z_velocity", "theta"), ("u (m/s)", "w (m/s)", "θ (K)"), ("viridis", "RdBu_r", "magma"), (None, 3.0, None)):
        f = sec[v]; kw = dict(vmin=-lim, vmax=lim) if lim else {}
        m = a.pcolormesh(X / 1000, Z, f, cmap=cmap, shading="auto", **kw); fig.colorbar(m, ax=a, pad=0.01, label=lab)
        a.fill_between(x / 1000, 0, h, color="#8A7355"); a.set_ylim(0, min(ext[2], 1500)); a.set_ylabel("z (m)")
    ax[-1].set_xlabel("x (km)"); fig.suptitle(f"{name}: x–z section at y = {(j + 0.5) * p.dx[1] / 1000:.1f} km, t = {p.time / 3600:.1f} h (ridge hmax {hmax:.0f} m, L {L:.0f} m)", fontsize=11)
    fig.tight_layout(); out = f"{H}/Results Summary/figs/terrain_{name}_section.png"; fig.savefig(out, dpi=140); print("wrote", out)
