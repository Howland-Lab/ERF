#!/usr/bin/env python3
"""
gz_filter_compare.py -- the fair test for the grey-zone scan.

A model with grid spacing Delta should reproduce the fine LES averaged over Delta-boxes (that is all a Delta grid
can represent).  So compare, at 0.5 zi over 3-4 h:
  model    : ERF run at Delta (gz_cbl_dx050 ... dx400), w on its own grid
  target   : the 25 m ERF runs (6 realisations: cases/cbl_sd25, cases_pvc/cbl_sd25, cases_pvc/ens_cbl25_s1..4),
             w averaged over boxes of side Delta (2, 4, 8, 16 cells)
Statistics: sigma_w (absolute, m/s), skewness, P(w' > 3 sigma) with sigma the field's own std.
The spread of the target across the 6 realisations is the sampling uncertainty.
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from erf_planes import Plotfile

ROOT = os.path.join(HERE, "..")
FIG = os.path.join(HERE, "figs")
GUIDE = "/work2/11859/hsna/stampede3/HowlandLab/Results Summary/figs"
ZI = 1040.0; K = int(round(0.5 * ZI / 20.0 - 0.5))          # cell index of 0.5 zi on the 20 m vertical grid
WIN = (10800, 14400)


def stats(samples):
    x = np.concatenate([s.ravel() for s in samples]); x = x - x.mean(); s = x.std()
    return dict(std=float(s), skew=float(np.mean(x ** 3) / s ** 3), p3=float(np.mean(x > 3 * s)), n=int(x.size))


def fields(path):
    pls = [Plotfile(os.path.join(path, p)) for p in sorted(os.listdir(path)) if p.startswith("plt")]
    pls = [p for p in pls if WIN[0] - 1 <= p.time <= WIN[1] + 1]
    return [p.plane("z_velocity", K) for p in pls]


fine_runs = ["cases/cbl_sd25", "cases_pvc/cbl_sd25"] + [f"cases_pvc/ens_cbl25_s{s}" for s in range(1, 5)]
boxes = {50: 2, 100: 4, 200: 8, 400: 16}
target = {L: [] for L in [25] + list(boxes)}
for r in fine_runs:
    W = fields(os.path.join(ROOT, r))
    target[25].append(stats(W))
    for L, n in boxes.items():
        Wb = []
        for w in W:
            ny, nx = w.shape; my, mx = (ny // n) * n, (nx // n) * n
            Wb.append(w[:my, :mx].reshape(my // n, n, mx // n, n).mean(axis=(1, 3)))
        target[L].append(stats(Wb))
    print("fine", r, "done", flush=True)
model = {}
for L in boxes:
    model[L] = stats(fields(os.path.join(ROOT, "cases", f"gz_cbl_dx{L:03d}")))
model[25] = target[25][0]

Ls = [25, 50, 100, 200, 400]
fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.9))
for a, key, lab, logy in ((ax[0], "std", r"$\sigma_w$ at 0.5 $z_i$ (m/s)", False), (ax[1], "skew", "skewness of w", False),
                          (ax[2], "p3", r"P($w' > 3\sigma_w$)", True)):
    tm = np.array([np.mean([t[key] for t in target[L]]) for L in Ls]); ts = np.array([np.std([t[key] for t in target[L]], ddof=1) for L in Ls])
    a.fill_between(np.array(Ls) / ZI, tm - ts, tm + ts, color="#BBBBBB", alpha=0.5, label="25 m LES averaged over Δ-boxes (±1 std, 6 runs)")
    a.plot(np.array(Ls) / ZI, tm, "-", color="#555555")
    a.plot(np.array(Ls) / ZI, [model[L][key] for L in Ls], "o-", color="#B2182B", lw=1.8, label="ERF run with grid Δ (same LES closure)")
    a.set_xscale("log"); a.set_xlabel(r"$\Delta / z_i$"); a.set_ylabel(lab); a.grid(alpha=0.25, which="both")
    a.axvspan(0.1, 0.5, color="#F3E9D2", zorder=0)
    if logy:
        a.set_yscale("log")
ax[0].legend(frameon=False, fontsize=7.5, loc="lower left")
fig.suptitle("Fair test: does a coarse grid reproduce what the fine LES looks like at that grid? (convective, 0.5 zi, 3–4 h)", fontsize=11)
fig.tight_layout(); fig.savefig(f"{FIG}/w2_gz_fairtest.png", dpi=160); plt.close(fig)
os.system(f"cp '{FIG}/w2_gz_fairtest.png' '{GUIDE}/'")
out = {str(L): dict(model=model[L], target_mean={k: float(np.mean([t[k] for t in target[L]])) for k in ("std", "skew", "p3")},
                    target_std={k: float(np.std([t[k] for t in target[L]], ddof=1)) for k in ("std", "skew", "p3")}) for L in Ls}
json.dump(out, open(os.path.join(HERE, "gz_fairtest.json"), "w"), indent=1)
for L in Ls:
    o = out[str(L)]
    print(f"Delta={L:3d} m: model std {o['model']['std']:.3f} skew {o['model']['skew']:.2f} p3 {o['model']['p3']:.1e} | "
          f"filtered LES std {o['target_mean']['std']:.3f}±{o['target_std']['std']:.3f} skew {o['target_mean']['skew']:.2f}±{o['target_std']['skew']:.2f} "
          f"p3 {o['target_mean']['p3']:.1e}±{o['target_std']['p3']:.1e}")
