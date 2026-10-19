#!/usr/bin/env python3
"""
erf_planes.py -- distributions of w, u, theta on horizontal planes from ERF plotfiles, computed exactly like
postproc_highres.py does for the WRF planes (same heights relative to zi / h, same Moments accumulator,
same histogram bins, tail counts beyond 2 and 3 sigma of the first frame), then compared with WRF.

    python3 erf_planes.py          # -> planes.json, figs/planes_*.png, printed table

ERF plotfiles (plot_per = 600 s) hold cell-centred x_velocity, y_velocity, z_velocity, theta.  Single level,
no AMR.  Last hour only: CBL 3-4 h (7 files), SBL 8-9 h (7 files).  WRF uses 31 frames (every 2 min), so the
ERF samples are about 4x smaller; moments up to kurtosis are still well converged at these sample sizes
(10^5-10^6 points per level), the 3-sigma tails less so.
"""
import json, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ERF = os.path.join(HERE, "..", "cases")
WRFP = "/work2/11859/hsna/stampede3/HowlandLab/WRF-4.6.1/test/abl_cases/postproc"
sys.path.insert(0, WRFP)
FIG = os.path.join(HERE, "figs"); os.makedirs(FIG, exist_ok=True)
HBINS = dict(w=np.linspace(-8, 8, 321), u=np.linspace(-10, 10, 401), th=np.linspace(-4, 4, 321))
ZREL = {"cbl": (0.1, 0.5, 0.9), "sbl": (0.25, 0.5, 0.75)}
ZHUB = {"cbl": 100.0, "sbl": 50.0}
CASES = {
    # ERF run       case   window (s)        WRF runs to compare (valid ones first)
    "cbl_sd25":    ("cbl", (10800, 14400), ["cbl_dx30", "cbl_dx10", "cbl_sd25"]),
    "sbl_dx6p25":  ("sbl", (28800, 32400), ["sbl_dx6p25"]),
    "sbl_dx3p125": ("sbl", (28800, 32400), ["sbl_dx3p125"]),
}


class Moments:                                            # same as postproc_highres.py
    def __init__(self):
        self.n = 0; self.s = np.zeros(4); self.tail = np.zeros(4); self.hist = None

    def add(self, x, bins, sigma_ref=None):
        x = x.ravel().astype(np.float64)
        self.n += x.size
        self.s += [x.sum(), (x ** 2).sum(), (x ** 3).sum(), (x ** 4).sum()]
        h, _ = np.histogram(x, bins=bins)
        self.hist = h if self.hist is None else self.hist + h
        if sigma_ref is not None:
            self.tail += [np.sum(x > 2 * sigma_ref), np.sum(x < -2 * sigma_ref), np.sum(x > 3 * sigma_ref), np.sum(x < -3 * sigma_ref)]

    def result(self):
        n = max(self.n, 1); m1 = self.s[0] / n; m2 = self.s[1] / n - m1 ** 2
        m3 = self.s[2] / n - 3 * m1 * self.s[1] / n + 2 * m1 ** 3
        m4 = self.s[3] / n - 4 * m1 * self.s[2] / n + 6 * m1 ** 2 * self.s[1] / n - 3 * m1 ** 4
        sd = np.sqrt(max(m2, 1e-30))
        return dict(n=n, mean=m1, std=sd, skew=m3 / sd ** 3, kurt=m4 / sd ** 4, hist=self.hist, tail=self.tail / n)


class Plotfile:
    """Minimal reader for a single-level AMReX plotfile (HyperCLaw-V1.1, 8-byte little-endian reals)."""
    def __init__(self, path):
        self.path = path
        L = open(f"{path}/Header").read().split("\n")
        nv = int(L[1]); self.vars = L[2:2 + nv]; i = 2 + nv
        self.time = float(L[i + 1]); i += 3
        self.lo = [float(x) for x in L[i].split()]; self.hi = [float(x) for x in L[i + 1].split()]
        dom = re.findall(r"\((\d+),(\d+),(\d+)\) \((\d+),(\d+),(\d+)\)", L[i + 3])[0]
        self.n = [int(dom[3]) + 1, int(dom[4]) + 1, int(dom[5]) + 1]
        self.dx = [(self.hi[d] - self.lo[d]) / self.n[d] for d in range(3)]
        H = open(f"{path}/Level_0/Cell_H").read()
        self.boxes = [tuple(map(int, b)) for b in re.findall(r"\(\((\d+),(\d+),(\d+)\) \((\d+),(\d+),(\d+)\) \(0,0,0\)\)", H)]
        self.fabs = re.findall(r"FabOnDisk: (\S+) (\d+)", H)
        assert len(self.boxes) == len(self.fabs)

    def plane(self, var, k):
        """2-D field (ny, nx) of variable var at cell index k"""
        c = self.vars.index(var); out = np.empty((self.n[1], self.n[0]))
        for (i0, j0, k0, i1, j1, k1), (fn, off) in zip(self.boxes, self.fabs):
            if not (k0 <= k <= k1):
                continue
            with open(f"{self.path}/Level_0/{fn}", "rb") as f:
                f.seek(int(off)); hdr = f.readline().decode()
                ncomp = int(hdr.split()[-1]); assert "(8 7 6 5 4 3 2 1)" in hdr, hdr
                nx, ny, nz = i1 - i0 + 1, j1 - j0 + 1, k1 - k0 + 1
                start = f.tell() + 8 * (c * nx * ny * nz + (k - k0) * nx * ny)
                f.seek(start); a = np.fromfile(f, dtype="<f8", count=nx * ny).reshape(ny, nx)
            out[j0:j1 + 1, i0:i1 + 1] = a
        return out


def erf_stats(name, case, win):
    d = os.path.join(ERF, name)
    pls = sorted((Plotfile(f"{d}/{p}") for p in os.listdir(d) if p.startswith("plt")), key=lambda p: p.time)
    pls = [p for p in pls if win[0] - 1 <= p.time <= win[1] + 1]
    summ = json.load(open(os.path.join(HERE, "summary.json")))["pairs"][name]["ERF"]
    zi = summ["h_minflux"] if case == "cbl" else summ["h_stress"]
    dz = pls[0].dx[2]; zc = (np.arange(pls[0].n[2]) + 0.5) * dz
    zlist = [(f"{f:.2f}", f * zi) for f in ZREL[case]] + [("hub", ZHUB[case])]
    klist = [int(np.argmin(abs(zc - z))) for _, z in zlist]
    sig = {}
    for (lab, _), k in zip(zlist, klist):
        w = pls[0].plane("z_velocity", k); sig[lab] = (w - w.mean()).std()
    mom = {lab: {v: Moments() for v in ("w", "u", "th")} for lab, _ in zlist}
    for p in pls:
        for (lab, _), k in zip(zlist, klist):
            w = p.plane("z_velocity", k); w -= w.mean()
            u = p.plane("x_velocity", k); u -= u.mean()
            th = p.plane("theta", k); th -= th.mean()
            mom[lab]["w"].add(w, HBINS["w"], sig[lab]); mom[lab]["u"].add(u, HBINS["u"]); mom[lab]["th"].add(th, HBINS["th"])
        print(f"  {name} t={p.time / 3600:.2f} h done", flush=True)
    out = dict(zlist=[(lab, float(z), float(zc[k])) for (lab, z), k in zip(zlist, klist)], zi=zi, nframes=len(pls))
    for lab, _ in zlist:
        for v in ("w", "u", "th"):
            out[f"{lab}_{v}"] = mom[lab][v].result()
    return out


def wrf_planes(name):
    d = np.load(f"{WRFP}/data_hr/{name}.npz", allow_pickle=True)
    return d["planes"].item()


def main():
    res = {}
    lines = []
    for name, (case, win, wruns) in CASES.items():
        E = erf_stats(name, case, win)
        Ws = {w: wrf_planes(w) for w in wruns}
        res[name] = dict(ERF=E, WRF={w: {k: v for k, v in P.items() if k[:4] in ("0.10", "0.25", "0.50", "0.75", "0.90", "hub_")} for w, P in Ws.items()})
        labs = [l for l, _, _ in E["zlist"]]
        lines.append(f"== ERF {name} ({E['nframes']} frames)  vs  WRF {', '.join(wruns)}")
        for lab in labs:
            for v in ("w", "u"):
                e = E[f"{lab}_{v}"]
                row = f"   {lab:4s} {v}: ERF std {e['std']:.3f} skew {e['skew']:+.2f} kurt {e['kurt']:.2f}"
                if v == "w":
                    row += f" P(>3s) {e['tail'][2]:.1e} P(<-3s) {e['tail'][3]:.1e}"
                for w, P in Ws.items():
                    q = P.get(f"{lab}_{v}")
                    if q is not None:
                        row += f" | {w.split('_')[1]}: {q['std']:.3f} {q['skew']:+.2f} {q['kurt']:.2f}"
                        if v == "w":
                            row += f" {q['tail'][2]:.1e}"
                lines.append(row)
        # figure: w pdfs at the levels, ERF vs WRF valid runs
        fig, ax = plt.subplots(1, len(labs), figsize=(3.6 * len(labs), 3.4), sharey=True)
        c = 0.5 * (HBINS["w"][1:] + HBINS["w"][:-1]); dw = HBINS["w"][1] - HBINS["w"][0]
        cols = {"cbl_dx30": "#555555", "cbl_dx10": "#222222", "cbl_sd25": "#BBBBBB", "sbl_dx6p25": "#222222", "sbl_dx3p125": "#222222"}
        for j, lab in enumerate(labs):
            e = E[f"{lab}_w"]; s = e["std"]
            ax[j].semilogy(c / s, e["hist"] / e["n"] / dw * s, color="#9B2335", lw=1.6, label=f"ERF ({name.split('_')[1]})")
            for w, P in Ws.items():
                q = P.get(f"{lab}_w")
                if q is None:
                    continue
                ls = ":" if w == "cbl_sd25" else "-"
                ax[j].semilogy(c / q["std"], q["hist"] / q["n"] / dw * q["std"], ls, color=cols[w], lw=1.0,
                               label=f"WRF {w.split('_')[1]}" + (" (tiled, invalid)" if w == "cbl_sd25" else ""))
            g = np.exp(-0.5 * c ** 2) / np.sqrt(2 * np.pi); ax[j].semilogy(c, g, "--", color="#999999", lw=0.8, label="Gaussian")
            ax[j].set_xlim(-5, 7); ax[j].set_ylim(1e-6, 1); ax[j].grid(alpha=0.25)
            ax[j].set_title(f"{'z/zi' if case == 'cbl' else 'z/h'} = {lab}" if lab != "hub" else f"z = {ZHUB[case]:.0f} m", fontsize=10)
            ax[j].set_xlabel(r"$w'/\sigma_w$")
        ax[0].set_ylabel("PDF"); ax[0].legend(frameon=False, fontsize=7)
        fig.suptitle(f"Vertical-velocity distributions, {name}, last hour", fontsize=11)
        fig.tight_layout(); fig.savefig(f"{FIG}/planes_w_{name}.png", dpi=150); plt.close(fig)
    txt = "\n".join(lines); print(txt)
    open(os.path.join(HERE, "planes_summary.txt"), "w").write(txt + "\n")

    def clean(o):
        if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [clean(v) for v in o]
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.generic): return o.item()
        return o
    json.dump(clean(res), open(os.path.join(HERE, "planes.json"), "w"))


if __name__ == "__main__":
    main()
