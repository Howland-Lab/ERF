#!/usr/bin/env python3
"""ar_test_analysis.py -- grid aspect-ratio test of the VLES scan (runs made by make_ar_runs.py, submitted 22 Oct 2026).
Mike (22 Oct): the LES SGS model is isotropic, so dx/dz matters.  The scan result "sigma_w +5/+11/+18 % at 100/200/400 m"
was obtained with dz = 20 m, i.e. aspect ratios 5/10/20.  Here dx = 100 m is held and dz = 20, 50, 100 m compared
(ratio 5, 2, 1).  If sigma_w at dx = 100 m changes with dz, part of the scan's excess is an aspect-ratio effect.

Reads the plane-averaged profile logs (resolved <w'w'> from covar_profiles, SGS TKE and tau_33 from the logs) of
  cases/gz_cbl_dx100          dz = 20 m   (the scan run)
  cases/gz_cbl_dx100_dz050    dz = 50 m
  cases/gz_cbl_dx100_dz100    dz = 100 m
and the fair-test target (25 m LES box-averaged to 100 m) from gz_fairtest.json ('100' -> target_mean.std, target_std.std).
Prints a table and writes Results Summary/figs/ar_test_sigw.png + Results Summary/ar_test.json.  Runs that have not
finished are skipped with a message, so the script can be run as the results come in.
"""
import json, os
import numpy as np

H = "/work2/11859/hsna/stampede3/HowlandLab"
RUNS = [("gz_cbl_dx100", 20), ("gz_cbl_dx100_dz050", 50), ("gz_cbl_dx100_dz100", 100)]
WIN = (10800, 14400); G, TH = 9.81, 300.0
COL = {20: "#0072B2", 50: "#009E73", 100: "#D55E00"}


def load(path):
    a = np.loadtxt(path); z = np.unique(a[:, 1]); nz = z.size; t = a[::nz, 0]
    return t, z, a.reshape(t.size, nz, -1)


def analyse(d):
    tc, z, C = load(f"{d}/covar_profiles.dat"); ts, _, S = load(f"{d}/sfs_profiles.dat"); tm, _, M = load(f"{d}/mean_profiles.dat")
    if tc[-1] < WIN[1] - 1: raise RuntimeError(f"{os.path.basename(d)}: only reached t = {tc[-1]:.0f} s")
    wc, ws, wm = (tc >= WIN[0]), (ts >= WIN[0]), (tm >= WIN[0])
    n = z.size - 1; dz = np.diff(z); zc = z[:-1] + 0.5 * dz
    ww_f = C[wc][:, :, 7].mean(0); ww = 0.5 * (ww_f[:-1] + ww_f[1:])                    # resolved, faces -> centres
    wth_tot = C[wc][:, :, 10].mean(0) + S[ws][:, :, 8].mean(0)
    zi = float(z[np.argmin(wth_tot)]); q0 = wth_tot[0]; wstar = (G / TH * q0 * zi) ** (1 / 3)
    ksgs = M[wm][:, :n, 7].mean(0); tau33 = 0.5 * (S[ws][:, :, 7].mean(0)[:-1] + S[ws][:, :, 7].mean(0)[1:])
    ml = (zc > 0.2 * zi) & (zc < 0.8 * zi)
    sig_res = np.sqrt(np.maximum(ww, 0)); sig_tot = np.sqrt(np.maximum(ww + 2 / 3 * ksgs, 0))                # isotropic share of SGS TKE
    uw_tot = C[wc][:, :, 4].mean(0) + S[ws][:, :, 4].mean(0); vw_tot = C[wc][:, :, 6].mean(0) + S[ws][:, :, 6].mean(0)
    ustar = float((uw_tot[0] ** 2 + vw_tot[0] ** 2) ** 0.25)
    return dict(z=zc.tolist(), zi=zi, wstar=float(wstar), ustar=ustar, sig_res=sig_res.tolist(), sig_tot=sig_tot.tolist(),
                ml_sig_res=float(sig_res[ml].mean()), ml_sig_tot=float(sig_tot[ml].mean()), ml_sgs_share=float((2 / 3 * ksgs[ml] / (ww[ml] + 2 / 3 * ksgs[ml])).mean()),
                ml_sig_res_over_wstar=float(sig_res[ml].mean() / wstar), dz=float(dz[0]))


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    ft = json.load(open(f"{H}/ERF_runs/analysis/gz_fairtest.json")); tgt = ft["100"]["target_mean"]["std"]; tgt_sd = ft["100"]["target_std"]["std"]
    res = {}
    for name, dz in RUNS:
        try:
            res[dz] = analyse(f"{H}/ERF_runs/cases/{name}")
        except Exception as e:
            print(f"skip {name}: {e}"); continue
        r = res[dz]
        print(f"dz = {dz:3d} m (dx/dz = {100 // dz:2d}): zi {r['zi']:.0f} m, w* {r['wstar']:.2f}, u* {r['ustar']:.3f} | mixed-layer sigma_w resolved {r['ml_sig_res']:.3f} m/s "
              f"({100 * (r['ml_sig_res'] / tgt - 1):+.1f} % vs the fair-test reference {tgt:.3f} ± {tgt_sd:.3f}, a level-based number, so compare the three dz runs with each other first), total {r['ml_sig_tot']:.3f}, SGS share {100 * r['ml_sgs_share']:.1f} %")
    if res:
        fig, ax = plt.subplots(1, 2, figsize=(8.6, 4.2))
        for dz, r in res.items():
            zz = np.array(r["z"]) / r["zi"]
            ax[0].plot(r["sig_res"], zz, color=COL[dz], lw=2, label=f"dz = {dz} m (dx/dz = {100 // dz})")
            ax[1].plot(np.array(r["sig_tot"]), zz, color=COL[dz], lw=2)
        ax[0].axvline(tgt, color="#555555", lw=1, ls="--", label="fair-test reference: 25 m LES box-averaged to 100 m")
        ax[0].set_title("resolved $\\sigma_w$, dx = 100 m", fontsize=10); ax[1].set_title("resolved + SGS (isotropic share of $k_{sgs}$)", fontsize=10)
        for a in ax: a.set_ylim(0, 1.3); a.set_xlabel("$\\sigma_w$ (m/s)"); a.set_ylabel("$z/z_i$"); a.grid(color="#EEEEEE", lw=0.6)
        ax[0].legend(frameon=False, fontsize=8, loc="upper right")
        fig.suptitle("Grid aspect-ratio test: same dx = 100 m, dz = 20 / 50 / 100 m (3–4 h mean)", fontsize=11); fig.tight_layout()
        fig.savefig(f"{H}/Results Summary/figs/ar_test_sigw.png", dpi=150)
        json.dump({str(k): {kk: vv for kk, vv in v.items() if kk not in ("z", "sig_res", "sig_tot")} | {"target_les_box100": tgt, "target_sd": tgt_sd} for k, v in res.items()},
                  open(f"{H}/Results Summary/ar_test.json", "w"), indent=1)
        print("wrote Results Summary/figs/ar_test_sigw.png and ar_test.json")
