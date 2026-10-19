#!/usr/bin/env python3
"""
sd25_fix_check.py -- did the non-tiling rerun of WRF cbl_sd25 (14x15 ranks) remove the 600-900 m anomaly?

The rerun lost its wrfout (pnetcdf could not create the file at t = 0; the run itself finished), so everything
here is computed from the last-hour 3-D plane output (U, V, W, T every 2 min, 31 frames) and compared on a
RESOLVED-only basis:
  WRF sd25_fix : planes_d01 (31 frames, 3-4 h)
  WRF sd25     : npz from wrfout (7 frames, 3-4 h; resolved parts uw_res, wth_res, tke_res)  -- the tiled run
  WRF dx10     : npz from wrfout (valid, 10 m)
  ERF sd25     : covar_profiles (resolved covariances), mean over 3-4 h
Heights: WRF levels from the sd25 npz (identical eta list); T in planes is theta - 300 (WRF perturbation theta).
Also w moments on the same planes as postproc_highres.py (z/zi = 0.1, 0.5, 0.9 and 100 m).
"""
import glob, json, os, sys
import numpy as np, netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from erf_planes import Moments, HBINS
from compare_wrf_erf import erf_profiles, PAIRS

WRF = "/work2/11859/hsna/stampede3/HowlandLab/WRF-4.6.1/test/abl_cases"
WRFP = WRF + "/postproc"
FIG = os.path.join(HERE, "figs")


def planes_profiles(run, zh):
    d = nc.Dataset(sorted(glob.glob(f"{WRF}/{run}/planes_d01_*"))[0]); d.set_auto_mask(False)
    nt = len(d.dimensions["Time"]); nz = len(d.dimensions["bottom_top"])
    acc = {k: np.zeros(nz) for k in ("u", "v", "th", "uw", "vw", "wth", "tke", "wvar")}
    for it in range(nt):
        for k in range(nz):
            U = np.asarray(d["U"][it, k], float); u = 0.5 * (U[:, 1:] + U[:, :-1])
            V = np.asarray(d["V"][it, k], float); v = 0.5 * (V[1:, :] + V[:-1, :])
            w = 0.5 * (np.asarray(d["W"][it, k], float) + np.asarray(d["W"][it, k + 1], float))
            th = np.asarray(d["T"][it, k], float) + 300.0
            up, vp, wp, tp = u - u.mean(), v - v.mean(), w - w.mean(), th - th.mean()
            acc["u"][k] += u.mean(); acc["v"][k] += v.mean(); acc["th"][k] += th.mean()
            acc["uw"][k] += np.mean(up * wp); acc["vw"][k] += np.mean(vp * wp); acc["wth"][k] += np.mean(wp * tp)
            acc["tke"][k] += 0.5 * np.mean(up ** 2 + vp ** 2 + wp ** 2); acc["wvar"][k] += np.mean(wp ** 2)
        print(f"  {run} frame {it + 1}/{nt}", flush=True)
    out = {k: v / nt for k, v in acc.items()}; out["zh"] = zh[:nz]
    out["zi"] = float(out["zh"][2 + int(np.argmin(out["wth"][2:]))])
    # w moments at the standard planes
    zlist = [("0.10", 0.1 * out["zi"]), ("0.50", 0.5 * out["zi"]), ("0.90", 0.9 * out["zi"]), ("hub", 100.0)]
    mom = {}
    for lab, z in zlist:
        k = int(np.argmin(abs(out["zh"] - z))); m = Moments(); sig = None
        for it in range(nt):
            w = 0.5 * (np.asarray(d["W"][it, k], float) + np.asarray(d["W"][it, k + 1], float)); w -= w.mean()
            if sig is None:
                sig = w.std()
            m.add(w, HBINS["w"], sig)
        mom[lab] = m.result()
    out["mom"] = mom
    return out


def main():
    old = np.load(f"{WRFP}/data_hr/cbl_sd25.npz", allow_pickle=True)
    Mo = old["mean"].item(); Po = old["planes"].item()
    ten = np.load(f"{WRFP}/data_hr/cbl_dx10.npz", allow_pickle=True)
    Mt = ten["mean"].item(); Pt = ten["planes"].item()
    F = planes_profiles("cbl_sd25_fix", Mo["zh"])
    E, _, _ = erf_profiles("cbl_sd25", PAIRS["cbl_sd25"])
    ustE = E["ust"]
    rows = [("WRF 25 m, 14×16 ranks (tiled)", Mo["zh"], np.hypot(Mo["u"], Mo["v"]), Mo["th"], np.hypot(Mo["uw_res"], Mo["vw_res"]), Mo["wth_res"], Mo["tke_res"], "#BBBBBB", ":"),
            ("WRF 25 m, 14×15 ranks (rerun)", F["zh"], np.hypot(F["u"], F["v"]), F["th"], np.hypot(F["uw"], F["vw"]), F["wth"], F["tke"], "#222222", "-"),
            ("WRF 10 m", Mt["zh"], np.hypot(Mt["u"], Mt["v"]), Mt["th"], np.hypot(Mt["uw_res"], Mt["vw_res"]), Mt["wth_res"], Mt["tke_res"], "#777777", "-"),
            ("ERF 25 m", None, None, None, None, None, None, "#9B2335", "-")]
    fig, ax = plt.subplots(1, 4, figsize=(13, 4.2), sharey=True)
    for lab, z, spd, th, tau, wth, tke, col, ls in rows:
        if lab.startswith("ERF"):
            ax[0].plot(np.hypot(E["u"], E["v"]), E["zc"], ls, color=col, label=lab, lw=1.6)
            ax[1].plot(E["th"], E["zc"], ls, color=col, lw=1.6)
            ax[2].plot(np.hypot(E["uw"] - E["uw_sgs"], E["vw"] - E["vw_sgs"]), E["zf"], ls, color=col, lw=1.6)
            ax[3].plot(E["tke_res"], E["zc"], ls, color=col, lw=1.6)
            continue
        ax[0].plot(spd, z, ls, color=col, label=lab, lw=1.2); ax[1].plot(th, z, ls, color=col, lw=1.2)
        ax[2].plot(tau, z, ls, color=col, lw=1.2); ax[3].plot(tke, z, ls, color=col, lw=1.2)
    ax[0].set_xlabel("wind speed (m/s)"); ax[1].set_xlabel(r"$\theta$ (K)"); ax[1].set_xlim(302, 312)
    ax[2].set_xlabel(r"resolved $|\overline{u'w'}|$ (m$^2$/s$^2$)"); ax[3].set_xlabel(r"resolved TKE (m$^2$/s$^2$)")
    ax[0].set_ylabel("z (m)"); ax[0].legend(frameon=False, fontsize=8, loc="upper left")
    for a in ax:
        a.set_ylim(0, 1500); a.grid(alpha=0.25)
    fig.suptitle("Convective case, mean over 3–4 h: the rerun without tiling vs ERF and the 10 m run", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/sd25_fix_profiles.png", dpi=150); plt.close(fig)

    lines = [f"WRF sd25_fix zi(min resolved heat flux) = {F['zi']:.0f} m"]
    for lab in ("0.10", "0.50", "0.90", "hub"):
        f = F["mom"][lab]; o = Po[f"{lab}_w"]; t = Pt[f"{lab}_w"]
        lines.append(f"  w at {lab}: std/skew/kurt  fix {f['std']:.3f} {f['skew']:+.2f} {f['kurt']:.2f} P>3s {f['tail'][2]:.1e} | "
                     f"tiled {o['std']:.3f} {o['skew']:+.2f} {o['kurt']:.2f} {o['tail'][2]:.1e} | dx10 {t['std']:.3f} {t['skew']:+.2f} {t['kurt']:.2f} {t['tail'][2]:.1e}")
    k = (F["zh"] > 200) & (F["zh"] < 800)
    lines.append(f"  mixed-layer (200-800 m) wind speed: fix {np.mean(np.hypot(F['u'], F['v'])[k]):.2f}, tiled {np.mean(np.hypot(Mo['u'], Mo['v'])[k]):.2f}, "
                 f"ERF {np.mean(np.hypot(E['u'], E['v'])[(E['zc'] > 200) & (E['zc'] < 800)]):.2f} m/s")
    lines.append(f"  speed range 200-800 m (max-min): fix {np.ptp(np.hypot(F['u'], F['v'])[k]):.2f}, tiled {np.ptp(np.hypot(Mo['u'], Mo['v'])[k]):.2f}, "
                 f"ERF {np.ptp(np.hypot(E['u'], E['v'])[(E['zc'] > 200) & (E['zc'] < 800)]):.2f} m/s")
    txt = "\n".join(lines); print(txt); open(os.path.join(HERE, "sd25_fix_summary.txt"), "w").write(txt + "\n")
    # PDFs
    fig, ax = plt.subplots(1, 4, figsize=(14, 3.4), sharey=True)
    c = 0.5 * (HBINS["w"][1:] + HBINS["w"][:-1]); dw = HBINS["w"][1] - HBINS["w"][0]
    for j, lab in enumerate(("0.10", "0.50", "0.90", "hub")):
        for q, col, ls, name in ((Po[f"{lab}_w"], "#BBBBBB", ":", "WRF 25 m tiled"), (Pt[f"{lab}_w"], "#777777", "-", "WRF 10 m"),
                                 (F["mom"][lab], "#222222", "-", "WRF 25 m rerun")):
            ax[j].semilogy(c / q["std"], q["hist"] / q["n"] / dw * q["std"], ls, color=col, lw=1.1, label=name)
        ax[j].set_xlim(-5, 7); ax[j].set_ylim(1e-6, 1); ax[j].grid(alpha=0.25); ax[j].set_xlabel(r"$w'/\sigma_w$")
        ax[j].set_title(f"z/zi = {lab}" if lab != "hub" else "z = 100 m", fontsize=10)
    ax[0].set_ylabel("PDF"); ax[0].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/sd25_fix_pdfs.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
