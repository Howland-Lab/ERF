#!/usr/bin/env python3
"""
week2_analysis.py -- analysis of the runs submitted on 15 Oct 2026 (make_week2_runs.py). Runs whatever is finished, skips the rest.

    python3 week2_analysis.py            # -> figs/w2_*.png, week2_summary.json, week2_summary.txt

A. Stable ensemble (cases/ens_sbl6_s1..5, plus the original cases/sbl_dx6p25 as member 0 where useful)
   A1 time series of u* and h for every member: how far the members drift apart (internal variability)
   A2 profiles at 9 h: ensemble mean and spread vs the 8-9 h time mean of one member
   A3 w at 0.5 h over 8-9 h: sigma, skewness, kurtosis and P(|w'| > 3 sigma) per member -> ensemble mean and the
      member-to-member spread = the confidence interval of a single-run estimate
B. Convective ensemble (cases_pvc/ens_cbl25_s1..4 plus cases/cbl_sd25 and cases_pvc/cbl_sd25): same as A at 3-4 h, 0.5 zi
C. Grey-zone scan (cases/gz_cbl_dx050..400, with cases/cbl_sd25 as the 25 m point on an 8 km domain)
   C1 mean profiles; C2 resolved share of TKE and of the heat flux at 0.5 zi against Delta/zi (Honnert-type partition)
   C3 w distribution at 0.5 zi (sigma, skewness, 3-sigma tail) against Delta/zi
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_wrf_erf as C
from erf_planes import Plotfile, Moments, HBINS

ROOT = os.path.join(HERE, "..")
FIG = os.path.join(HERE, "figs")
GUIDE = "/work2/11859/hsna/stampede3/HowlandLab/Results Summary/figs"
COL = ["#B2182B", "#EF8A62", "#2166AC", "#67A9CF", "#1B7837", "#762A83"]
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
OUT = {}


def done(path):
    log = os.path.join(path, "log.erf")
    return os.path.exists(log) and "finalized" in open(log, errors="ignore").read()[-3000:]


def profiles(path, dz, win):
    C.ERF = os.path.dirname(path)
    meta = dict(case="x", dz=dz, label="", win=win)
    mean, ts, n = C.erf_profiles(os.path.basename(path), meta)
    return mean, ts


def frames_at(path, times):
    pls = [Plotfile(os.path.join(path, p)) for p in sorted(os.listdir(path)) if p.startswith("plt")]
    return sorted([p for p in pls if any(abs(p.time - t) < 1 for t in times)], key=lambda p: p.time)


def w_stats(pls, k):
    m = Moments(); sig = None; ws = []
    for p in pls:
        w = p.plane("z_velocity", k); w = w - w.mean()
        if sig is None:
            sig = w.std()
        m.add(w, HBINS["w"], sig); ws.append(w)
    r = m.result()
    allw = np.concatenate([x.ravel() for x in ws]); s = allw.std()
    r["p3"] = float(np.mean(np.abs(allw) > 3 * s)); r["p3up"] = float(np.mean(allw > 3 * s))
    return r


def ensemble(kind):
    if kind == "sbl":
        members = [os.path.join(ROOT, "cases", f"ens_sbl6_s{s}") for s in range(1, 6)]
        dz, win, frac, tw, label = 6.25, (8.0, 9.0), 0.5, (28800, 32400), "stable 6.25 m"
    else:
        members = [os.path.join(ROOT, "cases_pvc", f"ens_cbl25_s{s}") for s in range(1, 5)] + \
                  [os.path.join(ROOT, "cases", "cbl_sd25"), os.path.join(ROOT, "cases_pvc", "cbl_sd25")]
        dz, win, frac, tw, label = 20.0, (3.0, 4.0), 0.5, (10800, 14400), "convective 25 m"
    members = [m for m in members if done(m)]
    if len(members) < 2:
        print(f"[{kind}] ensemble: {len(members)} finished members, skipped"); return
    res = []
    for m in members:
        mean, ts = profiles(m, dz, win)
        h = mean["h_stress"] if kind == "sbl" else mean["h_minflux"]
        zc = mean["zc"]; k = int(np.argmin(abs(zc - frac * h)))
        step = 120 if kind == "sbl" else 600
        pls = frames_at(m, np.arange(tw[0], tw[1] + 1, step))
        ws = w_stats(pls, k) if pls else None
        res.append(dict(name=os.path.relpath(m, ROOT), mean=mean, ts=ts, h=h, k=k, w=ws, nfr=len(pls)))
        print(f"[{kind}] {res[-1]['name']}: u*={mean['ust']:.3f} h={h:.0f} frames={len(pls)}" + (f" skew={ws['skew']:.2f} P3={ws['p3']:.1e}" if ws else ""))
    # A1/B1 time series
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    for i, r in enumerate(res):
        ax[0].plot(r["ts"]["ust_t"], r["ts"]["ust"], color=COL[i % 6], lw=0.8, label=r["name"].split("/")[-1])
        ax[1].plot(r["ts"]["t"], r["ts"]["h_stress"] if kind == "sbl" else r["ts"]["h_stress"], color=COL[i % 6], lw=0.8)
    ax[0].set_xlabel("time (h)"); ax[0].set_ylabel("u* (m/s)"); ax[1].set_xlabel("time (h)"); ax[1].set_ylabel("h, stress method (m)")
    for a in ax:
        a.axvspan(*win, color="#EEEEEE", zorder=0); a.grid(alpha=0.25)
    ax[0].legend(frameon=False, fontsize=7)
    fig.suptitle(f"ERF {label}: every ensemble member (grey band = averaging window)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/w2_ens_{kind}_timeseries.png", dpi=160); plt.close(fig)
    # A2/B2 profiles: member time means, ensemble mean and spread
    zc = res[0]["mean"]["zc"]; n = min(len(r["mean"]["zc"]) for r in res)
    spd = np.array([np.hypot(r["mean"]["u"][:n], r["mean"]["v"][:n]) for r in res])
    th = np.array([r["mean"]["th"][:n] for r in res])
    fig, ax = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
    for a, X, xl in ((ax[0], spd, "wind speed (m/s)"), (ax[1], th, r"$\theta$ (K)")):
        for i, x in enumerate(X):
            a.plot(x, zc[:n], color=COL[i % 6], lw=0.7, alpha=0.8)
        a.plot(X.mean(0), zc[:n], color="k", lw=1.8, label="ensemble mean")
        a.fill_betweenx(zc[:n], X.mean(0) - X.std(0, ddof=1), X.mean(0) + X.std(0, ddof=1), color="#999999", alpha=0.35, label="±1 std across members")
        a.set_xlabel(xl); a.grid(alpha=0.25)
    ax[0].set_ylim(0, 1.4 * np.mean([r["h"] for r in res])); ax[0].set_ylabel("z (m)"); ax[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"ERF {label}: time-mean profile of each member ({win[0]:.0f}–{win[1]:.0f} h)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/w2_ens_{kind}_profiles.png", dpi=160); plt.close(fig)
    # A3/B3 statistics with spread
    keys = ["ust", "h"]
    summ = {"members": [r["name"] for r in res]}
    for key in keys:
        v = np.array([r["mean"]["ust"] if key == "ust" else r["h"] for r in res])
        summ[key] = dict(mean=float(v.mean()), std=float(v.std(ddof=1)), cv=float(v.std(ddof=1) / v.mean()))
    W = [r["w"] for r in res if r["w"]]
    if W:
        for key in ("std", "skew", "kurt", "p3", "p3up"):
            v = np.array([w[key] for w in W])
            summ["w_" + key] = dict(mean=float(v.mean()), std=float(v.std(ddof=1)), values=v.tolist())
        fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
        c = 0.5 * (HBINS["w"][1:] + HBINS["w"][:-1]); dw = HBINS["w"][1] - HBINS["w"][0]
        for i, w in enumerate(W):
            ax[0].semilogy(c / w["std"], w["hist"] / w["n"] / dw * w["std"], color=COL[i % 6], lw=0.9)
        ax[0].semilogy(c, np.exp(-0.5 * c ** 2) / np.sqrt(2 * np.pi), "--", color="#999999", lw=0.8)
        ax[0].set_xlim(-6, 7); ax[0].set_ylim(1e-6, 1); ax[0].set_xlabel(r"$w'/\sigma_w$ at mid-layer"); ax[0].set_ylabel("PDF"); ax[0].set_title("each member", fontsize=10)
        for a, key, lab in ((ax[1], "skew", "skewness of w"), (ax[2], "p3", r"P($|w'|>3\sigma$)")):
            v = np.array([w[key] for w in W])
            a.plot(np.arange(len(v)), v, "o", color="k"); a.axhline(v.mean(), color="k", lw=1)
            a.axhspan(v.mean() - v.std(ddof=1), v.mean() + v.std(ddof=1), color="#DDDDDD", zorder=0)
            a.set_xlabel("member"); a.set_ylabel(lab); a.grid(alpha=0.25); a.set_title(f"spread: ±{v.std(ddof=1) / abs(v.mean()) * 100:.0f}% of the mean", fontsize=10)
        if kind == "sbl":
            ax[2].axhline(0.0027, color="#999999", ls="--", lw=0.8)
        fig.suptitle(f"ERF {label}: how much one run's statistics vary from member to member", fontsize=11)
        fig.tight_layout(); fig.savefig(f"{FIG}/w2_ens_{kind}_wstats.png", dpi=160); plt.close(fig)
    # non-stationarity (stable): ensemble mean at 8 h vs at 9 h vs member spread
    if kind == "sbl":
        def prof_at(path, t):
            C.ERF = os.path.dirname(path)
            Mf = C.frames(C.read_log(os.path.join(path, "mean_profiles.dat"), 16))
            tt = min(Mf, key=lambda x: abs(x - t)); b = Mf[tt]
            return np.hypot(b[:, 2], b[:, 3]), b[:, 6], b[:, 1] + 0.5 * dz
        P8 = [prof_at(os.path.join(ROOT, r["name"]), 28800) for r in res]
        P9 = [prof_at(os.path.join(ROOT, r["name"]), 32400) for r in res]
        z = P8[0][2]; m = min(len(p[1]) for p in P8 + P9)
        th8 = np.array([p[1][:m] for p in P8]); th9 = np.array([p[1][:m] for p in P9])
        k = int(np.argmin(abs(z[:m] - 0.5 * summ["h"]["mean"])))
        summ["nonstat"] = dict(dtheta_8to9_mid=float(th9[:, k].mean() - th8[:, k].mean()), spread_mid=float(th9[:, k].std(ddof=1)),
                               dtheta_8to9_low=float(th9[:, 1].mean() - th8[:, 1].mean()))
    OUT[f"ens_{kind}"] = summ


def greyzone():
    runs = [("cases/cbl_sd25", 25.0, 600), ("cases/gz_cbl_dx050", 50.0, 300), ("cases/gz_cbl_dx100", 100.0, 300),
            ("cases/gz_cbl_dx200", 200.0, 300), ("cases/gz_cbl_dx400", 400.0, 300)]
    res = []
    for rel, dx, step in runs:
        path = os.path.join(ROOT, rel)
        if not done(path):
            print(f"[gz] {rel} not finished, skipped"); continue
        mean, ts = profiles(path, 20.0, (3.0, 4.0))
        zi = mean["h_minflux"]; zc = mean["zc"]; k = int(np.argmin(abs(zc - 0.5 * zi)))
        res_tke = mean["tke_res"][k]; sgs_tke = mean["ksgs"][k]
        kf = int(np.argmin(abs(mean["zf"] - 0.5 * zi)))
        wth_tot = mean["wth"][kf]; wth_sgs = mean["wth_sgs"][kf]
        pls = frames_at(path, np.arange(10800, 14401, step))
        ws = w_stats(pls, k) if pls else None
        res.append(dict(rel=rel, dx=dx, zi=zi, ratio=dx / zi, mean=mean, tke_frac=res_tke / (res_tke + sgs_tke),
                        flux_frac=1 - wth_sgs / wth_tot if wth_tot else np.nan, w=ws, ust=mean["ust"]))
        print(f"[gz] dx={dx:.0f}: zi={zi:.0f} dx/zi={dx / zi:.3f} resolved TKE share={res[-1]['tke_frac']:.2f} resolved heat-flux share={res[-1]['flux_frac']:.2f}"
              + (f" w skew={ws['skew']:.2f} P3up={ws['p3up']:.1e}" if ws else ""))
    if len(res) < 2:
        return
    fig, ax = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for i, r in enumerate(res):
        m = r["mean"]; col = plt.cm.viridis(i / max(1, len(res) - 1))
        ax[0].plot(np.hypot(m["u"], m["v"]), m["zc"], color=col, label=f"Δx = {r['dx']:.0f} m")
        ax[1].plot(m["th"], m["zc"], color=col)
        ax[2].plot(m["tke_res"], m["zc"], color=col); ax[2].plot(m["ksgs"], m["zc"], ":", color=col, lw=0.9)
    ax[0].set_xlabel("wind speed (m/s)"); ax[1].set_xlabel(r"$\theta$ (K)"); ax[1].set_xlim(302, 312); ax[2].set_xlabel("TKE: resolved (solid), subgrid (dotted)")
    ax[0].set_ylabel("z (m)"); ax[0].legend(frameon=False, fontsize=8)
    for a in ax:
        a.set_ylim(0, 1600); a.grid(alpha=0.25)
    fig.suptitle("Grey-zone scan, convective case, mean over 3–4 h (same LES closure at every Δx)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/w2_gz_profiles.png", dpi=160); plt.close(fig)
    R = np.array([r["ratio"] for r in res])
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    ax[0].plot(R, [r["tke_frac"] for r in res], "o-", color="k", label="TKE")
    ax[0].plot(R, [r["flux_frac"] for r in res], "s--", color="#B2182B", label="heat flux")
    ax[0].set_ylabel("resolved share at 0.5 $z_i$"); ax[0].legend(frameon=False, fontsize=8); ax[0].set_ylim(0, 1.05)
    W = [(r["ratio"], r["w"]) for r in res if r["w"]]
    ax[1].plot([a for a, _ in W], [w["skew"] for _, w in W], "o-", color="k"); ax[1].set_ylabel("skewness of w at 0.5 $z_i$")
    ax[2].semilogy([a for a, _ in W], [max(w["p3up"], 1e-6) for _, w in W], "o-", color="k"); ax[2].set_ylabel(r"P($w' > 3\sigma_w$) at 0.5 $z_i$")
    ax[2].axhline(0.00135, color="#999999", ls="--", lw=0.8)
    for a in ax:
        a.set_xscale("log"); a.set_xlabel(r"$\Delta x / z_i$"); a.grid(alpha=0.25, which="both")
        a.axvspan(0.1, 1.0, color="#F3E9D2", zorder=0)
    fig.suptitle("What the grid resolves as Δx approaches the boundary-layer depth (shaded: grey zone, Δx/zi ≈ 0.1–1)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/w2_gz_partition.png", dpi=160); plt.close(fig)
    OUT["greyzone"] = [dict(dx=r["dx"], zi=r["zi"], ratio=r["ratio"], ust=r["ust"], tke_frac=r["tke_frac"], flux_frac=r["flux_frac"],
                            w_std=r["w"]["std"] if r["w"] else None, w_skew=r["w"]["skew"] if r["w"] else None,
                            w_p3up=r["w"]["p3up"] if r["w"] else None) for r in res]


def clean(o):
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, np.generic): return o.item()
    return o


if __name__ == "__main__":
    for f in (lambda: ensemble("sbl"), lambda: ensemble("cbl"), greyzone):
        try:
            f()
        except Exception as e:
            import traceback; traceback.print_exc()
    json.dump(clean(OUT), open(os.path.join(HERE, "week2_summary.json"), "w"), indent=1)
    os.makedirs(GUIDE, exist_ok=True)
    for fn in os.listdir(FIG):
        if fn.startswith("w2_"):
            os.system(f"cp '{FIG}/{fn}' '{GUIDE}/'")
    print("written week2_summary.json")
