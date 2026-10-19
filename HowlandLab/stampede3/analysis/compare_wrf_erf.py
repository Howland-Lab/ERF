#!/usr/bin/env python3
"""
compare_wrf_erf.py -- first-order comparison of the ERF CPU runs with the matching WRF-LES runs.

    python3 compare_wrf_erf.py            # -> figs/*.png, summary.txt, summary.json

Pairs (same grid, same closure):  sbl_dx6p25 (6.25 m, Smagorinsky), sbl_dx3p125 (3.125 m, Smagorinsky),
cbl_sd25 (25 m x 20 m, 1.5-order TKE).  The 1-km MYNN RANS runs of the first report are drawn for context.

ERF profile logs (erf.destag_profiles = false -> ERF_Write1DProfiles_stag.cpp):
  mean_profiles.dat   t z u v w rho theta ksgs Kmv Khv qv qc qr qi qs qg
                      cell-centred values written against the LOWER FACE height -> plotted at z + dz/2
  covar_profiles.dat  t z u'u' u'v' u'w' v'v' v'w' w'w' u'th' v'th' w'th' th'th' ...   (fluxes on faces z;
                      u'u', v'v' are cell-centred values written at the lower face)
  sfs_profiles.dat    t z tau11 tau12 tau13 tau22 tau23 tau33 sgshfx sgsq1 sgsq2 sgsdiss  (faces z)
  surf_hist.dat       t u_star t_star L
Checked 13 Oct 2026: tau13/tau23 and sgshfx are kinematic (no rho) and have the same sign as u'w', w'th'
(at the wall |tau| = u*^2 within 6 %, sgshfx = -u* t* exactly).
WRF numbers come from postproc/data_hr/<run>.npz (postproc_highres.py): last-hour plane means, total
fluxes = resolved + SGS at cell centres, h_stress = height where |tau| = 5 % of its lowest value, / 0.95.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ERF = os.path.join(HERE, "..", "cases")
WRFP = "/work2/11859/hsna/stampede3/HowlandLab/WRF-4.6.1/test/abl_cases/postproc"
FIG = os.path.join(HERE, "figs"); os.makedirs(FIG, exist_ok=True)

PAIRS = {
    "sbl_dx6p25":  dict(case="sbl", dz=6.25,  label="SBL 6.25 m",  win=(8.0, 9.0)),
    "sbl_dx3p125": dict(case="sbl", dz=3.125, label="SBL 3.125 m", win=(8.0, 9.0)),
    "cbl_sd25":    dict(case="cbl", dz=20.0,  label="CBL 25 m",    win=(3.0, 4.0)),
}
C_WRF, C_ERF, C_RANS = "#222222", "#9B2335", "#999999"


def read_log(path, ncol):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) != ncol:
                continue
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                pass
    return np.array(rows)


def frames(a):
    """split a log array into {time: block}, blocks sorted by height"""
    out = {}
    for t in np.unique(a[:, 0]):
        b = a[a[:, 0] == t]
        out[t] = b[np.argsort(b[:, 1])]
    return out


def h_stress(tau, z):                                   # identical to postproc_highres.py
    thr = 0.05 * tau[0]
    for k in range(1, len(z)):
        if tau[k] < thr:
            f = (tau[k - 1] - thr) / max(tau[k - 1] - tau[k], 1e-12)
            return float((z[k - 1] + f * (z[k] - z[k - 1])) / 0.95)
    return np.nan


def erf_profiles(name, meta):
    d = os.path.join(ERF, name)
    M = frames(read_log(f"{d}/mean_profiles.dat", 16))
    C = frames(read_log(f"{d}/covar_profiles.dat", 22))
    S = frames(read_log(f"{d}/sfs_profiles.dat", 12))
    sh = read_log(f"{d}/surf_hist.dat", 4)
    t0, t1 = (3600 * w for w in meta["win"])
    times = sorted(t for t in M if t0 - 1 <= t <= t1 + 1 and t in C and t in S)

    def one(t):
        m, c, s = M[t], C[t], S[t]
        nz = min(len(m), len(c), len(s)) - 1               # drop the top-face row
        zf = c[:nz, 1]; zc = m[:nz, 1] + 0.5 * meta["dz"]
        uw = c[:nz, 4] + s[:nz, 4]; vw = c[:nz, 6] + s[:nz, 6]; wth = c[:nz, 10] + s[:nz, 8]
        r = dict(zc=zc, zf=zf, u=m[:nz, 2], v=m[:nz, 3], th=m[:nz, 6], ksgs=m[:nz, 7],
                 uw=uw, vw=vw, wth=wth, uw_sgs=s[:nz, 4], vw_sgs=s[:nz, 6], wth_sgs=s[:nz, 8],
                 tke_res=0.5 * (c[:nz, 2] + c[:nz, 5] + c[:nz, 7]), w_var=c[:nz, 7])
        k = slice(1, None)                                 # faces above the wall, as WRF starts at its first level
        r["h_stress"] = h_stress(np.hypot(uw[k], vw[k]), zf[k])
        r["h_minflux"] = float(zf[1 + int(np.argmin(wth[1:]))])
        return r

    fr = [one(t) for t in times]
    mean = {k: np.mean([r[k] for r in fr], axis=0) for k in fr[0]}
    mean["h_stress"] = h_stress(np.hypot(mean["uw"][1:], mean["vw"][1:]), mean["zf"][1:])
    mean["h_minflux"] = float(mean["zf"][1 + int(np.argmin(mean["wth"][1:]))])
    w = (sh[:, 0] >= t0 - 1) & (sh[:, 0] <= t1 + 1)
    mean["ust"] = float(np.mean(sh[w, 1])); mean["tst"] = float(np.mean(sh[w, 2])); mean["L"] = float(np.mean(sh[w, 3]))
    mean["wth0"] = -mean["ust"] * mean["tst"]
    allt = sorted(t for t in M if t in C and t in S and t > 0)
    ts = dict(t=np.array(allt) / 3600.0, h_stress=np.array([one(t)["h_stress"] for t in allt]))
    ts["ust_t"] = sh[:, 0] / 3600.0; ts["ust"] = sh[:, 1]
    return mean, ts, len(times)


def wrf_profiles(name):
    d = np.load(f"{WRFP}/data_hr/{name}.npz", allow_pickle=True)
    m = d["mean"].item(); ts = d["ts"].item()
    return m, ts


def rans(case):
    d = np.load(f"{WRFP}/data/{case}_rans.npz", allow_pickle=True)
    g = lambda k: d["mean_" + k]
    return dict(zc=g("zh"), u=g("u"), v=g("v"), th=g("th"), uw=g("uw"), vw=g("vw"), wth=g("wth"),
                zflux=g("zflux"), ust=float(g("ust")), h_stress=float(g("h_stress")))


def llj(zc, spd, top):
    k = zc < top
    i = int(np.argmax(spd[k])); return float(spd[k][i]), float(zc[k][i])


def figure(name, meta, E, W, R):
    case = meta["case"]
    hW = W["h_stress"] if case == "sbl" else W["h_minflux"]
    ztop = 1.4 * hW
    fig, ax = plt.subplots(1, 5, figsize=(15, 4.2), sharey=True)
    spdE, spdW = np.hypot(E["u"], E["v"]), np.hypot(W["u"], W["v"])
    ax[0].plot(spdW, W["zh"], color=C_WRF, label="WRF-LES"); ax[0].plot(spdE, E["zc"], color=C_ERF, label="ERF-LES")
    if R is not None:
        ax[0].plot(np.hypot(R["u"], R["v"]), R["zc"], "--", color=C_RANS, label="WRF 1-km MYNN")
    ax[0].set_xlabel("wind speed (m/s)"); ax[0].legend(frameon=False, fontsize=8)
    ax[1].plot(W["th"], W["zh"], color=C_WRF); ax[1].plot(E["th"], E["zc"], color=C_ERF)
    if R is not None:
        ax[1].plot(R["th"], R["zc"], "--", color=C_RANS)
    ax[1].set_xlabel(r"$\theta$ (K)")
    ax[2].plot(np.hypot(W["uw"], W["vw"]) / W["ust"] ** 2, W["zh"], color=C_WRF)
    ax[2].plot(np.hypot(E["uw"], E["vw"]) / E["ust"] ** 2, E["zf"], color=C_ERF)
    ax[2].plot(np.hypot(E["uw_sgs"], E["vw_sgs"]) / E["ust"] ** 2, E["zf"], ":", color=C_ERF, lw=1, label="ERF SGS part")
    ax[2].set_xlabel(r"$|\tau|/u_*^2$"); ax[2].legend(frameon=False, fontsize=8)
    ax[3].plot(W["wth"] / abs(W["wth0"]), W["zh"], color=C_WRF); ax[3].plot(E["wth"] / abs(E["wth0"]), E["zf"], color=C_ERF)
    ax[3].plot(E["wth_sgs"] / abs(E["wth0"]), E["zf"], ":", color=C_ERF, lw=1)
    ax[3].set_xlabel(r"$\overline{w'\theta'}/|\overline{w'\theta'}_0|$")
    ax[4].plot(W["tke_res"], W["zh"], color=C_WRF, label="resolved"); ax[4].plot(E["tke_res"], E["zc"], color=C_ERF)
    if case == "cbl":
        ax[4].plot(W["tke_sgs"], W["zh"], ":", color=C_WRF, lw=1, label="SGS"); ax[4].plot(E["ksgs"], E["zc"], ":", color=C_ERF, lw=1)
        ax[4].legend(frameon=False, fontsize=8)
    ax[4].set_xlabel(r"TKE (m$^2$/s$^2$)")
    for a in ax:
        a.set_ylim(0, ztop); a.grid(alpha=0.25)
    ax[0].set_ylabel("z (m)")
    fig.suptitle(f"{meta['label']}: ERF vs WRF, mean over {meta['win'][0]:.0f}-{meta['win'][1]:.0f} h", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/{name}_profiles.png", dpi=150); plt.close(fig)


def timeseries(results):
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    for j, (name, case) in enumerate([("sbl_dx6p25", "sbl"), ("sbl_dx3p125", "sbl"), ("cbl_sd25", "cbl")]):
        E, W, Ets, Wts = results[name]
        a = ax[j]
        a.plot(Wts["t"], Wts["ust"], color=C_WRF, label="WRF u*"); a.plot(Ets["ust_t"], Ets["ust"], color=C_ERF, lw=0.8, label="ERF u*")
        a.set_xlabel("time (h)"); a.set_ylabel("u* (m/s)"); a.set_title(PAIRS[name]["label"], fontsize=10); a.grid(alpha=0.25)
        if case == "sbl":
            b = a.twinx(); b.plot(Wts["t"], Wts["h_stress"], "--", color=C_WRF); b.plot(Ets["t"], Ets["h_stress"], "--", color=C_ERF, lw=0.8)
            b.set_ylabel("h (m), dashed"); b.set_ylim(0, 300)
        if j == 0:
            a.legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout(); fig.savefig(f"{FIG}/timeseries_ust_h.png", dpi=150); plt.close(fig)


def convergence(results):
    fig, ax = plt.subplots(1, 3, figsize=(11, 4), sharey=True)
    for name, ls in (("sbl_dx6p25", "--"), ("sbl_dx3p125", "-")):
        E, W, _, _ = results[name]; lab = PAIRS[name]["dz"]
        ax[0].plot(np.hypot(W["u"], W["v"]), W["zh"], ls, color=C_WRF, label=f"WRF {lab} m")
        ax[0].plot(np.hypot(E["u"], E["v"]), E["zc"], ls, color=C_ERF, label=f"ERF {lab} m")
        ax[1].plot(W["th"], W["zh"], ls, color=C_WRF); ax[1].plot(E["th"], E["zc"], ls, color=C_ERF)
        ax[2].plot(np.hypot(W["uw"], W["vw"]), W["zh"], ls, color=C_WRF); ax[2].plot(np.hypot(E["uw"], E["vw"]), E["zf"], ls, color=C_ERF)
    ax[0].set_xlabel("wind speed (m/s)"); ax[1].set_xlabel(r"$\theta$ (K)"); ax[2].set_xlabel(r"$|\tau|$ (m$^2$/s$^2$)")
    ax[0].set_ylabel("z (m)"); ax[0].legend(frameon=False, fontsize=8)
    for a in ax:
        a.set_ylim(0, 300); a.grid(alpha=0.25)
    fig.suptitle("SBL grid convergence: 6.25 m (dashed) vs 3.125 m (solid), mean over 8-9 h", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/sbl_convergence.png", dpi=150); plt.close(fig)


def main():
    results, table = {}, {}
    rr = {"sbl": rans("sbl"), "cbl": rans("cbl")}
    for name, meta in PAIRS.items():
        E, Ets, nE = erf_profiles(name, meta)
        W, Wts = wrf_profiles(name)
        results[name] = (E, W, Ets, Wts)
        figure(name, meta, E, W, rr[meta["case"]])
        row = {}
        for tag, P, zc in (("WRF", W, W["zh"]), ("ERF", E, E["zc"])):
            spd = np.hypot(P["u"], P["v"])
            d = dict(ust=P["ust"], wth0=P["wth0"], L=P["L"], h_stress=P["h_stress"], h_minflux=P["h_minflux"])
            if meta["case"] == "sbl":
                d["llj_speed"], d["llj_z"] = llj(zc, spd, 350.0)
                d["th_bottom"] = float(P["th"][0])
            else:
                k = (zc > 0.2 * P["h_minflux"]) & (zc < 0.8 * P["h_minflux"])
                d["ml_speed"] = float(np.mean(spd[k])); d["ml_theta"] = float(np.mean(P["th"][k]))
                d["tke_res_max"] = float(np.max(P["tke_res"]))
            row[tag] = d
        row["ERF_frames"] = nE
        table[name] = row
    timeseries(results); convergence(results)
    R = {c: dict(ust=rr[c]["ust"], h_stress=rr[c]["h_stress"]) for c in rr}
    with open(os.path.join(HERE, "summary.json"), "w") as f:
        json.dump(dict(pairs=table, rans=R), f, indent=1)
    lines = []
    for name, row in table.items():
        lines.append(f"== {PAIRS[name]['label']}  ({row['ERF_frames']} ERF frames in window)")
        keys = [k for k in row["WRF"]]
        for k in keys:
            w, e = row["WRF"][k], row["ERF"][k]
            rel = (e - w) / abs(w) * 100 if w not in (0,) and np.isfinite(w) and np.isfinite(e) else np.nan
            lines.append(f"   {k:11s} WRF {w:10.4g}   ERF {e:10.4g}   diff {rel:+6.1f} %")
    lines.append(f"== WRF 1-km MYNN: sbl u*={R['sbl']['ust']:.3f} h={R['sbl']['h_stress']:.0f} | cbl u*={R['cbl']['ust']:.3f}")
    txt = "\n".join(lines); print(txt)
    open(os.path.join(HERE, "summary.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
