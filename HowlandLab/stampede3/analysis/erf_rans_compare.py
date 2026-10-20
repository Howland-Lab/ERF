#!/usr/bin/env python3
"""erf_rans_compare.py -- starting-guide item "run an NWP version in WRF/ERF": ERF 1 km MYNN2.5 (+ 2-D Smagorinsky
horizontally) against WRF 1 km MYNN2.5 (km_opt=4), with the LES as reference (WRF and ERF LES at the finest grid).
Window: SBL 8-9 h, CBL 3-4 h.   -> figs/rans1km_{sbl,cbl}.png, rans1km_summary.json"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from compare_wrf_erf import erf_profiles, wrf_profiles, WRFP, FIG, HERE

CASES = {"sbl": dict(erf="rans_sbl_1km", dz=400 / 64, win=(8.0, 9.0), les_w="sbl_dx3p125", les_e="sbl_dx3p125", les_dz=3.125, label="Stable (GABLS1)"),
         "cbl": dict(erf="rans_cbl_1km", dz=3500 / 120, win=(3.0, 4.0), les_w="cbl_dx10", les_e="cbl_sd25", les_dz=20.0, label="Convective (Shin & Dudhia)")}
C_LES, C_ERFLES, C_WRF1, C_ERF1 = "#222222", "#9B2335", "#2166AC", "#D9822B"
out = {}
for case, C in CASES.items():
    E1, _, nE = erf_profiles(C["erf"], dict(case=case, dz=C["dz"], win=C["win"]))
    d = np.load(f"{WRFP}/data/{case}_rans.npz", allow_pickle=True); g = lambda k: d["mean_" + k]
    W1 = dict(zc=g("zh"), u=g("u"), v=g("v"), th=g("th"), uw=g("uw"), vw=g("vw"), wth=g("wth"), zf=g("zflux"), k=g("tke_sgs"),
              ust=float(g("ust")), h_stress=float(g("h_stress")), h_minflux=float(g("h_minflux")), wth0=float(g("wth0")))
    WL, _ = wrf_profiles(C["les_w"])
    EL, _, _ = erf_profiles(C["les_e"], dict(case=case, dz=C["les_dz"], win=C["win"]))
    ztop = 1.5 * (WL["h_stress"] if case == "sbl" else WL["h_minflux"])
    fig, ax = plt.subplots(1, 5, figsize=(15, 4.3), sharey=True)
    rows = [("LES (WRF)", WL, WL["zh"], WL["zh"], C_LES, "-", 1.8), ("LES (ERF)", EL, EL["zc"], EL["zf"], C_ERFLES, "-", 1.0),
            ("1 km MYNN-EDMF, WRF (mass flux on, default)", W1, W1["zc"], W1["zf"], C_WRF1, "--", 1.6), ("1 km MYNN2.5, ERF", E1, E1["zc"], E1["zf"], C_ERF1, "-", 1.6)]
    W0 = None
    f0 = f"{WRFP}/data/{case}_rans_noedmf.npz"
    if os.path.exists(f0):                              # 20 Oct: WRF with bl_mynn_edmf=0 (mass flux off) = like-for-like with ERF MYNN25
        d0 = np.load(f0, allow_pickle=True); g0 = lambda k: d0["mean_" + k]
        W0 = dict(zc=g0("zh"), u=g0("u"), v=g0("v"), th=g0("th"), uw=g0("uw"), vw=g0("vw"), wth=g0("wth"), zf=g0("zflux"), k=g0("tke_sgs"),
                  ust=float(g0("ust")), h_stress=float(g0("h_stress")), h_minflux=float(g0("h_minflux")), wth0=float(g0("wth0")))
        rows.append(("1 km MYNN, WRF, mass flux off", W0, W0["zc"], W0["zf"], "#6BAED6", ":", 1.8))
    for lab, R, zc, zf, col, ls, lw in rows:
        ust = R["ust"]
        ax[0].plot(np.hypot(R["u"], R["v"]), zc, ls, color=col, lw=lw, label=lab)
        ax[1].plot(R["th"], zc, ls, color=col, lw=lw)
        if not (W0 is not None and R is W0):           # its flux diagnostic (K x gradient only) misses the flux in a well-mixed layer
            ax[2].plot(np.hypot(R["uw"], R["vw"]) / ust ** 2, zf, ls, color=col, lw=lw)
            w0 = abs(R["wth0"]) if abs(R["wth0"]) > 0 else 1.0
            ax[3].plot(R["wth"] / w0, zf, ls, color=col, lw=lw)
        if R is WL:   tke = WL["tke_res"] + WL["tke_sgs"]
        elif R is EL: tke = EL["tke_res"] + EL["ksgs"]
        elif R is W1: tke = W1["k"]
        elif W0 is not None and R is W0: tke = W0["k"]
        else:         tke = E1["ksgs"] + E1["tke_res"]
        ax[4].plot(tke, zc, ls, color=col, lw=lw)
    for a, xl in zip(ax, ["wind speed (m/s)", r"$\theta$ (K)", r"$|\tau|/u_*^2$", r"$\overline{w'\theta'}/|\overline{w'\theta'}_0|$", r"TKE (m$^2$/s$^2$)"]):
        a.set_xlabel(xl); a.set_ylim(0, ztop); a.grid(alpha=0.25)
    ax[0].set_ylabel("z (m)")
    h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=len(l), frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(f"{C['label']}: 1 km MYNN in ERF and in WRF vs the LES, mean over {C['win'][0]:.0f}-{C['win'][1]:.0f} h", fontsize=11)
    fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(f"{FIG}/rans1km_{case}.png", dpi=150); plt.close(fig)
    def summ(R, zc, k):
        S = dict(ust=float(R["ust"]), h_stress=float(R["h_stress"]), h_minflux=float(R["h_minflux"]))
        spd = np.hypot(R["u"], R["v"])
        if case == "sbl":
            j = int(np.argmax(np.where(zc < 300, spd, 0))); S.update(llj_speed=float(spd[j]), llj_z=float(zc[j]))
        else:
            zi = S["h_minflux"]; m = (zc > 0.2 * zi) & (zc < 0.6 * zi); S.update(ml_speed=float(spd[m].mean()), ml_theta=float(R["th"][m].mean()))
        if k is not None:
            for zz in ((50.0, 97.0) if case == "sbl" else (100.0, 544.0)):
                S[f"sigw_nwp_{zz:.0f}m"] = float(np.sqrt(2.0 / 3.0 * np.interp(zz, zc, k)))
        return S
    out[case] = {"LES_WRF": summ(WL, WL["zh"], None), "LES_ERF": summ(EL, EL["zc"], None),
                 "MYNN_WRF": summ(W1, W1["zc"], W1["k"]), "MYNN_ERF": summ(E1, E1["zc"], E1["ksgs"]), "ERF_frames": nE}
    if W0 is not None:
        out[case]["MYNN_WRF_noMF"] = summ(W0, W0["zc"], W0["k"])
json.dump(out, open(os.path.join(HERE, "rans1km_summary.json"), "w"), indent=1)
for case, v in out.items():
    print(case)
    for k2, s in v.items():
        if isinstance(s, dict):
            print("  %-9s " % k2 + "  ".join(f"{a}={b:.3f}" for a, b in s.items()))
