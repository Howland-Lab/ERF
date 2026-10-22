#!/usr/bin/env python3
"""
gz_scheme_analysis.py -- the "first tests at the same Delta": PBL schemes and VLES closures on grey-zone grids
(convective case, 16 km domain, dz = 20 m, 3-4 h) against the fair-test reference = 25 m LES box-averaged to Delta
(gz_fairtest.json, 6 realisations).  Runs whatever has finished, skips the rest.

WRF (abl_cases/gz_{scheme}_dx{100,200,400}):  mynn_edmf, mynn_nomf, ysu, shinhong, 3dtke, deardorff
ERF (ERF_runs/cases/gzpbl_{scheme}_dx{...}):   mynn25, mynnedmf, mrf, mrfblend ; plus the VLES scan gz_cbl_dx{...}

For every run: plane-mean profiles (u, theta, resolved fluxes, resolved + modelled TKE) and, at the reference level
0.5 zi (zi = 1040 m from the LES), the resolved w statistics (sigma, skewness, P(w'>3 sigma)) from the 2-min planes.
"Total" sigma_w^2 = resolved + (2/3) k_model is the isotropic estimate a 1-D scheme allows (see the results note).
Outputs: figs/gz_schemes_fairtest.png, figs/gz_schemes_profiles.png, gz_schemes.json
"""
import glob, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from erf_planes import Plotfile, Moments, HBINS
from compare_wrf_erf import erf_profiles
WRFROOT = "/work2/11859/hsna/stampede3/HowlandLab/WRF-4.6.1/test/abl_cases"
sys.path.insert(0, os.path.join(WRFROOT, "postproc"))
import postproc_highres as PH
ERFC = os.path.join(HERE, "..", "cases"); FIG = os.path.join(HERE, "figs")
ZI = 1040.0; WIN = (10800.0, 14400.0); K = int(round(0.5 * ZI / 20.0 - 0.5)); ZREF = (K + 0.5) * 20.0
DXS = [100, 200, 400]
WRF_SCHEMES = {"mynn_edmf": "MYNN-EDMF (WRF)", "mynn_nomf": "MYNN, no mass flux (WRF)", "ysu": "YSU (WRF)",
               "shinhong": "Shin-Hong scale-aware (WRF)", "3dtke": "SMS-3DTKE (WRF)", "deardorff": "Deardorff VLES (WRF)"}
ERF_SCHEMES = {"mynn25": "MYNN2.5 (ERF)", "mynnedmf": "MYNN-EDMF (ERF)", "mrf": "MRF (ERF)", "mrfblend": "MRF + Boutle blend (ERF)",
               "vles": "Deardorff VLES (ERF)"}
COLS = {"mynn_edmf": "#2166AC", "mynn_nomf": "#6BAED6", "ysu": "#762A83", "shinhong": "#1B7837", "3dtke": "#E6AB02", "deardorff": "#222222",
        "mynn25": "#D9822B", "mynnedmf": "#B2182B", "mrf": "#9E9E9E", "mrfblend": "#4D4D4D", "vles": "#9B2335"}


def wstats(samples):
    x = np.concatenate([s.ravel() for s in samples]); x = x - x.mean(); s = x.std()
    return dict(std=float(s), skew=float(np.mean(x ** 3) / s ** 3), p3up=float(np.mean(x > 3 * s)), n=int(x.size))


def level_mean(d, var, it, nz):
    a = np.asarray(d[var][it]); m = a.mean(axis=(1, 2))
    return 0.5 * (m[:-1] + m[1:]) if len(m) == nz + 1 else m[:nz]


def wrf_run(scheme, dx):
    name = f"gz_{scheme}_dx{dx}"; path = os.path.join(WRFROOT, name)
    fs = sorted(glob.glob(f"{path}/wrfout_d01_*"))
    if not fs:
        return None
    d = nc.Dataset(fs[0]); d.set_auto_mask(False)
    t = np.asarray(d["XTIME"][:]) * 60.0
    sel = [i for i in range(len(t)) if WIN[0] - 1 <= t[i] <= WIN[1] + 1]
    if not sel or t[-1] < WIN[1] - 1:
        print(f"  {name}: incomplete ({t[-1] / 3600:.2f} h)"); return None
    nz = len(d.dimensions["bottom_top"])
    fr = []
    for it in sel:
        f = PH.frame_profiles(d, it, nz)
        f["tke_model"] = (level_mean(d, "TKE_PBL", it, nz) if "TKE_PBL" in d.variables else np.zeros(nz)) + (f["tke_sgs"] if scheme in ("3dtke", "deardorff") else 0.0)
        fr.append(f)
    m = {k: np.mean([f[k] for f in fr], axis=0) for k in fr[0] if np.ndim(fr[0][k]) > 0}
    for k in ("ust", "wth0"):
        m[k] = float(np.mean([f[k] for f in fr]))
    m["h_grad"] = float(PH.h_gradient(m["th"], m["zh"]))
    # resolved w at the reference level from the 2-min planes
    pf = glob.glob(f"{path}/planes_d01_*"); W = []
    if pf:
        p = nc.Dataset(pf[0]); p.set_auto_mask(False); tp = np.asarray(p["XTIME"][:]) * 60.0
        kk = int(np.argmin(abs(m["zh"] - ZREF)))
        for it in range(len(tp)):
            if WIN[0] - 1 <= tp[it] <= WIN[1] + 1:
                w = 0.5 * (np.asarray(p["W"][it, kk]) + np.asarray(p["W"][it, kk + 1])); W.append(w - w.mean())
    ws = wstats(W) if W else None
    kref = float(np.interp(ZREF, m["zh"], m["tke_model"]))
    return dict(name=name, label=WRF_SCHEMES[scheme], dx=dx, scheme=scheme, code="WRF", zc=m["zh"], zf=m["zf"],
                u=m["u"], v=m["v"], th=m["th"], tau=np.hypot(m["uw"], m["vw"]), wth=m["wth"], tke_res=m["tke_res"], tke_model=m["tke_model"],
                ust=m["ust"], wth0=m["wth0"], h_grad=m["h_grad"], w=ws, k_ref=kref,
                sig_tot=float(np.sqrt(ws["std"] ** 2 + 2.0 / 3.0 * kref)) if ws else None, sig_iso=float(np.sqrt(2.0 / 3.0 * kref)))


def erf_run(scheme, dx):
    name = f"gzpbl_{scheme}_dx{dx}" if scheme != "vles" else f"gz_cbl_dx{dx:03d}"; path = os.path.join(ERFC, name)
    log = os.path.join(path, "log.erf")
    if not (os.path.exists(log) and "Total Time" in open(log, errors="ignore").read()[-4000:]):
        return None
    E, _, nE = erf_profiles(name, dict(case="cbl", dz=20.0, win=(3.0, 4.0)))
    pls = [Plotfile(os.path.join(path, q)) for q in sorted(os.listdir(path)) if q.startswith("plt")]
    W = [q.plane("z_velocity", K) for q in pls if WIN[0] - 1 <= q.time <= WIN[1] + 1]
    W = [w - w.mean() for w in W]
    ws = wstats(W) if W else None
    kmod = E["ksgs"]; kref = float(np.interp(ZREF, E["zc"], kmod))
    return dict(name=name, label=ERF_SCHEMES[scheme], dx=dx, scheme=scheme, code="ERF", zc=E["zc"], zf=E["zf"],
                u=E["u"], v=E["v"], th=E["th"], tau=np.hypot(E["uw"], E["vw"]), wth=E["wth"], tke_res=E["tke_res"], tke_model=kmod,
                ust=E["ust"], wth0=E["wth0"], h_grad=float(PH.h_gradient(E["th"], E["zc"])), w=ws, k_ref=kref,
                sig_tot=float(np.sqrt(ws["std"] ** 2 + 2.0 / 3.0 * kref)) if ws else None, sig_iso=float(np.sqrt(2.0 / 3.0 * kref)))


def main():
    T = json.load(open(os.path.join(HERE, "gz_fairtest.json")))
    runs = []
    for dx in DXS:
        for s in WRF_SCHEMES:
            r = wrf_run(s, dx); runs.append(r) if r else None
        for s in ERF_SCHEMES:
            r = erf_run(s, dx); runs.append(r) if r else None
    print("finished runs:", [r["name"] for r in runs])
    if not runs:
        return
    # ---- figure 1: profiles per dx, all schemes
    fig, ax = plt.subplots(len(DXS), 4, figsize=(14, 3.6 * len(DXS)), sharey=True, squeeze=False)
    for i, dx in enumerate(DXS):
        rs = [r for r in runs if r["dx"] == dx]
        for r in rs:
            ls = "-" if r["code"] == "ERF" else "--"; c = COLS[r["scheme"]]
            ax[i, 0].plot(np.hypot(r["u"], r["v"]), r["zc"], ls, color=c, lw=1.4, label=r["label"])
            ax[i, 1].plot(r["th"], r["zc"], ls, color=c, lw=1.4)
            ax[i, 2].plot(r["tau"] / r["ust"] ** 2, r["zf"], ls, color=c, lw=1.4)
            ax[i, 3].plot(r["tke_res"] + r["tke_model"], r["zc"], ls, color=c, lw=1.4)
            ax[i, 3].plot(r["tke_model"], r["zc"], ":", color=c, lw=0.9)
        for a, xl in zip(ax[i], ["wind speed (m/s)", r"$\theta$ (K)", r"resolved $|\tau|/u_*^2$", "TKE: total (solid), modelled (dotted)"]):
            a.set_xlabel(xl); a.set_ylim(0, 1600); a.grid(alpha=0.25); a.axhline(ZI, color="#BBBBBB", lw=0.8)
        ax[i, 0].set_ylabel(f"dx = {dx} m\nz (m)"); ax[i, 1].set_xlim(302, 312)
        if rs:
            ax[i, 0].legend(frameon=False, fontsize=7)
    fig.suptitle("Grey-zone closure scan, convective case, mean over 3-4 h (dashed = WRF, solid = ERF)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/gz_schemes_profiles.png", dpi=140); plt.close(fig)
    # ---- figure 2: fair test at 0.5 zi
    Ls = [25, 50, 100, 200, 400]; x = np.array(Ls) / ZI
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.0))
    for a, key, lab, logy in ((ax[0], "std", r"$\sigma_w$ at 0.5 $z_i$ (m/s)", False), (ax[1], "skew", "skewness of resolved w", False), (ax[2], "p3", r"P($w'>3\sigma_w$), resolved", True)):
        tm = np.array([T[str(L)]["target_mean"][key] for L in Ls]); ts = np.array([T[str(L)]["target_std"][key] for L in Ls])
        a.fill_between(x, tm - ts, tm + ts, color="#CCCCCC", alpha=0.6, label="25 m LES averaged over Δ-boxes (±1 std, 6 runs)")
        a.plot(x, tm, "-", color="#555555")
        a.axvspan(0.1, 0.5, color="#F3E9D2", zorder=0); a.set_xscale("log"); a.set_xlabel(r"$\Delta/z_i$"); a.set_ylabel(lab); a.grid(alpha=0.25, which="both")
        if logy:
            a.set_yscale("log")
    for sch in list(WRF_SCHEMES) + list(ERF_SCHEMES):
        rs = sorted([r for r in runs if r["scheme"] == sch and r["w"]], key=lambda r: r["dx"])
        if not rs:
            continue
        xx = np.array([r["dx"] for r in rs]) / ZI; c = COLS[sch]; mk = "o" if rs[0]["code"] == "ERF" else "s"
        ax[0].plot(xx, [r["w"]["std"] for r in rs], mk + "-", color=c, lw=1.2, ms=5, label=rs[0]["label"] + " (resolved)")
        ax[0].plot(xx, [r["sig_tot"] for r in rs], mk + ":", color=c, lw=1.0, ms=4, mfc="white")
        ax[1].plot(xx, [r["w"]["skew"] for r in rs], mk + "-", color=c, lw=1.2, ms=5)
        ax[2].plot(xx, [max(r["w"]["p3up"], 1e-6) for r in rs], mk + "-", color=c, lw=1.2, ms=5)
    ax[0].legend(frameon=False, fontsize=6.5, loc="lower left")
    ax[0].text(0.02, 0.98, "dotted/open = resolved + (2/3) modelled TKE", transform=ax[0].transAxes, fontsize=7, va="top")
    fig.suptitle("Fair test at 0.5 zi: PBL schemes and VLES closures vs the box-averaged 25 m LES (convective, 3-4 h)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/gz_schemes_fairtest.png", dpi=150); plt.close(fig)
    out = {r["name"]: dict(label=r["label"], dx=r["dx"], code=r["code"], ust=r["ust"], h_grad=r["h_grad"], k_ref=r["k_ref"], sig_iso=r["sig_iso"],
                           sig_tot=r["sig_tot"], w=r["w"], ml_speed=float(np.mean(np.hypot(r["u"], r["v"])[(r["zc"] > 0.2 * ZI) & (r["zc"] < 0.6 * ZI)])))
           for r in runs}
    out["_target"] = {L: T[str(L)]["target_mean"] for L in Ls}
    json.dump(out, open(os.path.join(HERE, "gz_schemes.json"), "w"), indent=1)
    for n, v in out.items():
        if n.startswith("_"):
            continue
        w = v["w"] or {}
        print(f"{n:24s} u*={v['ust']:.3f} zi(grad)={v['h_grad']:6.0f} U_ml={v['ml_speed']:.2f} | sig_res={w.get('std', float('nan')):.3f} "
              f"sig_tot={v['sig_tot'] if v['sig_tot'] else float('nan'):.3f} skew={w.get('skew', float('nan')):+.2f} p3={w.get('p3up', float('nan')):.1e} | target sig={T[str(v['dx'])]['target_mean']['std']:.3f}")


if __name__ == "__main__":
    main()
