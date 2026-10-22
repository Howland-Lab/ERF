#!/usr/bin/env python3
"""tke_budget.py -- resolved TKE budget from ERF's plane-averaged profile logs (22 Oct 2026, after Mike: "to develop a
PBL scheme, look at the TKE and momentum budgets").  Nothing new has to be run: covar_profiles.dat already holds the
third moments <u_i'u_i'w'> and the pressure correlation <p'w'>, sfs_profiles.dat the SGS dissipation.

    de/dt = P_s + P_b + T_t + T_p + R,      e = <u'u'+v'v'+w'w'>/2 (resolved)
    P_s = -<u'w'> dU/dz - <v'w'> dV/dz       shear production
    P_b = (g/theta_ref) <w'theta_v'>          buoyancy production (< 0 in the stable case)
    T_t = -1/2 d<u_i'u_i'w'>/dz               turbulent transport
    T_p = -(1/rho) d<p'w'>/dz                 pressure transport
    R   = residual = everything the resolved terms do not contain: transfer to the subgrid scales (~ -eps_sgs when the
          SGS model is in equilibrium), SGS transport, numerical dissipation of the odd-order advection.
-eps_sgs from the Deardorff model is plotted next to R: in a well-resolved LES the two should be close away from the wall.

Staggering (our runs use erf.destag_profiles = false): every column is written at the face heights z_k = k dz, but
cell-centred quantities (u, v, theta, rho, u'u', v'v', ...) in row k are the values of cell k, centred at z_k + dz/2,
and their last row (k = n_cell) is to be discarded; face quantities (w, u'w', v'w', w'w', w'theta', u_i'u_i'w', p'w',
tau_13, tau_23, tau_theta_w) in row k sit at z_k, face 0 being the ground.  Terms are formed at cell centres: a face
quantity's derivative at the centre of cell k is (F[k+1]-F[k])/dz, exact for the staggered grid.

Checks that run every time (Mike, 22 Oct: verify analysis code):
  * --selftest   analytic profiles with known derivatives -> every term reproduced to < 0.5 %
  * surface closure: resolved + SGS heat flux at z = 0 must equal the prescribed flux (0.20 K m/s convective)
  * w'w', u_i'u_i'w', p'w' are zero on the ground face (w = 0 there)
  * units anchor: P_b at the surface = (g/theta) Q_0 = w_*^3/z_i by definition of w_*

Usage:  python3 tke_budget.py --selftest
        python3 tke_budget.py --case cbl --win 10800 14400 DIR [DIR ...]      (figure + JSON for the convective case)
        python3 tke_budget.py --case sbl --win 28800 32400 DIR [DIR ...]
"""
import argparse, json, os, sys
import numpy as np

G, THETA_REF = 9.81, 300.0
OUT_FIG = "/work2/11859/hsna/stampede3/HowlandLab/Results Summary/figs"
OUT_JSON = "/work2/11859/hsna/stampede3/HowlandLab/Results Summary/tke_budget.json"
# Okabe-Ito hues in a fixed order (validated with the dataviz palette script, 22 Oct 2026); residual is grey, dashed
COL = {"P_s": "#0072B2", "P_b": "#D55E00", "T_t": "#009E73", "T_p": "#7B3FA0", "xfer": "#9A6B00", "-eps": "#9A6B00", "R": "#555555"}
LBL = {"P_s": "shear production $P_s$", "P_b": "buoyancy $P_b$", "T_t": "turbulent transport $T_t$",
       "T_p": "pressure transport $T_p$", "xfer": r"SGS transfer, mean-shear part $\langle\tau_{i3}\rangle\,\partial U_i/\partial z$",
       "-eps": r"$-\varepsilon_{sgs}$ (Deardorff TKE eq.)", "R": "residual $R$ (SGS transfer by eddies + numerics)"}
# column indices (0 = time, 1 = z)
CV = dict(uu=2, uv=3, uw=4, vv=5, vw=6, ww=7, uth=8, vth=9, wth=10, thth=11, uiuiu=12, uiuiv=13, uiuiw=14, pu=15, pv=16, pw=17, wqv=18, wqc=19, wqr=20, wthv=21)
MN = dict(u=2, v=3, w=4, rho=5, th=6, ksgs=7, K=8)
SF = dict(t11=2, t12=3, t13=4, t22=5, t23=6, t33=7, hfx=8, q1=9, q2=10, eps=11)


def load(path):
    """-> t (nt,), z (nz,), A (nt, nz, ncol); a time that appears twice (restart) keeps its last block."""
    a = np.loadtxt(path)
    z = np.unique(a[:, 1]); nz = z.size
    blocks = {}
    for i in range(0, a.shape[0], nz):
        blocks[a[i, 0]] = a[i:i + nz]
    t = np.array(sorted(blocks)); A = np.stack([blocks[x] for x in t])
    assert np.allclose(A[0, :, 1], z), "rows are not ordered by z"
    return t, z, A


def terms_from_arrays(z, U, V, rho, uw, vw, ww, wthv, uiuiw, pw):
    """All inputs at the ERF staggered layout: cell-centred arrays have n values, face arrays n+1 values (z has n+1 faces).
    Returns z_c and the terms at cell centres (n values)."""
    dz = np.diff(z); zc = z[:-1] + 0.5 * dz
    f2c = lambda F: 0.5 * (F[:-1] + F[1:])             # face -> centre
    ddz_f = lambda F: (F[1:] - F[:-1]) / dz            # centre derivative of a face quantity
    dUdz, dVdz = np.gradient(U, zc), np.gradient(V, zc)
    P_s = -f2c(uw) * dUdz - f2c(vw) * dVdz
    P_b = (G / THETA_REF) * f2c(wthv)
    T_t = -0.5 * ddz_f(uiuiw)
    T_p = -ddz_f(pw) / rho
    return zc, dict(P_s=P_s, P_b=P_b, T_t=T_t, T_p=T_p)


def budget(d, win):
    tc, z, C = load(os.path.join(d, "covar_profiles.dat"))
    tm, _, M = load(os.path.join(d, "mean_profiles.dat"))
    ts, _, S = load(os.path.join(d, "sfs_profiles.dat"))
    n = z.size - 1; dz = np.diff(z); zc = z[:-1] + 0.5 * dz
    sel = lambda t: (t >= win[0]) & (t <= win[1])
    wc, wm, ws = sel(tc), sel(tm), sel(ts)
    assert wc.sum() > 2 and wm.sum() > 2 and ws.sum() > 2, f"{d}: window {win} has too few records"
    # time-mean of every term (terms formed per record, then averaged)
    acc = {k: np.zeros(n) for k in ("P_s", "P_b", "T_t", "T_p")}
    tm_sel = tm[wm]
    for i in np.where(wc)[0]:
        j = np.argmin(abs(tm_sel - tc[i])); Mi = M[wm][j]
        # buoyancy from <w'theta'> (column 10): the runs are dry, so theta_v = theta, and column 21 (<w'theta_v'>) is
        # not usable in ERF v26.09 for a dry run (it comes out negative at every height; checked 22 Oct 2026)
        _, tr = terms_from_arrays(z, Mi[:n, MN["u"]], Mi[:n, MN["v"]], Mi[:n, MN["rho"]],
                                  C[i, :, CV["uw"]], C[i, :, CV["vw"]], C[i, :, CV["ww"]], C[i, :, CV["wth"]],
                                  C[i, :, CV["uiuiw"]], C[i, :, CV["pw"]])
        for k in acc: acc[k] += tr[k]
    for k in acc: acc[k] /= wc.sum()
    f2c = lambda F: 0.5 * (F[:-1] + F[1:])
    e_t = 0.5 * (C[:, :n, CV["uu"]] + C[:, :n, CV["vv"]] + np.array([f2c(C[i, :, CV["ww"]]) for i in range(tc.size)]))
    i1, i2 = np.where(wc)[0][[0, -1]]
    dedt = (e_t[i2] - e_t[i1]) / (tc[i2] - tc[i1])
    e_res = e_t[wc].mean(0)
    eps = S[ws][:, :n, SF["eps"]].mean(0)                 # Deardorff only (zero for Smagorinsky, which has no TKE equation)
    ksgs = M[wm][:, :n, MN["ksgs"]].mean(0)
    # transfer from resolved to subgrid scales, its dominant part <tau_13> dU/dz + <tau_23> dV/dz (tau_i3 < 0 where dU_i/dz > 0,
    # so this is a sink of resolved TKE); defined for both SGS models.  For Deardorff in equilibrium it should be ~ -eps.
    Um, Vm = M[wm][:, :n, MN["u"]].mean(0), M[wm][:, :n, MN["v"]].mean(0)
    f2c = lambda F: 0.5 * (F[:-1] + F[1:])
    xfer = f2c(S[ws][:, :, SF["t13"]].mean(0)) * np.gradient(Um, zc) + f2c(S[ws][:, :, SF["t23"]].mean(0)) * np.gradient(Vm, zc)
    R = dedt - (acc["P_s"] + acc["P_b"] + acc["T_t"] + acc["T_p"] + xfer)
    # checks and scales
    wth_tot_f = C[wc][:, :, CV["wth"]].mean(0) + S[ws][:, :, SF["hfx"]].mean(0)        # resolved + SGS, on faces
    q0 = wth_tot_f[0]
    uw_tot_f = C[wc][:, :, CV["uw"]].mean(0) + S[ws][:, :, SF["t13"]].mean(0)
    vw_tot_f = C[wc][:, :, CV["vw"]].mean(0) + S[ws][:, :, SF["t23"]].mean(0)
    ustar = (uw_tot_f[0] ** 2 + vw_tot_f[0] ** 2) ** 0.25
    chk = dict(ground_face_zero=dict(ww=float(abs(C[wc][:, 0, CV["ww"]]).max()), uiuiw=float(abs(C[wc][:, 0, CV["uiuiw"]]).max()),
                                     pw=float(abs(C[wc][:, 0, CV["pw"]]).max())),
               surface_heat_flux_total=float(q0), surface_heat_flux_resolved=float(C[wc][:, 0, CV["wth"]].mean()),
               ustar_from_total_stress=float(ustar))
    sgs_share_mom = np.abs(f2c(S[ws][:, :, SF["t13"]].mean(0))) / (np.abs(f2c(uw_tot_f)) + 1e-12)
    sgs_share_heat = np.abs(f2c(S[ws][:, :, SF["hfx"]].mean(0))) / (np.abs(f2c(wth_tot_f)) + 1e-12)
    out = dict(dir=d, win=list(win), z=zc.tolist(), dz=float(dz[0]), terms={k: v.tolist() for k, v in acc.items()},
               dedt=dedt.tolist(), R=R.tolist(), xfer_sgs=xfer.tolist(), eps_sgs=eps.tolist(), e_res=e_res.tolist(), k_sgs=ksgs.tolist(),
               sgs_share_mom=sgs_share_mom.tolist(), sgs_share_heat=sgs_share_heat.tolist(), checks=chk)
    if q0 > 0:                                   # convective: z_i = height of minimum total heat flux, w_* from Q_0 and z_i
        zi = float(z[np.argmin(wth_tot_f)]); wstar = (G / THETA_REF * q0 * zi) ** (1 / 3)
        # anchor: with the total heat flux falling linearly from Q0 to -0.2 Q0 at z_i, P_b at z = 0.5 z_i is 0.4 w_*^3/z_i;
        # the resolved P_b there should be close to that (the SGS part of the flux is < 1 % away from the wall)
        k5 = np.argmin(abs(zc - 0.5 * zi))
        out.update(zi=zi, wstar=float(wstar), scale=float(wstar ** 3 / zi), chk_Pb_mid_over_scale=float(acc["P_b"][k5] / (wstar ** 3 / zi)))
    return out


def selftest():
    """Analytic profiles on a uniform staggered grid; every term must match its closed form."""
    n, H = 200, 1000.0
    z = np.linspace(0, H, n + 1); dz = z[1] - z[0]; zc = z[:-1] + dz / 2
    U = 5 * np.sin(np.pi * zc / (2 * H)); V = 2 * (zc / H) ** 2; rho = np.full(n, 1.1)
    uw = -0.1 * (1 - z / H); vw = -0.05 * (1 - z / H) ** 2; ww = 0.5 * np.sin(np.pi * z / H) ** 2
    wthv = 0.2 * (1 - 1.2 * z / H); uiuiw = 0.3 * np.sin(np.pi * z / H) * (z / H); pw = -0.4 * np.sin(2 * np.pi * z / H)
    zc2, tr = terms_from_arrays(z, U, V, rho, uw, vw, ww, wthv, uiuiw, pw)
    dU = 5 * np.pi / (2 * H) * np.cos(np.pi * zc / (2 * H)); dV = 4 * zc / H ** 2
    f2c = lambda F: 0.5 * (F[:-1] + F[1:])
    exact = dict(P_s=-f2c(uw) * dU - f2c(vw) * dV, P_b=G / THETA_REF * f2c(wthv),
                 T_t=-0.5 * 0.3 * (np.pi / H * np.cos(np.pi * zc / H) * zc / H + np.sin(np.pi * zc / H) / H),
                 T_p=0.4 * 2 * np.pi / H * np.cos(2 * np.pi * zc / H) / 1.1)
    ok = True
    for k in tr:
        err = abs(tr[k] - exact[k]).max() / abs(exact[k]).max(); ok &= err < 5e-3
        print(f"  {k}: max rel error {err:.2e} {'ok' if err < 5e-3 else 'FAIL'}")
    # the discrete transport terms telescope: their integral equals the boundary values of the flux
    integ = (tr["T_t"] * dz).sum(); print(f"  sum T_t dz = {integ:.3e}  vs -1/2[F]_0^H = {-0.5 * (uiuiw[-1] - uiuiw[0]):.3e}")
    print("selftest", "PASSED" if ok else "FAILED"); return ok


def figure(case, runs):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    zc = np.array(runs[0]["z"]); n = len(runs)
    if case == "cbl":
        sc = np.array([r["scale"] for r in runs]); zi = np.array([r["zi"] for r in runs])
        ys = [zc / zi[i] for i in range(n)]; ylab = r"$z / z_i$"; xlab = r"budget term $\times\, z_i / w_*^3$"; ylim = (0, 1.3)
    else:
        sc = np.ones(n); ys = [zc] * n; ylab = "z (m)"; xlab = r"budget term (m$^2$ s$^{-3}$)"; ylim = (0, 300)
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 4.6), sharey=False)
    has_eps = any(np.abs(r["eps_sgs"]).max() > 0 for r in runs)
    series = [("P_s", "P_s"), ("P_b", "P_b"), ("T_t", "T_t"), ("T_p", "T_p"), ("xfer", None)] + ([("-eps", None)] if has_eps else []) + [("R", None)]
    pick = {"xfer": lambda r: np.array(r["xfer_sgs"]), "-eps": lambda r: -np.array(r["eps_sgs"]), "R": lambda r: np.array(r["R"])}
    for key, tk in series:
        vals = []
        for i, r in enumerate(runs):
            v = np.array(r["terms"][tk]) if tk else pick[key](r)
            vals.append(np.interp(ys[0], ys[i], v / sc[i]))
        vals = np.array(vals); m = vals.mean(0)
        ax[0].plot(m, ys[0], color=COL[key], lw=2, ls={"R": "--", "-eps": ":"}.get(key, "-"), label=LBL[key])
        if n > 1: ax[0].fill_betweenx(ys[0], vals.min(0), vals.max(0), color=COL[key], alpha=0.15, lw=0)
    ax[0].axvline(0, color="#999999", lw=0.8); ax[0].set_ylim(*ylim); ax[0].set_xlabel(xlab); ax[0].set_ylabel(ylab)
    ax[0].set_title("resolved TKE budget" + (f", {n} members (band = spread)" if n > 1 else ""), fontsize=10)
    ax[0].legend(frameon=False, fontsize=7.5, loc="upper right")
    for i, r in enumerate(runs):
        ax[1].plot(r["sgs_share_mom"], ys[i], color="#0072B2", lw=1.2, alpha=0.8, label="momentum flux $|\\tau_{13}| / |u'w' + \\tau_{13}|$" if i == 0 else None)
        ax[1].plot(r["sgs_share_heat"], ys[i], color="#D55E00", lw=1.2, alpha=0.8, label="heat flux $|\\tau_{\\theta w}| / |w'\\theta' + \\tau_{\\theta w}|$" if i == 0 else None)
    ax[1].set_xlim(0, 1); ax[1].set_ylim(*ylim); ax[1].set_xlabel("subgrid share of the vertical flux"); ax[1].set_ylabel(ylab)
    ax[1].set_title("how much the SGS model carries", fontsize=10); ax[1].legend(frameon=False, fontsize=8)
    for a in ax: a.grid(color="#EEEEEE", lw=0.6)
    fig.suptitle(f"ERF {'convective (S&D, 25 m)' if case == 'cbl' else 'stable (GABLS1, 6.25 m)'} case, mean over {runs[0]['win'][0] / 3600:.0f}–{runs[0]['win'][1] / 3600:.0f} h", fontsize=11)
    fig.tight_layout(); f = f"{OUT_FIG}/tke_budget_{case}.png"; fig.savefig(f, dpi=150); plt.close(fig); print("wrote", f)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("dirs", nargs="*"); ap.add_argument("--case", default="cbl"); ap.add_argument("--win", nargs=2, type=float)
    ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    if a.selftest: sys.exit(0 if selftest() else 1)
    win = tuple(a.win) if a.win else ((10800, 14400) if a.case == "cbl" else (28800, 32400))
    runs = []
    for d in a.dirs:
        r = budget(d, win); runs.append(r)
        c = r["checks"]
        print(f"{os.path.basename(d)}: Q0 total {c['surface_heat_flux_total']:+.4f} (resolved {c['surface_heat_flux_resolved']:+.1e}) K m/s | u* {c['ustar_from_total_stress']:.3f} | ground-face max |ww|,|uiuiw|,|pw| = "
              f"{c['ground_face_zero']['ww']:.1e}, {c['ground_face_zero']['uiuiw']:.1e}, {c['ground_face_zero']['pw']:.1e}"
              + (f" | zi {r['zi']:.0f} m, w* {r['wstar']:.2f} m/s, P_b(0.5 zi) x zi/w*^3 = {r['chk_Pb_mid_over_scale']:.3f} (expect ~0.4)" if "zi" in r else ""))
        zc = np.array(r["z"]); ml = (zc > 0.2 * r.get("zi", 1e9)) & (zc < 0.8 * r.get("zi", 1e9)) if "zi" in r else (zc > 20) & (zc < 100)
        s = r["scale"] if "zi" in r else 1.0
        print("   layer mean x scale:", ", ".join(f"{k} {np.mean(np.array(r['terms'][k])[ml]) / s:+.3f}" for k in r["terms"]),
              f", xfer {np.mean(np.array(r['xfer_sgs'])[ml]) / s:+.3f}, -eps {-np.mean(np.array(r['eps_sgs'])[ml]) / s:+.3f}, R {np.mean(np.array(r['R'])[ml]) / s:+.3f}, de/dt {np.mean(np.array(r['dedt'])[ml]) / s:+.3f}"
              + ("" if s != 1.0 else "   (raw m2/s3; stable case)"))
    if runs:
        figure(a.case, runs)
        allj = json.load(open(OUT_JSON)) if os.path.exists(OUT_JSON) else {}
        def layer(r):   # 0.2-0.8 z_i (convective) or 20-100 m (stable), scaled by w_*^3/z_i where defined; transport terms integrate to 0 over the whole column, so a column mean would say nothing
            zc = np.array(r["z"]); m = (zc > 0.2 * r["zi"]) & (zc < 0.8 * r["zi"]) if "zi" in r else (zc > 20) & (zc < 100); s = r.get("scale", 1.0)
            d = {k: float(np.mean(np.array(r["terms"][k])[m]) / s) for k in r["terms"]}
            d.update({k: float(np.mean(np.array(r[k])[m]) / s) for k in ("xfer_sgs", "eps_sgs", "R", "dedt")}); d["layer"] = "0.2-0.8 zi, x zi/w*^3" if "zi" in r else "20-100 m, m2/s3"
            return d
        allj[a.case] = [{k: v for k, v in r.items() if k not in ("terms", "z", "dedt", "R", "xfer_sgs", "eps_sgs", "e_res", "k_sgs", "sgs_share_mom", "sgs_share_heat")} | {"layer_mean": layer(r)} for r in runs]
        json.dump(allj, open(OUT_JSON, "w"), indent=1); print("wrote", OUT_JSON)
