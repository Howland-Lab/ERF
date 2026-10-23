#!/usr/bin/env python3
"""ens_wrf_vs_erf.py -- WRF ensembles (WRF_SEED patch, 22-23 Oct 2026) against the ERF ensembles: member spread of u*, the
boundary-layer depth, and the resolved-w statistics at 0.5 z_i (convective) / 0.5 h (stable), last hour.

WRF members: postproc/data_hr/<name>.npz written by postproc_highres.run_one (planes stats at the ZREL levels, same bins and
window as the ERF side).  ERF members: analysis/week2_summary.json (ens_cbl: 4 seeds + 2 base runs; ens_sbl: 5 seeds).
Members that are not processed yet are skipped, so the script can be re-run as they finish.
Writes Results Summary/figs/ens_wrf_vs_erf.png and Results Summary/ens_wrf_vs_erf.json.
"""
import json, os, glob
import numpy as np

H = "/work2/11859/hsna/stampede3/HowlandLab"
DATA = f"{H}/WRF-4.6.1/test/abl_cases/postproc/data_hr"
WRF = {"cbl": ["cbl_sd25_fix", "ens_cbl25_s1", "ens_cbl25_s2", "ens_cbl25_s3", "ens_cbl25_s4"],
       "sbl": ["sbl_dx6p25", "ens_sbl6_s1", "ens_sbl6_s2", "ens_sbl6_s3", "ens_sbl6_s4", "ens_sbl6_s5"]}
HKEY = {"cbl": "h_minflux", "sbl": "h_stress"}
STATS = ["ust", "h", "w_std", "w_skew", "w_kurt", "w_p3"]
LBL = {"ust": "$u_*$ (m/s)", "h": "depth (m)", "w_std": "$\\sigma_w$ at 0.5 $z_i$ (m/s)", "w_skew": "skewness", "w_kurt": "kurtosis", "w_p3": "$P(|w'|>3\\sigma)$, window σ (both)"}


def wrf_member(name, case):
    f = f"{DATA}/{name}.npz"
    if not os.path.exists(f): return None
    d = np.load(f, allow_pickle=True); M = d["mean"].item(); P = d["planes"].item(); w = P["0.50_w"]
    # P(|w'| > 3σ) with σ = the window standard deviation, from the stored histogram (bins in m/s).  The stored "tail" counts
    # use the FIRST frame's σ as the threshold (postproc_highres and erf_planes alike); in the stable case that σ differs
    # by ±10 % between frames and the 3σ count then scatters 4x between members (0.004–0.022) although the window-σ
    # tails agree (0.005–0.007).  Found 23 Oct 2026.  The ERF ensemble values (week2_analysis.w_stats) already use the
    # window σ (p3 from allw.std()), so with this recomputation the two sides are like-for-like.
    b = np.linspace(-8, 8, 321); c = 0.5 * (b[:-1] + b[1:]); h = np.asarray(w["hist"], float)
    p3 = float(h[np.abs(c) > 3 * w["std"]].sum() / h.sum())
    return dict(name=name, ust=float(M["ust"]), h=float(M[HKEY[case]]), w_std=float(w["std"]), w_skew=float(w["skew"]), w_kurt=float(w["kurt"]), w_p3=float(p3))


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    E = json.load(open(f"{H}/ERF_runs/analysis/week2_summary.json")); out = {}
    fig, ax = plt.subplots(2, len(STATS), figsize=(14, 5.6))
    for r, case in enumerate(("cbl", "sbl")):
        W = [m for m in (wrf_member(n, case) for n in WRF[case]) if m]; e = E["ens_" + case]
        out[case] = dict(wrf_members=W, erf=e)
        print(f"== {case}: WRF {len(W)} members processed ({', '.join(m['name'] for m in W)}), ERF {len(e['members'])} members")
        for j, s in enumerate(STATS):
            wv = np.array([m[s] for m in W]) if W else np.array([])
            ev = np.array(e[s]["values"]) if isinstance(e[s], dict) and "values" in e[s] else None
            em, es = (e[s]["mean"], e[s]["std"]) if isinstance(e[s], dict) else (np.nan, np.nan)
            a = ax[r, j]
            if ev is not None: a.plot(np.zeros(len(ev)), ev, "o", color="#0072B2", ms=6, alpha=0.8, label="ERF members")
            else: a.errorbar([0], [em], yerr=[es], fmt="o", color="#0072B2", ms=6, capsize=4, label="ERF mean ± std")
            if len(wv): a.plot(np.ones(len(wv)), wv, "s", color="#222222", ms=6, alpha=0.8, label="WRF members")
            a.set_xticks([0, 1]); a.set_xticklabels(["ERF", "WRF"]); a.set_xlim(-0.6, 1.6); a.set_title(LBL[s], fontsize=9); a.grid(color="#EEEEEE", lw=0.6)
            if len(wv): print(f"   {s:7s}: ERF {em:.4g} ± {es:.2g} (cv {100 * es / abs(em) if em else 0:.1f} %) | WRF {wv.mean():.4g} ± {wv.std():.2g} (cv {100 * wv.std() / abs(wv.mean()):.1f} %, n={len(wv)})")
            if r == 0 and j == 0: a.legend(frameon=False, fontsize=7)
        ax[r, 0].set_ylabel("convective 25 m, 3–4 h" if case == "cbl" else "stable 6.25 m, 8–9 h", fontsize=9)
    fig.suptitle("Ensemble spread, member by member: ERF (random_seed) vs WRF (WRF_SEED patch, 22 Oct)", fontsize=11); fig.tight_layout()
    fig.savefig(f"{H}/Results Summary/figs/ens_wrf_vs_erf.png", dpi=150)
    json.dump(out, open(f"{H}/Results Summary/ens_wrf_vs_erf.json", "w"), indent=1, default=float); print("wrote figure + json")
