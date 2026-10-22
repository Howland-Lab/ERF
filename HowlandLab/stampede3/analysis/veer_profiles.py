#!/usr/bin/env python3
"""veer_profiles.py -- wind speed and direction (veer) profiles, WRF vs ERF, compared as ensemble means (Mike, 22 Oct 2026:
"plot wind direction / veer; compare ensemble averages alongside the PDFs").

ERF: mean_profiles.dat (<u>, <v> at cell centres, erf.destag_profiles = false -> row k is cell k at z_k + dz/2), last hour.
WRF: planes_d01_* (31 snapshots of the last hour, 3-D U, V on the C-grid) -> destaggered, plane-averaged, time-averaged;
     heights from ZCOL in colout_d01_* (flat ideal case, so one column's heights are everyone's).
Direction is the meteorological one, deg the wind blows FROM (u = +10, v = 0 -> 270 deg).  The geostrophic wind of both
cases is along +x, so the direction must tend to 270 deg at the top: printed as a check.  Veer = direction minus the
direction at the lowest level (positive = clockwise with height, the Northern-Hemisphere Ekman sense).

Usage: python3 veer_profiles.py      (writes Results Summary/figs/veer_profiles.png and Results Summary/veer_profiles.json)
"""
import json, os, glob
import numpy as np

H = "/work2/11859/hsna/stampede3/HowlandLab"
OUT_FIG = f"{H}/Results Summary/figs/veer_profiles.png"; OUT_JSON = f"{H}/Results Summary/veer_profiles.json"
CASES = {
    "cbl": dict(label="convective (S&D, 25 m), 3–4 h", win=(10800, 14400), zmax=1300, zi=1040.0,
                erf=sorted(glob.glob(f"{H}/ERF_runs/cases_pvc/ens_cbl25_s*")), wrf=f"{H}/WRF-4.6.1/test/abl_cases/cbl_sd25_fix",
                wrf_planes="planes_d01_0001-01-01_03_00_00"),
    "sbl": dict(label="stable (GABLS1, 6.25 m), 8–9 h", win=(28800, 32400), zmax=300, zi=None,
                erf=sorted(glob.glob(f"{H}/ERF_runs/cases/ens_sbl6_s*")), wrf=f"{H}/WRF-4.6.1/test/abl_cases/sbl_dx6p25",
                wrf_planes="planes_d01_0001-01-01_08_00_00"),
}
COL = {"ERF": "#0072B2", "WRF": "#222222"}


def met_dir(u, v):
    return (270.0 - np.degrees(np.arctan2(v, u))) % 360.0


def erf_profile(d, win):
    a = np.loadtxt(os.path.join(d, "mean_profiles.dat")); z = np.unique(a[:, 1]); nz = z.size
    t = a[::nz, 0]; A = a.reshape(t.size, nz, -1); w = (t >= win[0]) & (t <= win[1])
    u, v = A[w, :, 2].mean(0), A[w, :, 3].mean(0); n = nz - 1          # last row of a cell-centred column is unused
    zc = z[:-1] + 0.5 * np.diff(z)
    return zc, u[:n], v[:n]


def wrf_profile(d, planes):
    import netCDF4 as nc
    g = nc.Dataset(os.path.join(d, planes)); nt = len(g.dimensions["Time"]); nz = len(g.dimensions["bottom_top"])
    u = np.zeros(nz); v = np.zeros(nz)
    for it in range(nt):                                            # one snapshot at a time (login-node memory)
        U = g["U"][it]; V = g["V"][it]                                # (nz, ny, nx+1), (nz, ny+1, nx)
        u += 0.5 * (U[:, :, :-1] + U[:, :, 1:]).mean(axis=(1, 2)); v += 0.5 * (V[:, :-1, :] + V[:, 1:, :]).mean(axis=(1, 2))
    u /= nt; v /= nt
    c = nc.Dataset(os.path.join(d, "colout_d01_0001-01-01_00_00_00")); zc = np.asarray(c["ZCOL"][:]).reshape(-1, nz).mean(0)
    return zc, u, v


def summarize(z, u, v, zi):
    sp = np.hypot(u, v); dr = met_dir(u, v); veer = (dr - dr[0] + 180) % 360 - 180
    k = np.argmin(abs(z - (0.5 * zi if zi else 100.0)))
    return dict(dir_lowest=float(dr[0]), dir_top=float(dr[-1]), speed_max=float(sp.max()), z_speed_max=float(z[np.argmax(sp)]),
                veer_mid=float(veer[k]), z_mid=float(z[k]), veer_max=float(veer[np.argmax(abs(veer))]))


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(9.2, 7.6)); out = {}
    for row, (case, C) in enumerate(CASES.items()):
        prof = [erf_profile(d, C["win"]) for d in C["erf"]]; z = prof[0][0]
        U = np.array([np.interp(z, p[0], p[1]) for p in prof]); V = np.array([np.interp(z, p[0], p[2]) for p in prof])
        sp = np.hypot(U, V); dr = met_dir(U, V)
        zw, uw, vw = wrf_profile(C["wrf"], C["wrf_planes"]); spw, drw = np.hypot(uw, vw), met_dir(uw, vw)
        y = (lambda h: h / C["zi"]) if C["zi"] else (lambda h: h); ylab = "$z/z_i$" if C["zi"] else "z (m)"; ylim = (0, C["zmax"] / C["zi"] if C["zi"] else C["zmax"])
        a = ax[row, 0]; a.plot(sp.mean(0), y(z), color=COL["ERF"], lw=2, label=f"ERF, {len(prof)}-member mean (band = spread)")
        a.fill_betweenx(y(z), sp.min(0), sp.max(0), color=COL["ERF"], alpha=0.15, lw=0)
        a.plot(spw, y(zw), color=COL["WRF"], lw=2, ls="--", label="WRF (one run)"); a.set_xlabel("wind speed (m/s)"); a.set_ylabel(ylab); a.set_ylim(*ylim)
        a.set_title(C["label"], fontsize=10); a.legend(frameon=False, fontsize=8, loc="lower right")
        a = ax[row, 1]; a.plot(dr.mean(0), y(z), color=COL["ERF"], lw=2); a.fill_betweenx(y(z), dr.min(0), dr.max(0), color=COL["ERF"], alpha=0.15, lw=0)
        a.plot(drw, y(zw), color=COL["WRF"], lw=2, ls="--"); a.axvline(270, color="#999999", lw=0.8); a.set_xlabel("wind direction (deg, from); geostrophic = 270"); a.set_ylabel(ylab); a.set_ylim(*ylim)
        a.set_title("veer: direction turns with height", fontsize=10)
        for aa in ax[row]: aa.grid(color="#EEEEEE", lw=0.6)
        out[case] = dict(ERF=[summarize(z, U[i], V[i], C["zi"]) for i in range(len(prof))], WRF=summarize(zw, uw, vw, C["zi"]),
                         ERF_members=[os.path.basename(d) for d in C["erf"]], WRF_run=os.path.basename(C["wrf"]))
        e = out[case]["ERF"]; w = out[case]["WRF"]
        print(f"{case}: top direction ERF {np.mean([x['dir_top'] for x in e]):.1f}, WRF {w['dir_top']:.1f} (expect ~270) | lowest-level direction ERF "
              f"{np.mean([x['dir_lowest'] for x in e]):.1f} ± {np.std([x['dir_lowest'] for x in e]):.1f}, WRF {w['dir_lowest']:.1f} | veer at {e[0]['z_mid']:.0f} m: ERF "
              f"{np.mean([x['veer_mid'] for x in e]):+.1f} ± {np.std([x['veer_mid'] for x in e]):.1f}, WRF {w['veer_mid']:+.1f} deg | speed max ERF {np.mean([x['speed_max'] for x in e]):.2f} at {np.mean([x['z_speed_max'] for x in e]):.0f} m, WRF {w['speed_max']:.2f} at {w['z_speed_max']:.0f} m")
    fig.suptitle("Wind speed and direction, WRF vs ERF (last hour, plane and time mean)", fontsize=11); fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=150); json.dump(out, open(OUT_JSON, "w"), indent=1); print("wrote", OUT_FIG, OUT_JSON)
