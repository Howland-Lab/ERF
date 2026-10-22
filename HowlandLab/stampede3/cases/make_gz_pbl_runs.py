#!/usr/bin/env python3
"""make_gz_pbl_runs.py -- grey-zone closure scan in ERF (21 Oct 2026), the ERF counterpart of
WRF-4.6.1/test/abl_cases/make_gz_scheme_runs.py.  Convective case, 16 km x 16 km x 2 km, dz = 20 m, 4 h,
dx = 100, 200, 400 m, with the PBL schemes ERF offers (vertical) + 2-D Smagorinsky (horizontal):
  mynn25    : erf.pbl_type = MYNN25
  mynnedmf  : erf.pbl_type = MYNNEDMF            (mass flux, like WRF's default bl_mynn_edmf=1)
  mrf       : erf.pbl_type = MRF                 (Hong & Pan 1996 K-profile)
  mrfblend  : MRF + erf.pbl_blend_length = 750   (Boutle et al. 2014 scale-aware blending of K_h; ERF applies it to
                                                  theta and q_v only, momentum K is not blended)
Also cases/rans_cbl_1km_edmf (1 km, MYNNEDMF).  Base decks: cases/gz_cbl_dx{100,200,400} (VLES scan) and rans_cbl_1km.
MYNN trips amrex.fpe_trap_invalid at step 1, so the trap is off (as in ERF's own inputs_GABLS1_mynn25).
Usage: python3 make_gz_pbl_runs.py [--submit]
"""
import os, re, shutil, subprocess, sys

R = os.path.dirname(os.path.abspath(__file__)); C = os.path.join(R, "cases")
SCHEMES = {"mynn25": ("MYNN25", 0), "mynnedmf": ("MYNNEDMF", 0), "mrf": ("MRF", 0), "mrfblend": ("MRF", 750)}
MGS = {100: 20, 200: 10, 400: 10}            # boxes: 64, 64, 16
RANKS = {100: 64, 200: 64, 400: 16}; TLIM = {100: "03:30:00", 200: "01:30:00", 400: "00:45:00"}


def setkey(txt, key, val, comment=""):
    pat = re.compile(r"^(%s\s*=\s*)[^\n#]*(#.*)?$" % re.escape(key), re.M)
    if pat.search(txt):
        return pat.sub(lambda m: f"{m.group(1)}{val}" + (f"   # {comment}" if comment else ""), txt, count=1)
    return txt + f"\n{key} = {val}" + (f"   # {comment}" if comment else "") + "\n"


def make(dx, scheme):
    pbl, Lb = SCHEMES[scheme]; name = f"gzpbl_{scheme}_dx{dx}"; d = os.path.join(C, name); os.makedirs(d, exist_ok=True)
    src = os.path.join(C, f"gz_cbl_dx{dx:03d}")
    shutil.copy(os.path.join(src, "input_sounding"), d)
    t = open(os.path.join(src, "inputs")).read()
    t = t.replace("# ERF convective ABL (Shin & Dudhia 2016 set-up)", f"# {name}: grey-zone closure scan, {pbl}{' + Boutle blending L=750 m' if Lb else ''}, dx = {dx} m (make_gz_pbl_runs.py, 21 Oct 2026)\n# ERF convective ABL (Shin & Dudhia 2016 set-up)", 1)
    t = setkey(t, "amrex.fpe_trap_invalid", "0", "MYNN trips the trap at step 1; profiles checked for NaN")
    t = setkey(t, "amr.max_grid_size_x", str(MGS[dx])); t = setkey(t, "amr.max_grid_size_y", str(MGS[dx]))
    t = setkey(t, "erf.fixed_dt", f"{dx / 100:.2f}", "dt = dx/100 s, acoustic Courant ~0.6 with ratio 6 (same as WRF scan)")
    t = setkey(t, "erf.profile_int", str(int(60 / (dx / 100))), "every 60 s")
    t = setkey(t, "erf.les_type", '"Smagorinsky2D"', "horizontal 2-D Smagorinsky, like WRF km_opt=4")
    t = setkey(t, "erf.Cs", "0.25")
    t = setkey(t, "erf.pbl_type", pbl, "vertical transport by the PBL scheme")
    t = setkey(t, "prob.KE_0", "0.4"); t = setkey(t, "prob.KE_decay_height", "1000."); t = setkey(t, "prob.KE_decay_order", "3")
    if Lb:
        t = setkey(t, "erf.pbl_blend_length", str(Lb), "Boutle et al. 2014: f = dx^2/(L^2+dx^2) multiplies K_h (theta, qv)")
    t = setkey(t, "erf.plot_vars_1", "x_velocity y_velocity z_velocity theta rhoKE Kmv Khv")
    open(os.path.join(d, "inputs"), "w").write(t)
    return d


def make_1km_edmf():
    d = os.path.join(C, "rans_cbl_1km_edmf"); os.makedirs(d, exist_ok=True)
    shutil.copy(os.path.join(C, "rans_cbl_1km", "input_sounding"), d)
    t = open(os.path.join(C, "rans_cbl_1km", "inputs")).read()
    t = setkey(t, "erf.pbl_type", "MYNNEDMF", "mass-flux version, like WRF's default bl_mynn_edmf=1")
    open(os.path.join(d, "inputs"), "w").write("# rans_cbl_1km_edmf: as rans_cbl_1km but MYNN-EDMF (21 Oct 2026)\n" + t)
    return d


if __name__ == "__main__":
    submit = "--submit" in sys.argv
    chain = os.path.join(R, "submit_erf_chain.slurm")
    jobs = {dx: [make(dx, s) for s in SCHEMES] for dx in MGS}
    jobs["1km"] = [make_1km_edmf()]
    for key, dirs in jobs.items():
        print(key, [os.path.basename(x) for x in dirs])
        if submit:
            n, tl = (RANKS[key], TLIM[key]) if key in RANKS else (25, "00:40:00")
            r = subprocess.run(["sbatch", "-J", f"gzpbl_{key}", "-p", "spr", "-N", "1", "-n", str(n), "-t", tl, chain] + dirs, capture_output=True, text=True)
            print("  ", [l for l in r.stdout.splitlines() if "Submitted" in l] or r.stderr[-300:])
