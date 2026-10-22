#!/usr/bin/env python3
"""make_ar_runs.py -- grid aspect-ratio test for the VLES scan (22 Oct 2026, after Mike's remark that the LES SGS model is
isotropic, so dx/dz matters).  cases/gz_cbl_dx100 has dx = 100 m, dz = 20 m (aspect ratio 5); here the same case is run
with dz = 50 m (ratio 2) and dz = 100 m (ratio 1, cubic cells).  If sigma_w at 100 m moves with dz, part of the
"+5/+11/+18 % sigma_w" of the scan is an aspect-ratio effect, not a Delta/zi effect.  Everything else is unchanged
(16 km domain, 2 km top, Deardorff, dt = 0.25 s, 4 h).  Usage: python3 make_ar_runs.py [--submit]
"""
import os, re, shutil, subprocess, sys
R = os.path.dirname(os.path.abspath(__file__)); C = os.path.join(R, "cases")
NZ = {50: 40, 100: 20}


def make(dz):
    name = f"gz_cbl_dx100_dz{dz:03d}"; d = os.path.join(C, name); os.makedirs(d, exist_ok=True)
    shutil.copy(os.path.join(C, "gz_cbl_dx100", "input_sounding"), d)
    t = open(os.path.join(C, "gz_cbl_dx100", "inputs")).read()
    t = re.sub(r"^amr\.n_cell\s*=.*$", f"amr.n_cell = 160 160 {NZ[dz]}   # dx = 100 m, dz = {dz} m (aspect-ratio test)", t, count=1, flags=re.M)
    t = f"# {name}: aspect-ratio test, dx = 100 m, dz = {dz} m, dx/dz = {100 // dz} (make_ar_runs.py, 22 Oct 2026)\n" + t
    open(os.path.join(d, "inputs"), "w").write(t)
    return d


if __name__ == "__main__":
    dirs = [make(dz) for dz in NZ]
    print([os.path.basename(x) for x in dirs])
    if "--submit" in sys.argv:
        r = subprocess.run(["sbatch", "-J", "gz_ar", "-p", "spr", "-N", "1", "-n", "64", "-t", "02:30:00",
                            os.path.join(R, "submit_erf_chain.slurm")] + dirs, capture_output=True, text=True)
        print(r.stdout.strip() or r.stderr[-300:])
