#!/usr/bin/env python3
"""verify_analysis.py -- known-answer tests for the analysis code behind the 22 Oct deck and the results blog
(Mike, 22 Oct 2026: "a result from new code is only as good as the test that code has passed").

    python3 verify_analysis.py            # synthetic inputs with analytic answers; exit 1 if anything fails
    python3 verify_analysis.py --real     # + a few checks on real run output (needs the run directories)

What is tested, and where the tested code lives
  spectra.py              spectrum1d / spectrum2d / line_share / tile_mode_share (the tiling diagnosis)
  erf_planes.Moments      streaming mean/std/skewness/kurtosis/tails/histogram (same class in postproc_highres.py)
  planes_clean_fig        merged_pdf: merged histogram bins with masking of bins that hold < 5 samples
  nwp_vs_les_tails        exceed_from_hist, gauss_exceed
  gz_filter_compare       box_average (fair test: 25 m LES averaged to Delta); same function in postproc_highres
  postproc_highres        autocorr, integral_scale, block_stats, h_minflux, h_stress, h_gradient
  tke_budget              selftest() (analytic budget terms) + a synthetic profile-log case + the wthv column trap

Conventions fixed by these tests (see spectra.py): 1-D spectrum one-sided in k = 2 pi m / L with sum(E) dk = var;
2-D spectrum |FFT2|^2 / N^2 with sum = var and the DC mode zeroed.

Bugs found by this suite are listed in the CHANGELOG (ERF_lab/HowlandLab/stampede3/CHANGELOG.md, 2026-10-22).
"""
import inspect, os, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RS = "/work2/11859/hsna/stampede3/HowlandLab/Results Summary"
WRFP = "/work2/11859/hsna/stampede3/HowlandLab/WRF-4.6.1/test/abl_cases/postproc"
for p in (HERE, RS, WRFP):
    sys.path.insert(0, p)
import matplotlib
matplotlib.use("Agg")
from spectra import spectrum1d, spectrum2d, line_share, tile_mode_share
from erf_planes import Moments, HBINS
import postproc_highres as PH
from planes_clean_fig import merged_pdf
from nwp_vs_les_tails import exceed_from_hist, gauss_exceed
from gz_filter_compare import box_average, stats as filter_stats
import tke_budget as TB

RES = []
rng = np.random.default_rng(20261022)


def check(name, ok, detail=""):
    RES.append((name, bool(ok))); print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   ({detail})" if detail else ""), flush=True)


def section(t):
    print(f"\n== {t}")


# ----------------------------------------------------------------------------------------------- 1. spectra
def test_spectra():
    section("spectra.py")
    n, dx = 4096, 10.0; L = n * dx; x = np.arange(n) * dx
    m0 = 37; A = 1.7
    k, E = spectrum1d(A * np.sin(2 * np.pi * m0 * x / L + 0.3), dx)
    dk = k[1] - k[0]; peak = int(np.argmax(E))
    check("1-D sine: single peak at k0", abs(k[peak] - 2 * np.pi * m0 / L) < 1e-9 * k[peak] and np.delete(E, peak).max() < 1e-6 * E[peak],
          f"m0 {m0} -> bin {peak + 1}, off-peak/peak {np.delete(E, peak).max() / E[peak]:.1e}")
    check("1-D sine: Parseval, sum(E) dk = var = A^2/2", abs(E.sum() * dk - A ** 2 / 2) < 1e-9 * A ** 2)
    y = rng.normal(size=n) * 2.0 + 5.0; k, E = spectrum1d(y, dx)
    check("1-D white noise: Parseval (mean removed)", abs(E.sum() * dk - y.var()) < 1e-9 * y.var())
    lo, hi = E[: E.size // 2].mean(), E[E.size // 2:].mean()
    check("1-D white noise: flat (low-k mean / high-k mean within 10 %)", abs(lo / hi - 1) < 0.10, f"ratio {lo / hi:.3f}")
    # odd length: Parseval must still hold (Nyquist handling)
    k, E = spectrum1d(rng.normal(size=4097), dx); dk = k[1] - k[0]
    check("1-D odd n: Parseval", abs(E.sum() * dk - 1.0) < 0.2, f"sum E dk = {E.sum() * dk:.3f} vs var ~ 1 (sampling)")
    # 2-D: one wavevector, Parseval
    ny, nx = 160, 320; X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    a = 0.8 * np.cos(2 * np.pi * (5 * X / nx + 7 * Y / ny)) + 3.0
    P = spectrum2d(a); kx = (np.fft.fftfreq(nx) * nx).round().astype(int); ky = (np.fft.fftfreq(ny) * ny).round().astype(int)
    big = P > 1e-9 * P.max(); modes = sorted({(int(ky[j]), int(kx[i])) for j, i in zip(*np.where(big))})
    check("2-D cosine: energy only at (+-7, +-5)", modes == [(-7, -5), (7, 5)], f"modes {modes}")
    check("2-D cosine: Parseval, sum P = var = A^2/2", abs(P.sum() - 0.8 ** 2 / 2) < 1e-12)
    g = rng.normal(size=(ny, nx)); P = spectrum2d(g)
    check("2-D noise: Parseval", abs(P.sum() - g.var()) < 1e-9 * g.var())
    # the WRF tiling: 16 identical copies along y of a 20-row random field -> energy only on ky = n * 16
    base = rng.normal(size=(20, nx)); tiled = np.tile(base, (16, 1)); P = spectrum2d(tiled)
    s16 = line_share(P, 16, axis=0); s15 = line_share(P, 15, axis=0); s16z = line_share(P, 16, axis=0, include_zero=True)
    check("16 copies along y: all energy on ky = n*16 lines (incl. n = 0), nothing off them", abs(s16z - 1) < 1e-12 and s15 < 1e-12, f"share incl. 0: {s16z:.12f}, share(15) {s15:.1e}")
    check("16 copies along y: diagnosis share (n != 0) = 15/16", abs(s16 - 15 / 16) < 0.02, f"share(16) {s16:.4f}")
    shifted = np.vstack([np.roll(base, int(rng.integers(0, nx)), axis=1) for _ in range(16)]); P2 = spectrum2d(shifted)
    s16b = line_share(P2, 16, axis=0)
    check("same copies shifted at random in x: no comb (share ~ 1/16)", 0.02 < s16b < 0.15, f"share(16) {s16b:.3f}")
    check("tile_mode_share: pure (px,py) mode -> 1", abs(tile_mode_share(spectrum2d(np.cos(2 * np.pi * (23 * X / nx + 16 * Y / ny))), 23, 16) - 1) < 1e-9)


# ----------------------------------------------------------------------------------------------- 2. moments, PDFs, tails
def test_moments():
    section("erf_planes.Moments / postproc_highres.Moments (same code)")
    N = 1_000_000; g = rng.normal(size=N)
    for cls, name in ((Moments, "erf_planes"), (PH.Moments, "postproc_highres")):
        M = cls(); M.add(g[: N // 3], HBINS["w"], 1.0); M.add(g[N // 3: 2 * N // 3], HBINS["w"], 1.0); M.add(g[2 * N // 3:], HBINS["w"], 1.0)
        r = M.result()
        check(f"{name}: Gaussian mean 0, std 1 (0.5 %)", abs(r["mean"]) < 5e-3 and abs(r["std"] - 1) < 5e-3, f"mean {r['mean']:+.4f} std {r['std']:.4f}")
        check(f"{name}: Gaussian skewness 0 +- 0.01", abs(r["skew"]) < 0.01, f"{r['skew']:+.4f}")
        check(f"{name}: Gaussian kurtosis 3 +- 0.03", abs(r["kurt"] - 3) < 0.03, f"{r['kurt']:.4f}")
        check(f"{name}: P(w > 3 sigma) = 1.35e-3 +- 10 %", abs(r["tail"][2] / 1.3499e-3 - 1) < 0.10, f"{r['tail'][2]:.3e}")
        check(f"{name}: P(w < -3 sigma) = 1.35e-3 +- 10 %", abs(r["tail"][3] / 1.3499e-3 - 1) < 0.10, f"{r['tail'][3]:.3e}")
        check(f"{name}: histogram holds every sample", r["hist"].sum() == N and r["n"] == N)
        one = cls(); one.add(g, HBINS["w"], 1.0); r1 = one.result()
        check(f"{name}: streaming in 3 chunks == one shot", all(abs(r[k] - r1[k]) < 1e-12 for k in ("mean", "std", "skew", "kurt")))
    e = rng.exponential(size=N) - 1.0; M = Moments(); M.add(e, HBINS["w"]); r = M.result()
    check("Moments: exponential-1 -> skewness 2 +- 0.05, kurtosis 9 +- 0.4", abs(r["skew"] - 2) < 0.05 and abs(r["kurt"] - 9) < 0.4, f"skew {r['skew']:.3f} kurt {r['kurt']:.2f}")
    section("planes_clean_fig.merged_pdf")
    s = 1.3; w = rng.normal(size=N) * s; M = Moments(); M.add(w, HBINS["w"]); r = M.result()
    x, pdf = merged_pdf(r["hist"], r["n"], r["std"]); dwr = (HBINS["w"][1] - HBINS["w"][0]) * 4
    full = np.asarray(pdf.filled(0.0)) + np.where(pdf.mask, np.asarray(pdf.data), 0.0)
    check("merged_pdf: total probability = 1 (unmasked)", abs(full.sum() * dwr / r["std"] - 1) < 1e-12)
    check("merged_pdf: masked bins hide < 1e-4 of the probability", 1 - pdf.filled(0.0).sum() * dwr / r["std"] < 1e-4, f"hidden {1 - pdf.filled(0.0).sum() * dwr / r['std']:.1e}")
    mean_h = (x * full).sum() * dwr / r["std"]; var_h = (x ** 2 * full).sum() * dwr / r["std"]
    check("merged_pdf: mean and variance of the merged histogram unbiased (0.01, 1 %)", abs(mean_h) < 0.01 and abs(var_h - 1) < 0.01, f"mean {mean_h:+.4f} var {var_h:.4f}")
    cnt = np.asarray(r["hist"], float).reshape(-1, 4).sum(1)
    check("merged_pdf: mask is exactly 'fewer than 5 samples'", bool(np.all(pdf.mask == (cnt < 5))))
    section("nwp_vs_les_tails.exceed_from_hist / gauss_exceed")
    xe, pe = exceed_from_hist(r["hist"], HBINS["w"])
    i3 = int(np.argmin(abs(xe - 3 * s)))
    check("exceed_from_hist: P(w' > 3 sigma) within 10 % of 1.35e-3", abs(pe[i3] / 1.3499e-3 - 1) < 0.1, f"{pe[i3]:.3e} at x = {xe[i3]:.2f}")
    check("exceed_from_hist: starts at 1, non-increasing", abs(pe[0] - 1) < 1e-12 and np.all(np.diff(pe) <= 1e-15))
    check("gauss_exceed: P(>0) = 0.5, P(>3s) = 1.3499e-3", abs(gauss_exceed([0.0], 2.0)[0] - 0.5) < 1e-12 and abs(gauss_exceed([3 * s], s)[0] / 1.3499e-3 - 1) < 1e-4)


# ----------------------------------------------------------------------------------------------- 3. box filter
def test_box_filter():
    section("gz_filter_compare.box_average (= postproc_highres.box_average)")
    ny, nx, n = 160, 320, 8; X, Y = np.meshgrid(np.arange(nx) + 0.5, np.arange(ny) + 0.5)
    w = rng.normal(size=(ny, nx))
    inline = w.reshape(ny // n, n, nx // n, n).mean(axis=(1, 3))
    check("box_average == the inline expression of the scripts (exact)", np.array_equal(box_average(w, n), inline) and np.array_equal(PH.box_average(w, n), inline))
    check("constant in -> constant out", np.allclose(box_average(np.full((ny, nx), 2.5), n), 2.5, atol=1e-12))
    ramp = 0.01 * X + 0.002 * Y; xb = (np.arange(nx // n) + 0.5) * n; yb = (np.arange(ny // n) + 0.5) * n
    Xb, Yb = np.meshgrid(xb, yb)
    check("linear ramp -> ramp at the box centres (exact)", np.allclose(box_average(ramp, n), 0.01 * Xb + 0.002 * Yb, atol=1e-10))
    s1 = np.sin(2 * np.pi * X / n)
    check("sine with wavelength = Delta -> 0", abs(box_average(s1, n)).max() < 1e-10)
    s4 = np.sin(2 * np.pi * X / (4 * n)); k = 2 * np.pi / (4 * n)
    fac = (np.sin(k * n / 2) / (n * np.sin(k / 2))) ** 2            # discrete n-point box: Dirichlet-kernel factor on the variance
    r = box_average(s4, n).var() / s4.var()
    check("sine with wavelength 4 Delta: variance x Dirichlet factor", abs(r / fac - 1) < 1e-6, f"{r:.5f} vs {fac:.5f}")
    check("mean preserved", abs(box_average(w, n).mean() - w.mean()) < 1e-12)
    w2 = rng.normal(size=(ny + 3, nx + 5))
    check("ragged edge: trailing cells dropped, shape (ny//n, nx//n)", box_average(w2, n).shape == (ny // n, nx // n))
    st = filter_stats([rng.normal(size=(50, 50)) * 0.7 for _ in range(10)])
    check("filter stats: std of N(0,0.7) within 2 %", abs(st["std"] / 0.7 - 1) < 0.02, f"{st['std']:.4f}")


# ----------------------------------------------------------------------------------------------- 4. time series tools
def test_timeseries():
    section("postproc_highres: autocorr / integral_scale / block_stats / h_*")
    n = 2 ** 14; dt = 1.0; t = np.arange(n) * dt; om = 2 * np.pi / 64
    rho = PH.autocorr(np.sin(om * t), 200); lag = np.arange(201)
    check("autocorr of a sine = cos(omega tau) (1 - tau/n) within 0.02", abs(rho - np.cos(om * lag) * (1 - lag / n)).max() < 0.02)
    rho = PH.autocorr(rng.normal(size=n), 200)
    check("autocorr of white noise: rho(0) = 1, |rho(tau>0)| < 4/sqrt(n)", abs(rho[0] - 1) < 1e-12 and abs(rho[1:]).max() < 4 / np.sqrt(n))
    phi = 0.9; r = phi ** np.arange(0, 2000); T = PH.integral_scale(r, dt)
    check("integral_scale of rho = phi^tau: (1/(1-phi) - 1/2) dt", abs(T - (1 / (1 - phi) - 0.5) * dt) < 1e-6, f"{T:.4f}")
    r2 = np.cos(np.arange(100) * np.pi / 50)
    check("integral_scale stops at the first zero crossing", abs(PH.integral_scale(r2, dt) - (np.sum(r2[:25]) - 0.5) * dt) < 1e-12)
    N = 2 ** 18; g = rng.normal(size=N); sd, sk, nb = PH.block_stats(g, 16)
    check("block_stats: std of 16-means = 1/4 (5 %), skew 0 +- 0.06", abs(sd * 4 - 1) < 0.05 and abs(sk) < 0.06 and nb == N // 16, f"std {sd:.4f} skew {sk:+.3f} n {nb}")
    z = np.arange(10, 2000, 20.0); zi = 1050.0
    wth = np.where(z <= zi, 0.2 * (1 - 1.2 * z / zi), -0.04 * (1 - (z - zi) / 300)); wth[z > zi + 300] = 0
    h = 190.0; tau = 0.08 * np.clip(1 - z / h, 0, None)
    th = 300 + 3 / (1 + np.exp(-(z - 1050) / 25))
    r = PH.depths(dict(zh=z, th=th, uw=tau, vw=np.zeros_like(tau), wth=wth))        # the three boundary-layer depths of the WRF post-processing
    check("depths/h_minflux: V-shaped heat flux -> z_i within one cell", abs(r["h_minflux"] - zi) <= 20, f"{r['h_minflux']:.0f}")
    check("depths/h_stress: linear stress to 0 at h -> h within 2 %", abs(r["h_stress"] / h - 1) < 0.02, f"{r['h_stress']:.1f}")
    check("depths/h_gradient: max gradient at the inversion within 2 cells", abs(r["h_grad"] - 1050) <= 40, f"{r['h_grad']:.0f}")


# ----------------------------------------------------------------------------------------------- 5. budget and the profile reader
def test_budget():
    section("tke_budget")
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ok = TB.selftest()
    check("tke_budget --selftest (analytic P_s, P_b, T_t, T_p)", ok)
    src = inspect.getsource(TB.budget)
    check("budget() reads buoyancy from <w'theta'> (col 10), not the dry-run-broken <w'theta_v'> (col 21)", 'CV["wth"]' in src and 'CV["wthv"]' not in src)
    # synthetic profile logs on the ERF staggered layout: known w'theta' profile, garbage in column 21, a duplicated time block
    n, dz = 50, 20.0; z = np.arange(n + 1) * dz; zc = z[:-1] + dz / 2; times = [3600.0, 3660.0, 3660.0, 3720.0]
    with tempfile.TemporaryDirectory() as d:
        wth = 0.2 * (1 - 1.2 * z / 1000.0); wth[0] = 0.0; wth[z > 1000] = 0.0
        hfx = np.zeros(n + 1); hfx[0] = 0.2; hfx[1] = 0.05
        rows_c, rows_m, rows_s = [], [], []
        for it, t in enumerate(times):
            garbage = -0.3 * np.ones(n + 1) if it != 1 else 7.0 * np.ones(n + 1)      # column 21: wrong on purpose
            for k in range(n + 1):
                c = np.zeros(22); c[0], c[1] = t, z[k]; c[2] = c[5] = 0.4; c[7] = 0.6 * np.sin(np.pi * min(z[k], 1000) / 1000) ** 2
                c[10] = wth[k]; c[14] = 0.01 * z[k] * (1 - z[k] / 2000); c[17] = -0.2 * np.sin(2 * np.pi * z[k] / 2000); c[21] = garbage[k]
                rows_c.append(c)
                m = np.zeros(16); m[0], m[1] = t, z[k]; m[2] = 5.0; m[3] = 0.0; m[5] = 1.1; m[6] = 300.0; rows_m.append(m)
                s = np.zeros(12); s[0], s[1] = t, z[k]; s[8] = hfx[k]; s[4] = -0.25 if k == 0 else 0.0; rows_s.append(s)
        for name, rows in (("covar_profiles.dat", rows_c), ("mean_profiles.dat", rows_m), ("sfs_profiles.dat", rows_s)):
            np.savetxt(os.path.join(d, name), np.array(rows))
        t_, z_, A = TB.load(os.path.join(d, "covar_profiles.dat"))
        # the two 3660 s blocks differ in column 21 (7.0 first, -0.3 second): the reader must keep the second (a restart rewrites the time)
        check("load(): duplicated time keeps the last block, z sorted", t_.tolist() == [3600.0, 3660.0, 3720.0] and np.allclose(z_, z) and abs(A[1, 0, 21] + 0.3) < 1e-12)
        r = TB.budget(d, (3600, 3720))
        Pb = np.array(r["terms"]["P_b"]); exp = TB.G / TB.THETA_REF * 0.5 * (wth[:-1] + wth[1:])
        check("budget(): P_b = (g/theta) <w'theta'> at cell centres, immune to column 21", np.allclose(Pb, exp, atol=1e-12))
        check("budget(): surface heat flux total = SGS(0) + resolved(0) = 0.20", abs(r["checks"]["surface_heat_flux_total"] - 0.2) < 1e-12 and r["checks"]["surface_heat_flux_resolved"] == 0.0)
        check("budget(): z_i = height of minimum total heat flux (1000 m here)", abs(r["zi"] - 1000.0) < 1e-9, f"{r['zi']:.0f}")
        check("budget(): u* from total surface stress = 0.5", abs(r["checks"]["ustar_from_total_stress"] - 0.5) < 1e-12)
        check("budget(): shear production 0 for a uniform wind", abs(np.array(r["terms"]["P_s"])).max() < 1e-12)
        Tt = np.array(r["terms"]["T_t"]); exp_t = -0.5 * np.diff(0.01 * z * (1 - z / 2000)) / dz
        check("budget(): T_t = -1/2 d<u_i u_i w>/dz from faces (exact)", np.allclose(Tt, exp_t, atol=1e-12))


# ----------------------------------------------------------------------------------------------- 6. real data (optional)
def test_real():
    section("real data")
    from erf_planes import Plotfile
    d = os.path.join(HERE, "..", "cases", "gz_cbl_dx100")
    pls = sorted(p for p in os.listdir(d) if p.startswith("plt")) if os.path.isdir(d) else []
    if pls:
        pf = Plotfile(os.path.join(d, pls[-1])); w = pf.plane("z_velocity", 26); P = spectrum2d(w)
        check("gz_cbl_dx100: Parseval on a real w plane", abs(P.sum() - w.var()) < 1e-9 * w.var(), f"var {w.var():.4f}")
        ny, nx = w.shape; ky = (np.fft.fftfreq(ny) * ny).round().astype(int); kx = (np.fft.fftfreq(nx) * nx).round().astype(int)
        KX, KY = np.meshgrid(kx, ky); kr = np.hypot(KX, KY); ring = np.array([P[(kr >= r - 0.5) & (kr < r + 0.5)].sum() for r in range(1, 40)])
        kpk = int(np.argmax(ring)) + 1
        check("gz_cbl_dx100: spectral peak at a few cycles per 16 km domain (convective rolls/cells ~ 1.5 z_i)", 1 <= kpk <= 12, f"peak at {kpk} cycles/domain = {16000 / kpk:.0f} m")
    else:
        print("  (gz_cbl_dx100 plotfiles not found, skipped)")
    import glob, netCDF4 as nc
    # the tiling diagnosis on the real runs: 14 x 16 ranks on 320 cells = 20-cell patches in y, all seeded alike
    for run, tiled in (("cbl_sd25", True), ("cbl_sd25_fix", False)):
        f = glob.glob(f"/work2/11859/hsna/stampede3/HowlandLab/WRF-4.6.1/test/abl_cases/{run}/planes_d01_*")
        if not f:
            print(f"  ({run} planes not found, skipped)"); continue
        ds = nc.Dataset(f[0]); ds.set_auto_mask(False); k = 25
        w = 0.5 * (np.asarray(ds["W"][-1, k], float) + np.asarray(ds["W"][-1, k + 1], float))
        per = bool(np.allclose(w, np.roll(w, 20, axis=0))); sh = line_share(spectrum2d(w), 16, axis=0, include_zero=True)
        what = "tiled run: w at 4 h still exactly periodic every 20 cells in y, all energy on ky = n*16" if tiled else "re-run (14 x 15): not periodic, energy spread over all ky"
        check(f"WRF {run}: {what}", (per and sh > 0.999) if tiled else (not per and sh < 0.5), f"periodic {per}, share(16, incl. 0) {sh:.3f}")
    e = os.path.join(HERE, "..", "cases_pvc", "ens_cbl25_s1")
    if os.path.isfile(os.path.join(e, "covar_profiles.dat")):
        r = TB.budget(e, (10800, 14400)); c = r["checks"]
        check("ens_cbl25_s1: surface heat flux closes to 0.200 K m/s", abs(c["surface_heat_flux_total"] - 0.2) < 1e-3, f"{c['surface_heat_flux_total']:.4f}")
        check("ens_cbl25_s1: z_i 1000-1100 m, w* 1.8-2.0 m/s", 1000 <= r["zi"] <= 1100 and 1.8 <= r["wstar"] <= 2.0, f"zi {r['zi']:.0f} w* {r['wstar']:.2f}")
        check("ens_cbl25_s1: w'w', <u_i u_i w>, p'w' vanish on the ground face", max(c["ground_face_zero"].values()) == 0.0)
    else:
        print("  (ens_cbl25_s1 profiles not found, skipped)")


if __name__ == "__main__":
    test_spectra(); test_moments(); test_box_filter(); test_timeseries(); test_budget()
    if "--real" in sys.argv:
        test_real()
    nfail = sum(1 for _, ok in RES if not ok)
    print(f"\n{len(RES) - nfail} passed, {nfail} failed")
    sys.exit(1 if nfail else 0)
