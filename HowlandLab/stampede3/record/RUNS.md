# Run record (ERF, Oct 2026)

Every ERF run done so far, with the exact `inputs` / `input_sounding` used (`inputs/<folder>__<run>/`) and the summary
numbers from the analysis (`results/*.json`). **Outputs (plotfiles, profiles, logs; ~155 GB) are not in git.** They are
on Stampede3 under `/work2/11859/hsna/stampede3/HowlandLab/ERF_runs/<folder>/<run>/`.

All successful runs used ERF `26.09-538-ge50b1611bc98` (upstream erf-model/ERF, as printed in each `log.erf`).
Wall time = ERF "Total Time". PVC runs override the boxes on the command line: `amr.max_grid_size_x=80 amr.max_grid_size_y=160`.

| run (folder/run) | what it is | grid | platform, MPI ranks | simulated | wall | finished | status |
|---|---|---|---|---|---|---|---|
| cases/sbl_dx6p25 | GABLS1 stable, mirror of WRF-LES | 64³, 6.25 m | CPU SPR, 32 | 9 h | 4.7 h | 10-08 | ok |
| cases/sbl_dx3p125 | GABLS1 stable, mirror of WRF-LES | 128³, 3.125 m | CPU SPR, 64 | 9 h | 32.4 h | 10-09 | ok |
| cases/cbl_sd25 | Shin & Dudhia convective, mirror of WRF-LES | 320²×100, 25 m, 8 km | CPU SPR, 256 (4 nodes) | 4 h | 2.5 h | 10-08 | ok |
| cases_gpu/cbl_sd25 | same, NVIDIA H100 | 320²×100 | H100, 4 | — | — | 10-11 | failed: `mpirun: command not found` (exit 127) |
| cases_gpu/cbl_sd25_full | same, NVIDIA H100 | 320²×100 | H100, 4 | — | — | 10-17 | failed: `OFI fi_getinfo() failed` (exit 143) → `I_MPI_FABRICS=shm` |
| cases_pvc/cbl_sd25_test | PVC smoke test | 320²×100 | PVC, 8 | — | — | 10-14 | failed: no SYCL device (libumf) → fixed |
| cases_pvc/cbl_sd25 | same convective case, Intel PVC | 320²×100 | PVC 1 node, 8 tiles | 4 h | 1.9 h | 10-15 | ok, = CPU (u* −0.3 %) |
| cases_pvc/ens_cbl25_s1..s4 | convective ensemble, `erf.random_seed` 1–4 | 320²×100 | PVC 1 node each | 4 h | 1.9 h each | 10-15 | ok |
| cases/ens_sbl6_s1..s5 | stable ensemble, `erf.random_seed` 1–5 | 64³, 6.25 m | CPU SPR, 32 each | 9 h | 4.7 h each | 10-16 | ok |
| cases/gz_cbl_dx050 | grey-zone scan, Δ = 50 m, 16 km domain | 320²×100 | CPU SPR, 256 | 4 h | 2.4 h | 10-16 | ok |
| cases/gz_cbl_dx100 | grey-zone scan, Δ = 100 m | 160²×100 | CPU SPR, 64 | 4 h | 2.3 h | 10-16 | ok |
| cases/gz_cbl_dx200 | grey-zone scan, Δ = 200 m | 80²×100 | CPU SPR, 16 | 4 h | 1.9 h | 10-16 | ok |
| cases/gz_cbl_dx400 | grey-zone scan, Δ = 400 m | 40²×100 | CPU SPR, 16 | 4 h | 0.6 h | 10-16 | ok |
| cases/rans_sbl_1km | GABLS1 at 1 km, MYNN2.5 + Smagorinsky2D (mirror of WRF sbl_rans) | 100²×64, 1 km | CPU SPR, 25 | 9 h | 0.14 h | 10-20 | ok (`amrex.fpe_trap_invalid = 0` needed) |
| cases/rans_cbl_1km | Shin & Dudhia at 1 km, MYNN2.5 + Smagorinsky2D (mirror of WRF cbl_rans) | 100²×120, 1 km | CPU SPR, 25 | 4 h | 0.18 h | 10-20 | ok |
| cases_gpu/cbl_sd25 (3rd try) | same, NVIDIA H100 with `I_MPI_FABRICS=shm` | 320²×100 | H100, 4 GPUs | 600 s (test, 2400 steps) | 228 s | 10-21 | ok → ~1.5 h for the 4 h case |
| cases/rans_cbl_1km_edmf | Shin & Dudhia at 1 km, MYNN-EDMF (mass flux on, as WRF's default `bl_mynn_edmf=1`) | 100²×120, 1 km | CPU SPR, 25 | 4 h | 0.16 h | 10-21 | ok |
| cases/gzpbl_mrf_dx{100,200,400} | grey-zone PBL-scheme scan: MRF (K-profile) + Smagorinsky2D horizontally | 160²/80²/40² × 100, dz 20 m | CPU SPR, 64/64/16 | 4 h | 0.54 / 0.09 / 0.04 h | 10-21 | ok |
| cases/gzpbl_mrfblend_dx{100,200,400} | same, MRF + Boutle 2014 blending `pbl_blend_length = 750` | as above | as above | 4 h | 0.54 / 0.09 / 0.04 h | 10-21 | ok |
| cases/gzpbl_mynn25_dx{100,200,400} | same, MYNN2.5 (`fpe_trap_invalid = 0`) | as above | as above | — | — | 10-21 | failed at step 1–2: `RhoTheta is negative at (i,0,0) -nan` (MRF fine on the same deck; 1 km MYNN fine) |
| cases/gzpbl_mynnedmf_dx{100,200,400} | same, MYNN-EDMF | as above | as above | — | — | 10-21 | failed, same message |
| cases/gzpbl_mynn25_dx400_fpe | debug copy, 10 steps with `fpe_trap_invalid = 1` to get the backtrace of the first invalid op | 40²×100 | CPU SPR, 16 | 40 s | — | submitted 10-22 | pending |
| cases/gz_cbl_dx100_dz{050,100} | grid aspect-ratio test of the VLES scan: dz 50 and 100 m at dx 100 m (ratio 2 and 1, vs 5 in gz_cbl_dx100) | 160²×40, 160²×20 | CPU SPR, 64 | 4 h | — | submitted 10-22 | pending |
| cases/cbl_woa_dx100 | first terrain test: gz_cbl_dx100 over a Witch-of-Agnesi ridge in x (hmax 200 m, L 1 km), `StaticFittedMesh` | 160²×100 | CPU SPR, 64 | 4 h | — | submitted 10-22 | pending |

## Main results (details in `results/`)

- **ERF vs WRF-LES** (`summary.json`). Stable 3.125 m: u* 0.278 vs 0.277 m/s, h 191 vs 195 m. Convective 25 m:
  mixed-layer depth 1040 vs 1052 m (WRF re-run) and the same θ and heat-flux profiles. WRF's extra wind/TKE bump at 600–900 m came
  from WRF's identical per-rank initial noise (exact tiling with 14 × 16 ranks), and it disappears with 14 × 15 ranks.
- **CPU vs PVC**: the same answer. PVC on 1 node takes 1.9 h; CPU on 4 nodes takes 2.5 h.
- **Ensembles** (`week2_summary.json`):
  - Means vary by 1–3 % between seeds.
  - Tail probabilities vary more: P(w′ > 3σ) ±31 % (convective, 6 members) and ±13 % (stable, 5 members).
  - In the stable case θ still drifts by 0.2 K/h at 8–9 h, against a 0.02 K spread between members, so it is not
    stationary. The wind is steady.
- **1 km MYNN in ERF** (`rans1km_summary.json`). Stable: close to LES (h 197 vs 191–195 m, u* −6 %). Convective: zi 1050 m = LES;
  WRF's default `bl_pbl_physics=5` has the EDMF mass flux ON (`bl_mynn_edmf=1`) and gives zi 1150 m; with it off WRF also gives 1051 m.
  u* is 8–20 % too high in every 1 km run (momentum over-mixing).
- **Grey zone** (`week2_summary.json`, `gz_fairtest.json`). Compared with the 25 m LES box-averaged to the same Δ:
  - Coarse grids keep the means.
  - Coarse grids overestimate σw by 11–18 %.
  - u* rises from 0.517 (25 m) to 0.541 (400 m).
- **Resolved TKE budget** (`tke_budget.json`, `analysis/tke_budget.py`, from the covar/sfs profile logs; terms at cell centres on the
  staggered grid, self-test against analytic profiles, surface-flux and ground-face checks).  Convective 25 m, 0.2–0.8 z_i, in
  units of w_*³/z_i (4 members): buoyancy +0.42, turbulent transport +0.07, pressure transport −0.12, shear ≈ 0; the Deardorff
  dissipation is only −0.16 while the residual is −0.37, i.e. more than half of the loss of resolved TKE is numerical (the
  5th/3rd-order upwind advection) plus SGS transport.  Stable 6.25 m (Smagorinsky, 5 members, 20–100 m): shear production
  1.8e-3 m²/s³, of which the mean-shear SGS transfer takes 9 %; the rest is in the residual (fluctuating-strain SGS transfer and
  numerics cannot be separated without a ⟨τ_ij s_ij⟩ output).  Note: `covar_profiles` column 21 ⟨w'θ_v'⟩ is unusable for dry runs
  in this ERF version (negative at every height); the budget uses column 10 ⟨w'θ'⟩.
