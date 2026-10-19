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
- **Grey zone** (`week2_summary.json`, `gz_fairtest.json`). Compared with the 25 m LES box-averaged to the same Δ:
  - Coarse grids keep the means.
  - Coarse grids overestimate σw by 11–18 %.
  - u* rises from 0.517 (25 m) to 0.541 (400 m).
