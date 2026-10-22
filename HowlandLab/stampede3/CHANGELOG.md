# Changelog (branch `hsna`, HowlandLab/stampede3)

Newest first. Record every change to the build/run environment here, with the symptom that prompted it.

## 2026-10-22
- H100 works: `cases_gpu/cbl_sd25` ran 2400 steps (600 s) in 228 s on 4 GPUs with `I_MPI_FABRICS=shm` (job 3616319, 21 Oct);
  the 4 h case would take ~1.5 h, about the PVC node's 1.9 h.
- Grey-zone PBL-scheme scan (`cases/make_gz_pbl_runs.py`, `run/submit_cpu_chain.slurm`): MRF and MRF+blending ran at 100/200/400 m;
  MYNN2.5 and MYNN-EDMF stop at step 1–2 with `RhoTheta is negative ... -nan` in the lowest cell, although the same deck runs with
  MRF and MYNN2.5 runs at 1 km.  Debug deck `gzpbl_mynn25_dx400_fpe` (trap on, 10 steps) submitted to locate the first invalid op.
- `cases/make_ar_runs.py`: aspect-ratio test (dz 50 and 100 m at dx 100 m) for the VLES scan, after the remark that the LES SGS model
  is isotropic so dx/dz matters.  `cases/cbl_woa_dx100`: first terrain-fitted run (`prob.custom_terrain_type = WoA`).
- `analysis/tke_budget.py`: resolved TKE budget from the profile logs, with a self-test (`--selftest`) and run-time checks
  (Mike, 22 Oct: analysis code written with AI help must be verified).  Found that `covar_profiles` column 21 ⟨w'θ_v'⟩ is wrong for
  dry runs (use column 10), and that in the 25 m convective case the Deardorff dissipation accounts for well under half of the
  resolved-TKE loss; see `record/RUNS.md`.
- `analysis/verify_analysis.py` + `analysis/spectra.py`: known-answer tests (66 checks, one command) for the analysis code behind
  the 22 Oct deck and the results blog: spectra (sine -> one peak, Parseval, white noise flat, 16 copies -> comb, shifted copies
  -> no comb), Moments (Gaussian/exponential skewness, kurtosis, 3-sigma tails, streaming), merged/masked PDF bins, exceedance
  curves, the box filter of the fair test, autocorrelation/integral scale/block statistics, the three boundary-layer depths, and
  the budget (analytic terms, surface closure, duplicated-time blocks, the w'theta_v trap).  `--real` re-derives the tiling
  diagnosis on the WRF planes: the 14 x 16-rank run is still exactly periodic every 20 cells in y at 4 h.  No bug was found in
  the published numbers; `spectra.py` fixes the spectrum normalisation (1-D: sum E dk = var; 2-D: |FFT2|^2/N^2, sum = var).
  To make the functions importable, `planes_clean_fig.py` and `gz_filter_compare.py` (here) and `nwp_vs_les_tails.py`,
  `postproc_highres.py` (WRF tree) got `main()` guards and `merged_pdf()` / `box_average()`; their outputs were regenerated and
  are bit-identical to before.
- For the record (WRF side, not in this repo): WRF's `module_initialize_ideal.F` now takes an optional `WRF_SEED` environment
  variable that seeds the initial noise per rank and member, so WRF ensembles comparable to `ens_cbl25_s*` / `ens_sbl6_s*` were
  submitted on 22 Oct.

- `analysis/veer_profiles.py`: wind speed/direction (veer) profiles, WRF vs ERF ensemble means (Mike, 22 Oct). WRF ≈ ERF:
  convective surface-layer direction 255.0 ± 0.5 vs 255.6°, stable 236.3 ± 0.4 vs 237.5° with +10° of turning by 100 m and the
  LLJ at 166 (ERF) vs 184 m (WRF); both tend to 270° (geostrophic) at the top.
- `analysis/gz_scheme_analysis.py` fixed (planes stream has no XTIME → Times strings; WRF stress on cell levels) and run on the
  finished scan members (`record/results/gz_schemes.json`); WRF `mynn_edmf` ≡ `mynn_nomf` at 100–200 m because the hard-wired
  scale-aware Psig switches the mass flux off at Δ ≪ z_i (namelists differ, BL_MYNN_EDMF 1 vs 0 confirmed in namelist.output).
- `analysis/ar_test_analysis.py` and `analysis/terrain_first_look.py` prepared for the pending aspect-ratio and terrain runs
  (smoke-tested on gz_cbl_dx100).
- `record/draft_email_NCAR_3DPBL_20261022.md`: draft of the question to NCAR about a 3-D PBL port into ERF (not sent).
- TODO (22 Oct 2026, from the TKE-budget work): add an SGS-production column <tau_ij s_ij> to sfs_profiles.dat so the resolved-TKE
  budget can separate the transfer to the subgrid scales from numerical dissipation (needed for Smagorinsky, where no epsilon
  exists; a check of the equilibrium assumption for Deardorff).  Where: Source/IO/ERF_Write1DProfiles.cpp (derive_diag_profiles,
  next to h_avg_sgsdiss) using the stored Tau11..Tau33 MultiFabs and the strain from the velocity gradients already formed in
  Source/Diffusion/ERF_ComputeStrain.cpp.  Needs a rebuild (build job on skx/spr) and a check that the column integrates to the
  same mixed-layer value as -epsilon_sgs in the Deardorff convective case (tke_budget.py gives -0.16 w*^3/zi).

## 2026-10-20
- 1 km MYNN2.5 runs in ERF (`pbl_type = MYNN25`, `les_type = Smagorinsky2D`, `Cs = 0.25`), mirroring WRF's 1 km runs.
  MYNN25 trips the floating-point trap at step 1 (`ComputeDiffusivityMYNN25`), so these decks set `amrex.fpe_trap_invalid = 0`
  as ERF's own `inputs_GABLS1_mynn25` does; profiles checked for NaN (none).
- `analysis/erf_rans_compare.py` (ERF vs WRF 1 km vs LES), `analysis/planes_clean_fig.py` (PDFs without the invalid WRF runs,
  merged tail bins, sample count in the legend).

## 2026-10-19
- `record/`: run table (`RUNS.md`), exact inputs of all 20 runs, and analysis summary JSON. Outputs stay on Stampede3.
- Created branch `hsna` from Howland-Lab/ERF `development` (`3a542f5d4`) and moved all Stampede3 material into
  `HowlandLab/stampede3/` (env, build, run, cases, analysis). No ERF source changes.
- H100: `I_MPI_FABRICS=shm` in `env/h100.sh`. The job of 17 Oct aborted in `MPI_Init` with
  `OFI fi_getinfo() failed (No data available)`, exit 143.
- Paths generalised: every script reads `ERF_HOME`; build trees in `_build/<target>`, runs in `runs/` (git-ignored).
- Documented the comparison with Ethan's CPU install guide and the ekat SSH-submodule workaround.

## 2026-10-15/16
- `run/submit_pvc_multi.slurm`: several runs per pvc job (pvc allows 2 running jobs per user).
- `cases/make_variants.py`: grey-zone scan (Δx 50–400 m, 16 km domain) and seed ensembles via `erf.random_seed`.
- Results: PVC convective ensemble (4 seeds), CPU stable ensemble (5 seeds), grey-zone scan; analysis in `analysis/`.

## 2026-10-14
- PVC (Intel Max 1550, SYCL) build and run working.
  - `intel/25.1` instead of the default `intel/24.0` (AMReX needs `sycl::ext::intel::math::jn/yn`).
  - `libumf` path added to `LD_LIBRARY_PATH` (missing from TACC's intel/25.1 module; without it the SYCL runtime
    sees no device).
  - AOT build for `pvc` must run on a pvc node (`ocloc`); `build_pvc.slurm` builds and smoke-tests.
  - `ZE_FLAT_DEVICE_HIERARCHY=FLAT`, 8 ranks per node, boxes 80 x 160 for the 320 x 320 convective grid.
- First PVC production run: convective 25 m, 4 h simulated in 1.9 h on one node (CPU: 2.5 h on four nodes).

## 2026-10-13
- H100: `env/h100.sh` loads `impi/21.9` explicitly (+ fallback PATH/LD_LIBRARY_PATH); `submit_h100.slurm`
  checks MPI and libraries before running. The jobs of 11 and 12 Oct died after 3 s with `mpirun: command not found`.

## 2026-10-08
- H100 (CUDA 12.8, sm_90) build on the login node with gcc 13.2 + Intel MPI 21.9.

## 2026-10-07
- CPU build (Intel oneAPI 2024.0 + Intel MPI 2021.11, pinned NetCDF/HDF5, FFT) of upstream erf-model/ERF
  `e50b1611b`; login-node -O1 fallback for files that exceed the 8 GB vmem limit.
- Input decks mirroring the WRF-LES namelists (GABLS1 6.25/3.125 m, Shin & Dudhia convective 25 m); WRF/ERF default
  traps documented (`Pr_t`, `Ri_crit`, `Ce_wall`).
