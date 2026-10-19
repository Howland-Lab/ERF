# Changelog (branch `hsna`, HowlandLab/stampede3)

Newest first. Record every change to the build/run environment here, with the symptom that prompted it.

## 2026-10-19
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
