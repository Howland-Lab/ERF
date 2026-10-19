# ERF on TACC Stampede3: CPU, NVIDIA H100 and Intel PVC (branch `hsna`)

Maintainer: Hyunsoo Na (Howland Lab, MIT). Started 7 Oct 2026; see `CHANGELOG.md` for every change.

**The ERF source is not modified on this branch.** Everything here is build/run environment, input decks and
analysis scripts for Stampede3, collected in this folder. Build trees go to `_build/` and run directories to
`runs/`, both git-ignored.

```
HowlandLab/stampede3/
  env/       cpu.sh  h100.sh  pvc.sh      module environments (source one of them)
  build/     build_cpu.sh  build_cpu.slurm  build_h100.sh  build_pvc.slurm
  run/       submit_cpu.slurm  submit_h100.slurm  submit_pvc.slurm  submit_pvc_multi.slurm  smoke_test_login.sh
  cases/     sbl_dx6p25  sbl_dx3p125  cbl_sd25  (decks that mirror our WRF-LES runs)  make_variants.py
  analysis/  scripts used for the ERF-vs-WRF comparison, ensembles and grey-zone scan (paths at the top are site-specific)
```

## Quick start

```bash
git clone https://github.com/Howland-Lab/ERF.git ERF && cd ERF && git checkout hsna
# ekat lists its own submodules with SSH URLs; without a GitHub SSH key fetch them over HTTPS (from Ethan's guide):
git -c url."https://github.com/".insteadOf="git@github.com:" submodule update --init --recursive
git submodule status --recursive        # no line may start with "-"
export ERF_HOME=$PWD                    # every script and slurm file reads ERF_HOME

# CPU:   bash HowlandLab/stampede3/build/build_cpu.sh -j 8     (login node; or sbatch build/build_cpu.slurm)
# H100:  bash HowlandLab/stampede3/build/build_h100.sh -j 4    (login node, nvcc needs no GPU)
# PVC:   sbatch HowlandLab/stampede3/build/build_pvc.slurm     (must build ON a pvc node; includes a 10-min test run)

cd somewhere/case && cp $ERF_HOME/HowlandLab/stampede3/cases/cbl_sd25/* .
sbatch -J cbl25 -t 04:00:00 $ERF_HOME/HowlandLab/stampede3/run/submit_pvc.slurm "$PWD"
```
Write large output to `$SCRATCH` (plotfiles grow quickly; `$WORK` has a quota).

## The three targets side by side

| | CPU (SPR) | NVIDIA H100 | Intel Max 1550 (PVC) |
|---|---|---|---|
| Stampede3 partition | `spr` (112 cores) | `h100` (4 GPUs, 96 cores) | `pvc` (4 GPUs = 8 tiles, 96 cores) |
| charge (SU per node-hour) | 2 | 4 | 3 |
| compiler | Intel oneAPI 2024.0 `mpiicpx` (default `intel/24.0`) | gcc 13.2 + CUDA 12.8 (`nvcc`) | **Intel oneAPI 2025.1** `icpx -fsycl` (`intel/25.1`) |
| MPI | Intel MPI 2021.11 | Intel MPI 2021.9 (gcc wrappers) | Intel MPI 2021.15 |
| AMReX backend / CMake | none, `-DERF_ENABLE_MPI=ON` | `-DERF_ENABLE_CUDA=ON -DAMReX_CUDA_ARCH=9.0` | `-DERF_ENABLE_SYCL=ON -DAMReX_SYCL_AOT=ON -DAMReX_INTEL_ARCH=pvc` |
| NetCDF / FFT | NetCDF+HDF5 (pinned), FFT on | off / on | off / off |
| where to build | login (with -O1 fallback) or compute node | login node | **on a pvc node** (`ocloc` exists only there) |
| ranks per node | 32–64 | 4 (1 per GPU) | 8 (1 per tile, `ZE_FLAT_DEVICE_HIERARCHY=FLAT`) |
| boxes for 320 x 320 x 100 | 20 x 20 (256 boxes) | 160 x 160 (4 boxes) | 80 x 160 (8 boxes) |
| runtime environment quirks | pinned NetCDF/HDF5 modules | `I_MPI_FABRICS=shm`, impi loaded explicitly | **`libumf` on `LD_LIBRARY_PATH`** |
| status (Oct 2026) | production runs done | built; first runs failed in MPI init, fixed, not yet re-run | production runs done |

Convective 25 m case (320 x 320 x 100, 4 h, 57 600 steps): **CPU 4 nodes 2.5 h (~20 SU) vs PVC 1 node 1.9 h (~6 SU)**,
same answer (u* within 0.3 %, zi identical; `analysis/`).

## What had to be fixed, and why (symptom → cause → fix)

1. **CPU build dies on the login node** with `LLVM ERROR: out of memory` on four AMReX/ERF files.
   Login nodes cap virtual memory at 8 GB; icpx -O3 exceeds it on the biggest files.
   `build_cpu.sh` recompiles exactly those files at -O1 and resumes. Building on a compute node avoids it entirely
   (Ethan's guide: 6 min at -j40 on one skx node), which is the better default.
2. **H100 run: `mpirun: command not found`, exit 127 after 3 s** (13 Oct). After `module reset; module load gcc`
   the h100 nodes had no active Intel MPI, and the env file hid module messages with `> /dev/null`.
   `env/h100.sh` loads `impi/21.9` explicitly (with a fallback to the exact MPI the binary is linked against), and
   `submit_h100.slurm` prints the modules and aborts early if MPI or a library is missing.
3. **H100 run: `OFI fi_getinfo() failed ... No data available` in MPI_Init, exit 143** (17 Oct). Intel MPI 2021.9
   could not open its network fabric on the h100 node. The job is single-node, so `I_MPI_FABRICS=shm` (shared memory
   only) is enough. Note: Ethan saw the same message from a different cause (conda's Intel Python libfabric on
   `LD_LIBRARY_PATH`); our shells have no conda, but check that first if it appears elsewhere.
4. **PVC build fails after 1 min: `no member named 'jn' in namespace 'sycl::ext::intel::math'`** (14 Oct).
   The default `intel/24.0` (oneAPI 2024.0) is older than this AMReX expects. Use `intel/25.1` (`intel/26.0` also
   compiles).
5. **PVC run: `No device of requested type 'info::device_type::gpu' available`** on every rank (14 Oct), even
   `sycl-ls` lists 0 platforms with `intel/25.1`. TACC's module does not put `libumf.so.0` (Unified Memory
   Framework) on `LD_LIBRARY_PATH`, so the SYCL runtime cannot load its Level Zero / OpenCL adapters
   (`SYCL_UR_TRACE=1 sycl-ls` shows `libumf.so.0: cannot open shared object file`). Fix in `env/pvc.sh`:
   `export LD_LIBRARY_PATH=/scratch/projects/compilers/intel25.1/umf/0.10/lib:$LD_LIBRARY_PATH`.
   `submit_pvc.slurm` aborts early if no GPU is visible.
6. **PVC ahead-of-time compilation needs `ocloc`,** which is installed only on pvc compute nodes, so the PVC build
   is a batch job (`build_pvc.slurm`, ~40 min) that also runs a 10-minute smoke test.
7. **Recursive clone stops at `Submodules/ekat`**: ekat's own submodules use SSH URLs. Use the HTTPS rewrite in
   the quick start (unless a GitHub SSH key is set up).

## Input decks: how they mirror the WRF-LES namelists

The decks in `cases/` reproduce our WRF ideal LES runs (GABLS1 stable case, Shin & Dudhia 2016 convective case).
ERF defaults that differ from WRF and must be set explicitly:

| WRF | ERF |
|---|---|
| `km_opt = 3, c_s = 0.15` | `erf.les_type = Smagorinsky`, `erf.Cs = 0.15`, **`erf.Pr_t = 1/3` (ERF default 1.0)**, `erf.use_Ri_correction = true`, **`erf.Ri_crit = 1/3` (ERF default 0.25)** |
| `km_opt = 2, c_k = 0.10` | `erf.les_type = Deardorff`, `erf.Ck = 0.10`, `erf.Ce = 0.93`, **`erf.Ce_wall = 0`** |
| `time_step_sound = 6`, 5th/3rd-order advection | `erf.fixed_mri_dt_ratio = 6`, `Upwind_5th` / `Upwind_3rd` |
| `sf_sfclay_physics = 1`, `spec_sf_heatrate` / `tke_heat_flux` | `zlo.type = surface_layer`, `erf.most.surf_heating_rate` / `erf.most.surf_temp_flux` |
| `pert_coriolis`, geostrophic base wind | `erf.abl_driver_type = GeostrophicWind`, `erf.abl_geo_wind`, `erf.latitude`, `erf.rotational_time_period` |
| `damp_opt = 2, zdamp, dampcoef` | `erf.rayleigh_damp_{U,V,W,T} = true`, `erf.rayleigh_zdamp`, `erf.rayleigh_dampcoef` |
| random θ perturbation (ideal.exe) | `prob.T_0_Pert_Mag`, `prob.pert_ref_height`; **`erf.random_seed = N` gives independent ensemble members** |

`cases/make_variants.py` builds the grey-zone scan (convective, Δx = 50–400 m, 16 km domain) and the seed ensembles
(stable 6.25 m seeds 1–5, convective 25 m seeds 1–4) used in Oct 2026.

## Provenance of the Oct 2026 results

The CPU, PVC and ensemble results reported in Oct 2026 were produced with **upstream erf-model/ERF `e50b1611b`
(7 Oct 2026, AMReX 26.09-179)**, built with exactly these scripts. This branch starts from Howland-Lab/ERF
`development` at `3a542f5d4` (19 Oct 2026, AMReX 26.10). Rebuild from the branch before new production runs and
record the commit in each run's log (`git -C $ERF_HOME rev-parse HEAD`).

Validation: with the same grid and closure, ERF reproduces the WRF-LES mean profiles, fluxes and w distributions
(stable 3.125 m: u* within 0.5 %, h within 2 %; convective 25 m: mixed-layer depth 1040 vs 1072 m, same θ and
heat-flux profile; w skewness at mid-layer 0.70 vs 0.77 in the WRF re-run). The convective WRF run first showed a wind/TKE
bump at 600–900 m; that was a WRF artefact (ideal LES seeds every MPI rank identically, and with 14 x 16 ranks the
patch tiles the domain exactly), gone after re-running with 14 x 15 ranks. The CPU and PVC builds give the same answer
(u* −0.3 %, zi identical).

## Differences from Ethan's CPU guide (Oct 2026)

| | Ethan's guide | this branch |
|---|---|---|
| CPU compiler | gcc 13.2 + Intel MPI 21.9 ("default Intel does not work") | Intel oneAPI 2024.0 + Intel MPI 21.11, works (only the login-node OOM above) |
| NetCDF/HDF5 | `phdf5` (unpinned, 1.14.6), HDF5 off | pinned phdf5 1.14.4 / parallel-netcdf 4.9.2 / pnetcdf 1.12.3 (avoids the Oct 2026 TACC zlib crash) |
| build location | skx compute node, 6 min | login node with -O1 fallback, or compute node |
| executable | `cmake --install` → `install/bin/erf_exec` | `_build/<target>/Exec/erf_exec` |
| GPU | not covered | H100 and PVC (above) |
Either CPU recipe should give the same executable behaviour; for consistency in the lab, Ethan's gcc recipe can
replace `env/cpu.sh` if preferred.
