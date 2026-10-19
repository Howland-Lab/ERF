#!/bin/bash
# Build ERF for Intel PVC (SYCL, ahead-of-time for pvc) ON a pvc node -- ocloc, needed for AOT compilation, is only
# installed there -- then run a 10-minute smoke test of the convective case.  ~40 min.
#     export ERF_HOME=/path/to/ERF ; sbatch build_pvc.slurm
#SBATCH -J erf_build_pvc
#SBATCH -o erf_build_pvc.o%j
#SBATCH -p pvc
#SBATCH -N 1
#SBATCH --ntasks-per-node 8
#SBATCH -t 04:00:00
#SBATCH -A TG-ATM170028
: "${ERF_HOME:?export ERF_HOME=/path/to/this/ERF/clone before sbatch}"
source "$ERF_HOME/HowlandLab/stampede3/env/pvc.sh"
export HOME=${HOME:-/home1/$(id -u)}
echo "node: $(hostname)   started: $(date)"; module -t list 2>&1 | tr '\n' ' '; echo
echo "ocloc: $(command -v ocloc || echo MISSING)"; echo "SYCL Level Zero GPUs: $(sycl-ls 2>/dev/null | grep -c 'level_zero:gpu')"
AOT="-DAMReX_SYCL_AOT=ON -DAMReX_INTEL_ARCH=pvc"
command -v ocloc > /dev/null || { echo "no ocloc: JIT build (slow first start)"; AOT="-DAMReX_SYCL_AOT=OFF"; }
mkdir -p "$HL_BUILD"; rm -rf "$ERF_BUILD"
cmake -S "$ERF_HOME" -B "$ERF_BUILD" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx \
      -DERF_DIM=3 -DERF_ENABLE_MPI=ON -DERF_ENABLE_SYCL=ON $AOT -DERF_ENABLE_FFT=OFF -DERF_ENABLE_NETCDF=OFF \
      -DERF_ENABLE_TESTS=OFF -DERF_ENABLE_UNIT_TESTS=OFF -DERF_ENABLE_FCOMPARE=OFF -DERF_ENABLE_DOCUMENTATION=OFF \
      > "$HL_BUILD/configure_pvc.log" 2>&1
echo "configure exit $?"
cmake --build "$ERF_BUILD" -j 48 > "$HL_BUILD/build_pvc.log" 2>&1
echo "build exit $?   $(date)"; grep -E "error:" "$HL_BUILD/build_pvc.log" | sort -u | head -20
ls -la "$ERF_EXE" || exit 1
T=$HL_DIR/runs/pvc_smoke; mkdir -p "$T"; cp "$HL_DIR/cases/cbl_sd25/inputs" "$HL_DIR/cases/cbl_sd25/input_sounding" "$T/"; cd "$T"
ibrun "$ERF_EXE" inputs stop_time=600 erf.plot_per_1=300 erf.check_per=-1 amr.max_grid_size_x=80 amr.max_grid_size_y=160 > log.erf 2>&1
echo "smoke test exit $?   $(date)"; grep "Coarse STEP" log.erf | tail -1
