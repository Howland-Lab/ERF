#!/bin/bash
# Configure + build ERF for CPU (MPI, NetCDF, HDF5, FFT; no OpenMP) into _build/spr.   Usage: build_cpu.sh [-j N]
# On a login node icpx -O3 exceeds the 8 GB virtual-memory limit on a few AMReX/ERF files ("LLVM ERROR: out of memory");
# the loop below recompiles exactly those files at -O1 and resumes (needed 2 passes / 4 files in Oct 2026).
# Alternatively build on a compute node with build_cpu.slurm (no memory limit).
J=8; [ "$1" = "-j" ] && J=$2
source "$(dirname "$0")/../env/cpu.sh"
cmake -S "$ERF_HOME" -B "$ERF_BUILD" \
      -DCMAKE_CXX_COMPILER=mpiicpx -DCMAKE_C_COMPILER=mpiicx -DCMAKE_Fortran_COMPILER=mpiifx \
      -DCMAKE_BUILD_TYPE=Release -DERF_DIM=3 -DERF_ENABLE_MPI=ON -DERF_ENABLE_OPENMP=OFF \
      -DERF_ENABLE_NETCDF=ON -DERF_ENABLE_HDF5=ON -DNETCDF_DIR=$NETCDF_DIR -DERF_ENABLE_FFT=ON \
      -DERF_ENABLE_TESTS=OFF -DERF_ENABLE_UNIT_TESTS=OFF -DERF_ENABLE_FCOMPARE=ON \
      -DERF_ENABLE_DOCUMENTATION=OFF -DCMAKE_EXPORT_COMPILE_COMMANDS=ON > "$HL_BUILD/configure_spr.log" 2>&1
echo "configure exit $?  (log: $HL_BUILD/configure_spr.log)"
for it in 1 2 3 4 5 6; do
  cmake --build "$ERF_BUILD" -j $J > "$HL_BUILD/build_spr_pass$it.log" 2>&1
  [ -x "$ERF_EXE" ] && { echo "built: $ERF_EXE"; exit 0; }
  files=$(grep -B4 "out of memory" "$HL_BUILD/build_spr_pass$it.log" | grep -oE "[A-Za-z_0-9]+\.cpp\.o" | sort -u | sed 's/\.o$//')
  [ -z "$files" ] && { echo "build failed without OOM, see $HL_BUILD/build_spr_pass$it.log"; exit 1; }
  echo "pass $it: out of memory on: $files -> recompiling those at -O1"
  python3 -I - "$files" <<'PY'
import json, sys, subprocess, os
want = sys.argv[1].split()
for e in json.load(open(os.environ["ERF_BUILD"] + "/compile_commands.json")):
    if e["file"].split("/")[-1] in want:
        subprocess.run(e["command"].replace(" -O3 ", " -O1 "), shell=True, cwd=e["directory"], check=True)
PY
done
echo "gave up after 6 passes"; exit 1
