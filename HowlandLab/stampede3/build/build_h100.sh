#!/bin/bash
# Configure + build ERF for NVIDIA H100 (CUDA, sm_90) into _build/h100.  Runs on the login node (nvcc needs no GPU).
# Usage: build_h100.sh [-j N]    (-j 4: nvcc + AMReX headers use several GB per process)
J=4; [ "$1" = "-j" ] && J=$2
source "$(dirname "$0")/../env/h100.sh"
cmake -S "$ERF_HOME" -B "$ERF_BUILD" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CUDA_HOST_COMPILER=g++ \
      -DERF_DIM=3 -DERF_ENABLE_MPI=ON -DERF_ENABLE_CUDA=ON -DAMReX_CUDA_ARCH=9.0 \
      -DERF_ENABLE_FFT=ON -DERF_ENABLE_NETCDF=OFF \
      -DERF_ENABLE_TESTS=OFF -DERF_ENABLE_UNIT_TESTS=OFF -DERF_ENABLE_FCOMPARE=OFF \
      -DERF_ENABLE_DOCUMENTATION=OFF > "$HL_BUILD/configure_h100.log" 2>&1
echo "configure exit $?"
cmake --build "$ERF_BUILD" -j $J -- -k > "$HL_BUILD/build_h100.log" 2>&1
echo "build exit $?"; grep -E " error" "$HL_BUILD/build_h100.log" | sort -u | head; ls -la "$ERF_EXE"
