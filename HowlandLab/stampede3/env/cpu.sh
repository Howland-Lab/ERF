#!/bin/bash
# CPU build/run environment on TACC Stampede3 (SPR nodes): Intel oneAPI 2024.0 icpx + Intel MPI 2021.11.
#   source HowlandLab/stampede3/env/cpu.sh
# NetCDF/HDF5 are PINNED: since the late-Sep 2026 TACC maintenance an unpinned `module load netcdf` pulls an
# Intel-built zlib that crashes csh/ls and the build (same pins as the lab's WRF build).
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
module reset > /dev/null 2>&1                    # intel/24.0 + impi/21.11 + cmake/4.1.1
module load phdf5/1.14.4 parallel-netcdf/4.9.2 pnetcdf/1.12.3 fftw3/3.3.10
export NETCDF_DIR=/opt/apps/intel24/impi21/parallel-netcdf/4.9.2/x86_64
export HDF5_ROOT=/opt/apps/intel24/impi21/phdf5/1.14.4
export HDF5_DIR=$HDF5_ROOT
export FFTW_DIR=$TACC_FFTW3_DIR FFTW_ROOT=$TACC_FFTW3_DIR
export CMAKE_PREFIX_PATH=$NETCDF_DIR:$HDF5_ROOT:$FFTW_DIR:$CMAKE_PREFIX_PATH
export ERF_BUILD=$HL_BUILD/spr
export ERF_EXE=$ERF_BUILD/Exec/erf_exec
