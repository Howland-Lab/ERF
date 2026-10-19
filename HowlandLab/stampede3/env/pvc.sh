#!/bin/bash
# Intel Data Center GPU Max 1550 ("PVC", SYCL) environment on Stampede3 (pvc nodes: 4 GPUs = 8 tiles, 96 cores).
#   source HowlandLab/stampede3/env/pvc.sh
# intel/25.1 (oneAPI 2025.1 icpx with SYCL) + impi/21.15.  No NetCDF, no FFT.
# 14 Oct 2026: (1) the default intel/24.0 is too old for this AMReX: AMReX_Parser_Y.H uses sycl::ext::intel::math::jn/yn,
#   which oneAPI 2024.0 does not have (build dies after 1 min).  (2) TACC's intel/25.1 module does not put libumf.so.0
#   (Unified Memory Framework) on LD_LIBRARY_PATH; without it the SYCL runtime cannot load its Level Zero / OpenCL
#   adapters and sees NO device ("No device of requested type 'info::device_type::gpu'", even sycl-ls shows 0 platforms).
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
module reset > /dev/null 2>&1
module load intel/25.1
export LD_LIBRARY_PATH=/scratch/projects/compilers/intel25.1/umf/0.10/lib:$LD_LIBRARY_PATH
# Each Max 1550 has two tiles; FLAT exposes every tile as its own device (8 per node), one MPI rank per tile.
export ZE_FLAT_DEVICE_HIERARCHY=FLAT
export ERF_BUILD=$HL_BUILD/pvc
export ERF_EXE=$ERF_BUILD/Exec/erf_exec
