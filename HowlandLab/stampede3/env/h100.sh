#!/bin/bash
# NVIDIA H100 (CUDA) environment on Stampede3 (h100 nodes: 4 x H100, 96 cores).
#   source HowlandLab/stampede3/env/h100.sh
# gcc 13.2 + CUDA 12.8 + Intel MPI 21.9 (gcc wrappers).  No NetCDF (plotfiles only).
# 13 Oct 2026: after `module reset; module load gcc` Intel MPI was NOT active on the h100 nodes ("mpirun: command not
# found", exit 127 after 3 s).  impi/21.9 is now loaded explicitly, with a fallback to the exact MPI the binary was
# linked against, and module messages are no longer hidden.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
module reset > /dev/null 2>&1
module load gcc/13.2.0
module load cuda/12.8
module load impi/21.9 2>/dev/null || true
IMPI=/scratch/projects/compilers/intel23.1/mpi/2021.9.0
if ! command -v mpirun > /dev/null 2>&1 || [[ "$(command -v mpirun)" != "$IMPI"/* ]]; then
    export I_MPI_ROOT=$IMPI PATH=$IMPI/bin:$PATH
    export LD_LIBRARY_PATH=$IMPI/lib/release:$IMPI/lib:$IMPI/libfabric/lib:$LD_LIBRARY_PATH
    export FI_PROVIDER_PATH=$IMPI/libfabric/lib/prov
fi
# 19 Oct 2026: on the h100 nodes Intel MPI could not open its OFI fabric ("OFI fi_getinfo() failed ... No data
# available", MPI_Init abort, exit 143).  Single-node jobs only need shared memory.
export I_MPI_FABRICS=shm
export ERF_BUILD=$HL_BUILD/h100
export ERF_EXE=$ERF_BUILD/Exec/erf_exec
