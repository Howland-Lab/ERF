#!/bin/bash
# ERF on NVIDIA H100 (h100 node: 4 GPUs, one MPI rank per GPU).  From a case directory:
#     export ERF_HOME=/path/to/ERF ; sbatch -J name -t hh:mm:ss submit_h100.slurm "$PWD" [overrides]
# GPUs want few, large boxes: amr.max_grid_size_x=160 amr.max_grid_size_y=160 gives 4 boxes on a 320 x 320 grid.
#SBATCH -J erf_h100
#SBATCH -o erf_h100.o%j
#SBATCH -p h100
#SBATCH -N 1
#SBATCH --ntasks-per-node 4
#SBATCH -t 00:30:00
#SBATCH -A TG-ATM170028
: "${ERF_HOME:?export ERF_HOME=/path/to/this/ERF/clone before sbatch}"
cd "$1" || exit 1; shift
source "$ERF_HOME/HowlandLab/stampede3/env/h100.sh"
echo "run dir: $(pwd)  node: $(hostname)  ranks: $SLURM_NTASKS  started: $(date)"; nvidia-smi -L
module -t list 2>&1 | tr '\n' ' '; echo; echo "mpirun: $(command -v mpirun || echo MISSING)"
if ! command -v mpirun > /dev/null || ldd "$ERF_EXE" | grep -q "not found"; then
    echo "ABORT: MPI or a shared library is missing:"; ldd "$ERF_EXE" | grep "not found"; exit 2; fi
ibrun "$ERF_EXE" inputs amrex.the_arena_is_managed=0 amr.max_grid_size_x=160 amr.max_grid_size_y=160 "$@" >& log.erf
echo "exit: $?  finished: $(date)"
