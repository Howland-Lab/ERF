#!/bin/bash
# ERF on Intel PVC (pvc node: 4 x Max 1550 = 8 tiles, one MPI rank per tile).  From a case directory:
#     export ERF_HOME=/path/to/ERF ; sbatch -J name -t hh:mm:ss submit_pvc.slurm "$PWD" [overrides]
# 8 boxes for 8 tiles on a 320 x 320 grid: amr.max_grid_size_x=80 amr.max_grid_size_y=160.
#SBATCH -J erf_pvc
#SBATCH -o erf_pvc.o%j
#SBATCH -p pvc
#SBATCH -N 1
#SBATCH --ntasks-per-node 8
#SBATCH -t 04:00:00
#SBATCH -A TG-ATM170028
: "${ERF_HOME:?export ERF_HOME=/path/to/this/ERF/clone before sbatch}"
cd "$1" || exit 1; shift
source "$ERF_HOME/HowlandLab/stampede3/env/pvc.sh"
echo "run dir: $(pwd)  node: $(hostname)  ranks: $SLURM_NTASKS  started: $(date)"
module -t list 2>&1 | tr '\n' ' '; echo
[ -x "$ERF_EXE" ] || { echo "ABORT: $ERF_EXE not built"; exit 2; }
NGPU=$(sycl-ls 2>/dev/null | grep -c "level_zero:gpu"); echo "SYCL Level Zero GPU devices: $NGPU"
[ "$NGPU" -ge 1 ] || { echo "ABORT: no GPU visible to SYCL (see env/pvc.sh, libumf)"; SYCL_UR_TRACE=1 sycl-ls 2>&1 | head -8; exit 3; }
ibrun "$ERF_EXE" inputs amr.max_grid_size_x=80 amr.max_grid_size_y=160 "$@" >& log.erf
echo "exit: $?  finished: $(date)"
