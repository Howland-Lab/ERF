#!/bin/bash
# Several ERF runs one after another on ONE pvc node (the pvc queue allows only 2 running jobs per user, so ensemble
# members are packed into jobs):   export ERF_HOME=... ; sbatch -J name -t hh:mm:ss submit_pvc_multi.slurm DIR1 DIR2 ...
#SBATCH -J erf_pvc_multi
#SBATCH -o erf_pvc_multi.o%j
#SBATCH -p pvc
#SBATCH -N 1
#SBATCH --ntasks-per-node 8
#SBATCH -t 05:00:00
#SBATCH -A TG-ATM170028
: "${ERF_HOME:?export ERF_HOME=/path/to/this/ERF/clone before sbatch}"
source "$ERF_HOME/HowlandLab/stampede3/env/pvc.sh"
echo "node: $(hostname)  started: $(date)  SYCL GPUs: $(sycl-ls 2>/dev/null | grep -c 'level_zero:gpu')"
for d in "$@"; do
    cd "$d" || { echo "missing $d"; continue; }
    echo "== $d start $(date)"
    ibrun "$ERF_EXE" inputs amr.max_grid_size_x=80 amr.max_grid_size_y=160 >& log.erf
    echo "== $d exit $? end $(date)"; grep "Coarse STEP" log.erf | tail -1
done
