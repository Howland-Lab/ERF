#!/bin/bash
# Run several ERF cases one after another in one job (CPU).  Usage:
#     sbatch -J name -p spr -N 1 -n ranks -t hh:mm:ss submit_erf_chain.slurm /path/case1 /path/case2 ...
# Each case directory holds "inputs" and "input_sounding"; output goes to log.erf in that directory.
#SBATCH -o /scratch/11859/hsna/output_files/erfchain_.o%J
#SBATCH -e /scratch/11859/hsna/output_files/erfchain_.e%J
#SBATCH -p spr
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -A TG-ATM170028
#SBATCH -J erfchain
: "${ERF_HOME:?export ERF_HOME=/path/to/this/ERF/clone before sbatch}"
source "$ERF_HOME/HowlandLab/stampede3/env/cpu.sh"
export HOME=/home1/11859/hsna
export FI_PROVIDER=tcp
echo "node    : $(hostname)"; echo "ranks   : $SLURM_NTASKS"; echo "started : $(date)"
for d in "$@"; do
    cd "$d" || { echo "no dir $d"; continue; }
    echo "== $d  $(date)"
    ibrun "$ERF_EXE" inputs >& log.erf
    echo "   exit $?  $(grep -c -i nan mean_profiles.dat 2>/dev/null) nan lines  $(grep 'Total Time' log.erf | tail -1)"
done
echo "finished: $(date)"
