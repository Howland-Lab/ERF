#!/bin/bash
# ERF on CPU (SPR).  From inside a case directory holding "inputs" and "input_sounding":
#     export ERF_HOME=/path/to/ERF ; sbatch -J name -t hh:mm:ss -N n --ntasks-per-node p submit_cpu.slurm "$PWD" [overrides]
# Trailing arguments override entries in inputs (e.g. stop_time=600 for a test).
#SBATCH -J erf_cpu
#SBATCH -o erf_cpu.o%j
#SBATCH -p spr
#SBATCH -N 1
#SBATCH --ntasks-per-node 32
#SBATCH -t 00:30:00
#SBATCH -A TG-ATM170028
: "${ERF_HOME:?export ERF_HOME=/path/to/this/ERF/clone before sbatch}"
cd "$1" || exit 1; shift
source "$ERF_HOME/HowlandLab/stampede3/env/cpu.sh"
echo "run dir: $(pwd)  node: $(hostname)  ranks: $SLURM_NTASKS  started: $(date)"
export FI_PROVIDER=tcp
ibrun "$ERF_EXE" inputs "$@" >& log.erf
echo "exit: $?  finished: $(date)"
