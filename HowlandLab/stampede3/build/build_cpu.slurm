#!/bin/bash
# CPU build on a compute node (no 8 GB vmem limit).  Configure first with build_cpu.sh on the login node
# (it may stop at the first OOM), then:   export ERF_HOME=/path/to/ERF ; sbatch build_cpu.slurm
#SBATCH -J erf_build_cpu
#SBATCH -o erf_build_cpu.o%j
#SBATCH -p spr
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:30:00
#SBATCH -A TG-ATM170028
: "${ERF_HOME:?export ERF_HOME=/path/to/this/ERF/clone before sbatch}"
source "$ERF_HOME/HowlandLab/stampede3/env/cpu.sh"
cmake --build "$ERF_BUILD" -j 48 2>&1 | grep -E "error|Linking|Built target erf_exec|out of memory"
ls -la "$ERF_EXE"
