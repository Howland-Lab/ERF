#!/bin/bash
# Tiny serial ERF run on the login node (32^3 GABLS1, 20 s of model time) to prove the executable works.
# Usage: ./smoke_test_login.sh      Output in ERF_runs/smoke/ (log.erf, plt*, *.dat).
source "$(dirname "$0")/../env/cpu.sh"; R=$HL_DIR
[ -x "$ERF_EXE" ] || { echo "no executable: $ERF_EXE"; exit 1; }
S=$HL_DIR/runs/smoke; rm -rf "$S"; mkdir -p "$S"; cd "$S" || exit 1
cp "$R/cases/sbl_dx6p25/inputs" "$R/cases/sbl_dx6p25/input_sounding" .
ulimit -s unlimited
export I_MPI_FABRICS=shm            # login node has no OFI fabric: singleton MPI needs shm
/usr/bin/time -f "wall %e s  maxRSS %M kB" timeout 600 "$ERF_EXE" inputs \
   amr.n_cell="32 32 32" amr.max_grid_size_x=16 amr.max_grid_size_y=16 \
   stop_time=20.0 erf.fixed_dt=0.1 erf.plot_per_1=10.0 erf.profile_int=50 erf.check_per=-1 \
   erf.sum_interval=50 > log.erf 2>&1
echo "exit $?"; grep -E "Coarse TimeStep|Total Time|ERROR|Error|error|abort|Adding random" log.erf | tail -8
ls -d plt* 2>/dev/null; ls *.dat 2>/dev/null; echo "mean_profiles lines: $(wc -l < mean_profiles.dat 2>/dev/null)"
