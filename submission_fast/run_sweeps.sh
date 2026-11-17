#!/usr/bin/env bash
# Sequential version: the same runs as run_sweeps_parallel.sh, but one job at a time
# (WORKERS=1). All the parameters of run_sweeps_parallel.sh can be given here too, e.g.
#   JOB_THREADS=32 nohup ./run_sweeps.sh > sweeps.log 2>&1 &
cd "$(dirname "$0")"
WORKERS=${WORKERS:-1} exec ./run_sweeps_parallel.sh "$@"
