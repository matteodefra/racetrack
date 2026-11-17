#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# DP (= ALGOWIN = Bekos et al.) on 1 thread against DP-parallel on JOB_THREADS
# threads, on the same instances. smax = floor(sqrt(3L)).
#
# Phase 1 (many jobs at the same time, each job on its own JOB_THREADS physical cores):
#   sweep "n": L fixed (one sweep per value of SWEEP_N_L), n grows
#   sweep "L": n fixed, L grows
# Phase 2 (after phase 1, one job at a time, machine otherwise idle):
#   sweep "T": scaling of DP-parallel in the number of threads, at (SCALE_N, SCALE_L)
#              (speed-up = t(DP) / t(DP-parallel)).
#
# Before: build the binary with the Makefile (make output/fulltest).
#   WARNING: "make clean" deletes output/*, so also the results.
#
# Usage, in the folder of the Makefile:
#   chmod +x run_sweeps_parallel.sh
#   nohup ./run_sweeps_parallel.sh > sweeps.log 2>&1 &
#   ./run_sweeps_parallel.sh merge        # partial results while it runs
#
# Result: output/sweeps_merged.csv, one row per job:
#   sweep,status,cores,n,L,smax,threads,trials,seed,good,mismatches,dp_us,dp_parallel_us
#   (status = ok | mismatch | timeout | error | running; times = mean us per instance)
#
# Cores (tuned for 2 x AMD EPYC 9654, 755 GB RAM):
#   physical cores 0-191 (192-383 are SMT siblings: NOT used); NUMA node0 = 0-95, node1 = 96-191.
#   Phase 1: WORKERS = 12 jobs x JOB_THREADS = 16 cores: 6 jobs in each NUMA node, each job
#   on 2 full CCDs (8 cores = one L3 cache), so no job shares a cache or crosses a node.
#   The sequential DP of a job runs on one of its cores, then DP-parallel uses all 16.
# ---------------------------------------------------------------------------
set -u
cd "$(dirname "$0")"

TRIALS=${TRIALS:-10}
SEED=${SEED:-1}
TIMEOUT=${TIMEOUT:-72h}               # max time of one job (all its trials)
BIN=${BIN:-$(pwd)/output/fulltest}

JOB_THREADS=${JOB_THREADS:-16}        # threads of DP-parallel in phase 1 (= cores per job)
WORKERS=${WORKERS:-12}                # phase 1 jobs at the same time
FIRST_CORE=${FIRST_CORE:-0}           # worker k uses cores FIRST_CORE + k*JOB_THREADS ...

SWEEP_N_L=${SWEEP_N_L:-"100 250 500"}    # one n sweep for EACH value of L
SWEEP_N_VALUES=${SWEEP_N_VALUES:-"10 15 20 25 30 35 40 50 60 80 100"}
SWEEP_L_N=${SWEEP_L_N:-50}
SWEEP_L_VALUES=${SWEEP_L_VALUES:-"100 150 200 250 300 350 400 450 500"}

SCALE_N=${SCALE_N:-50}                # phase 2 point
SCALE_L=${SCALE_L:-500}
SCALE_TRIALS=${SCALE_TRIALS:-5}
SCALE_THREADS=${SCALE_THREADS:-"1 2 4 8 16 32 48 64 96 128 192"}   # cores 0..T-1

JOBS=output/jobs
MERGED=output/sweeps_merged.csv
HEADER="n,L,smax,threads,trials,seed,good,mismatches,dp_us,dp_parallel_us"
EMPTY=",,,,"                          # good, mismatches, dp_us, dp_parallel_us

# ----------------------------------------------------------------- merge
merge () {
  echo "sweep,status,cores,$HEADER" > "$MERGED.tmp"
  for d in "$JOBS"/*/; do
    [ -f "$d/job" ] || continue
    read -r sweep n L thr trials < "$d/job"
    local cores status code
    cores=$(cat "$d/cores" 2>/dev/null || echo "")
    if [ ! -f "$d/exit_code" ]; then
      echo "$sweep,running,$cores,$n,$L,,$thr,$trials,$SEED$EMPTY" >> "$MERGED.tmp"; continue
    fi
    code=$(cat "$d/exit_code")
    if [ -f "$d/output/results.csv" ] && [ "$(wc -l < "$d/output/results.csv")" -ge 2 ]; then
      case $code in 0) status=ok ;; 1) status=mismatch ;; *) status=error ;; esac
      echo "$sweep,$status,$cores,$(tail -n 1 "$d/output/results.csv")" >> "$MERGED.tmp"
    else
      if [ "$code" = 124 ]; then status=timeout; else status=error; fi
      echo "$sweep,$status,$cores,$n,$L,,$thr,$trials,$SEED$EMPTY" >> "$MERGED.tmp"
    fi
  done
  mv "$MERGED.tmp" "$MERGED"
  echo "merged $(($(wc -l < "$MERGED") - 1)) rows in $MERGED"
}

if [ "${1:-}" = "merge" ]; then merge; exit 0; fi

# ----------------------------------------------------------------- setup
[ -x "$BIN" ] || { echo "binary $BIN not found: run 'make output/fulltest' first"; exit 1; }
mkdir -p output
[ -d "$JOBS" ] && mv "$JOBS" "output/jobs_old_$(date +%Y%m%d_%H%M%S)"   # keep old runs
mkdir -p "$JOBS"

PIN=1
command -v taskset >/dev/null || { echo "taskset not found: jobs are not pinned"; PIN=0; }
LAST=$((FIRST_CORE + WORKERS * JOB_THREADS - 1))
[ "$LAST" -lt 192 ] || echo "WARNING: phase 1 uses core $LAST (SMT sibling or missing core)"
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false

# run_job <sweep> <n> <L> <threads> <trials> <cores>
run_job () {
  local sweep=$1 n=$2 L=$3 thr=$4 trials=$5 cores=$6
  local d="$JOBS/${sweep}_n${n}_L${L}_thr${thr}"
  mkdir -p "$d"
  echo "$sweep $n $L $thr $trials" > "$d/job"
  echo "$cores" > "$d/cores"
  echo "[$(date +%F' '%T)] start $d (cores $cores)"
  (
    cd "$d"
    if [ $PIN = 1 ]; then
      timeout "$TIMEOUT" taskset -c "$cores" "$BIN" "$trials" "$n" "$L" heur "$thr" "$SEED" > log.txt 2>&1
    else
      timeout "$TIMEOUT" "$BIN" "$trials" "$n" "$L" heur "$thr" "$SEED" > log.txt 2>&1
    fi
    echo $? > exit_code
  )
  echo "[$(date +%F' '%T)] end   $d (exit $(cat "$d/exit_code"))"
}

# ----------------------------------------------------------------- phase 1
# queue, largest jobs first (cost ~ n * L^3); one line: sweep n L
{
  for LL in $SWEEP_N_L; do
    for n in $SWEEP_N_VALUES; do echo "$((n * LL * LL * LL / 1000)) n $n $LL"; done
  done
  for L in $SWEEP_L_VALUES; do echo "$((SWEEP_L_N * L * L * L / 1000)) L $SWEEP_L_N $L"; done
} | sort -k1,1nr | awk '{print $2, $3, $4}' > "$JOBS/queue"
echo 0 > "$JOBS/next"

next_job () {   # next line of the queue (empty if none), with a lock
  (
    flock 9
    local i total
    i=$(cat "$JOBS/next"); total=$(wc -l < "$JOBS/queue")
    if [ "$i" -lt "$total" ]; then
      echo $((i + 1)) > "$JOBS/next"
      sed -n "$((i + 1))p" "$JOBS/queue"
    fi
  ) 9> "$JOBS/lock"
}

worker () {     # worker <cores>
  local cores=$1
  while true; do
    local line; line=$(next_job)
    [ -z "$line" ] && break
    read -r sweep n L <<< "$line"
    run_job "$sweep" "$n" "$L" "$JOB_THREADS" "$TRIALS" "$cores"
  done
}

echo "[$(date +%F' '%T)] phase 1: $(wc -l < "$JOBS/queue") jobs, $WORKERS at a time, $JOB_THREADS cores each"
for ((w = 0; w < WORKERS; w++)); do
  c=$((FIRST_CORE + w * JOB_THREADS))
  worker "$c-$((c + JOB_THREADS - 1))" &
done
wait
merge

# ----------------------------------------------------------------- phase 2
echo "[$(date +%F' '%T)] phase 2: scaling in threads at n=$SCALE_N, L=$SCALE_L"
for T in $SCALE_THREADS; do
  run_job T "$SCALE_N" "$SCALE_L" "$T" "$SCALE_TRIALS" "0-$((T - 1))"
  merge
done
echo "[$(date +%F' '%T)] all jobs done"
