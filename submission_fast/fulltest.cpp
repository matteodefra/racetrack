/**
 * Random instances: the DP (ALGOWIN paper = Bekos et al., sequential) against the
 * same DP in parallel (dp_parallel, `threads` OpenMP threads). Both must return the
 * same optimal time on every instance.
 *
 * Usage: ./fulltest #trials n L mode threads [seed]
 *   #trials : number of random instances
 *   n       : number of cities (the last one is the first one again)
 *   L       : side of the box, coordinates are in [0, L)
 *   mode    : safe  -> smax = S/4 (S = time of the tour with |speed| <= 1);
 *                      the result is the true optimum.
 *             heur  -> smax = floor(sqrt(3L)) (heuristic choice);
 *                      the result is the optimum with |speed| <= smax only.
 *   threads : threads of DP-parallel (0 = all the cores). Compile with -fopenmp.
 *   seed    : seed of the random generator (default 1)
 *
 * Output: printed on screen, and one line appended to output/results.csv
 *         (header written if the file is new). Times: mean microseconds per instance,
 *         DP only (the graph and its filtering are built before the timers, and are
 *         the same for both DPs).
 *
 * Exit code: 0 if the two DPs agree on all the instances, 1 otherwise.
 */
#include <boost/heap/fibonacci_heap.hpp>
#include "utils/solvers.h"
#include "utils/utimer.hpp"
#include "utils/heuristic.h"   // only for the upper bound (pruning of the graph)
#include "utils/pruning.h"
#include <iostream>
#include <random>
#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// Bi-vectors allowed in every city. If `filter` is set, a bi-vector is kept only if
// it can belong to a solution of total time <= S (see filteringS in pruning.h).
auto bi_vectors_per_city(const vector<int>& xs, const vector<int>& ys, const vector<UniVec>& univecs,
                         bool filter, int S) {
  int n = xs.size();

  auto bivecs = bi_vectors(univecs);
  vector<vector<BiVec>> result(n);

  result[0].push_back({0,0});       // start at rest
  result[n-1].push_back({0,0});     // end at rest
  for (int i=1; i<n-1; i++){
    for (const auto& bivec: bivecs){
      // FIX: filteringS is true when the bi-vector CAN be part of a solution of time <= S,
      // i.e. when it must be KEPT. The old code kept it only if the test was false.
      // `cond` was also uninitialised when filter == false (undefined behaviour).
      bool keep = !filter || filteringS(bivec, univecs, xs, ys, i, S);
      if (keep) result[i].push_back(bivec);
    }
  }
  return result;
}


// Upper bound S on the optimal time: DP with the slowest vectors only (|speed| <= 1).
int trivialsolution(const vector<int> &xs, const vector<int> &ys) {
  auto univecs = uni_vectors(1);
  UniProfile X(xs, univecs);
  UniProfile Y(ys, univecs);
  vector<vector<BiVec>> bivecs = bi_vectors_per_city(xs, ys, univecs, false, 0);
  return solve(X, Y, bivecs, univecs, xs, ys);
}


// Better upper bound: optimal time with small speeds only (smax = 3, then 6).
// Every such tour is also a tour of the full model, so the result is a valid bound.
// Each step uses the bound of the step before to filter its own graph.
int upper_bound_small_speeds(const vector<int> &xs, const vector<int> &ys, int S, int smax) {
  int ub = S;
  for (int sm : {3, 6}) {
    if (sm >= smax) break;
    auto uv = uni_vectors(sm);
    UniProfile X(xs, uv), Y(ys, uv);
    auto bv = bi_vectors_per_city(xs, ys, uv, true, ub);
    auto h  = compute_heuristics(Heuristic::VectorTSP, bv, uv, xs, ys);
    ub = min(ub, solveastarfh(X, Y, bv, h, nullptr, ub));
  }
  return ub;
}





int main(int argc, const char * argv[]) {

  if (argc < 6 || argc > 7) {
    cout << "Usage ./executable #trials n L mode(safe|heur) threads [seed]" << endl;
    return 0;
  }

  int trials = std::atoi(argv[1]);   // number of experiments to run
  int n = std::atoi(argv[2]);        // number of cities to visit
  int L = std::atoi(argv[3]);        // dimension of the box
  string mode = argv[4];             // "safe" or "heur"
  int threads = std::atoi(argv[5]);  // threads of DP-parallel, 0 = all the cores
  unsigned seed = (argc == 7) ? std::atoi(argv[6]) : 1;
#ifdef _OPENMP
  if (threads <= 0) threads = omp_get_num_procs();
  omp_set_dynamic(0);                // exactly `threads` threads
#else
  cout << "compiled without -fopenmp: DP-parallel runs on one thread" << endl;
  threads = 1;
#endif
  if (mode != "safe" && mode != "heur") { cout << "mode must be safe or heur" << endl; return 1; }
  if (n < 4) { cout << "n must be at least 4" << endl; return 1; }
  if (L < 1) { cout << "L must be at least 1" << endl; return 1; }

  // heuristic smax: largest s with s*s <= 3L (integer square root, no rounding error)
  int smax_heur = 0;
  while ((long)(smax_heur + 1) * (smax_heur + 1) <= 3L * L) smax_heur++;

  std::mt19937 rng(seed);

  const int NALG = 2;
  const char* names[NALG] = {"DP", "DP-parallel"};
  long total_time[NALG] = {0, 0};
  long total_smax = 0;
  int mismatches = 0;

  for (int t=0; t<trials; t++) {

    vector<int> xs(n), ys(n);
    for (int i = 0; i<n; i++){
      xs[i] = rng() % L;
      ys[i] = rng() % L;
    }
    xs[n-1] = xs[0]; ys[n-1] = ys[0]; // setting last equal to first

    // S: time of the tour with |speed| <= 1 (always an upper bound on the optimal time)
    int S = trivialsolution( xs, ys );
    int smax = (mode == "safe") ? S/4 : smax_heur;
    total_smax += smax;

    // tighter upper bound UB <= S: used only to prune the graph (same graph for both DPs)
    int UB = upper_bound_small_speeds(xs, ys, S, smax);

    auto univecs = uni_vectors(smax);
    UniProfile X(xs, univecs);
    UniProfile Y(ys, univecs);
    vector<vector<BiVec>> bivecs = bi_vectors_per_city(xs, ys, univecs, true, UB);

    int cost[NALG];
    long time[NALG];
    { Utimer timer(names[0], &time[0]); cost[0] = solve( X, Y, bivecs, univecs, xs, ys ); }
    { Utimer timer(names[1], &time[1]); cost[1] = dp_parallel( X, Y, bivecs, threads ); }

    if (cost[0] != cost[1]) {
      mismatches++;
      cout << "MISMATCH on trial " << t << " (seed " << seed << ", n=" << n << ", L=" << L << "):"
           << "  DP=" << cost[0] << "  DP-parallel=" << cost[1] << endl;
    } else {
      for (int alg=0; alg<NALG; alg++) total_time[alg] += time[alg];
    }
  }

  int good = trials - mismatches;
  cout << "n=" << n << " L=" << L << " mode=" << mode << " threads=" << threads << " trials=" << trials
       << " seed=" << seed << "  mismatches=" << mismatches << endl;
  cout << "  mean smax " << (double)total_smax/trials << endl;
  if (good > 0) {
    for (int alg=0; alg<NALG; alg++)
      cout << "  " << names[alg] << ": mean time " << (double)total_time[alg]/good << " us" << endl;
    cout << "  speed-up " << (double)total_time[0]/total_time[1] << endl;
  }

  // one line per run in output/results.csv (all the runs of a sweep in one file)
  std::filesystem::create_directories("output");
  const string csv = "output/results.csv";
  bool is_new = !std::filesystem::exists(csv);
  ofstream out(csv, std::ios::app);
  if (is_new) out << "n,L,smax,threads,trials,seed,good,mismatches,dp_us,dp_parallel_us" << endl;
  out << n << "," << L << "," << (double)total_smax/trials << "," << threads << ","
      << trials << "," << seed << "," << good << "," << mismatches;
  for (int alg=0; alg<NALG; alg++) {
    out << ",";
    if (good) out << (double)total_time[alg]/good;
  }
  out << endl;

  return mismatches ? 1 : 0;
}
