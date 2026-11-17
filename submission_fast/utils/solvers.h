#pragma once
/**
 * Solvers on the layered graph of bi-vectors.
 *
 *   solvedijkstra : Dijkstra  (A* with h = 0), Fibonacci heap
 *   solveastarfh  : A*, Fibonacci heap, heuristic given as a table
 *                   heuristics[city][bivec_idx]  (see heuristic.h)
 *   solve         : the DP (ALGOWIN paper = Bekos et al.), sequential
 *   dp_parallel   : the same DP; in every layer the sources are split between
 *                   `threads` OpenMP threads, then a parallel min-reduction
 *
 * Layer i contains the bi-vectors allowed when crossing city i.
 * The weight of the edge (i, a) -> (i+1, b) is the smallest time tau that is
 * feasible on both axes: min_value(X.get(...), Y.get(...)).
 */
#include <climits>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include "profiles.cpp"
#include "fibonacciheap.h"

static const int INF_COST = INT_MAX / 2;   // "not reached yet" (no overflow when adding a weight)

struct SearchStats {
  long pushed    = 0;   // nodes inserted in the heap
  long popped    = 0;   // nodes extracted (= expanded)
  long decreased = 0;   // key decreases of nodes already in the heap
  long reopened  = 0;   // closed nodes inserted again (only if h is inconsistent)
};


/**
 * A* with a Fibonacci heap. If `heuristics` is nullptr, h = 0 (Dijkstra).
 *
 * Returns the optimal time, or INF_COST if the last city cannot be reached.
 * The result is optimal if h is admissible (h <= true time left).
 * UB (optional): a known upper bound on the optimal time. Nodes with g + h > UB are
 * not inserted (they cannot be on an optimal path if h is admissible).
 * If h is also consistent, every node is expanded at most once; if it is not,
 * a closed node is re-opened whenever a strictly shorter path to it is found,
 * so the result stays optimal as long as h is admissible.
 */
inline int astar_core(const UniProfile& X, const UniProfile& Y,
                      const std::vector<std::vector<BiVec>>& bivecs,
                      const std::vector<std::vector<int>>* heuristics,
                      SearchStats* stats, int UB = INF_COST) {
  const int n = X.n;
  SearchStats local;
  SearchStats& st = stats ? *stats : local;
  st = SearchStats();

  auto h_of = [&](int city, int idx) { return heuristics ? (*heuristics)[city][idx] : 0; };

  // handles[c][j] : handle of node (c,j) while it is in the heap, NOHANDLE otherwise
  // best_g[c][j]  : shortest time found so far to reach node (c,j)
  std::vector<std::vector<Handle>> handles(n);
  std::vector<std::vector<int>> best_g(n);
  for (int i = 0; i < n; i++) {
    handles[i].assign(bivecs[i].size(), NOHANDLE);
    best_g[i].assign(bivecs[i].size(), INF_COST);
  }

  FibHeap heap;
  best_g[0][0] = 0;
  handles[0][0] = heap.push(heapdata(Conff{0, 0, h_of(0, 0), 0, h_of(0, 0)}));
  st.pushed++;

  while (!heap.empty()) {
    const Conff from = heap.top().config;   // copy: the node is removed just below
    heap.pop();
    st.popped++;
    handles[from.city][from.bivec_idx] = NOHANDLE;

    if (from.city == n - 1) return from.distance;   // h(last) = 0, so f = g: optimal

    const BiVec& bv_from = bivecs[from.city][from.bivec_idx];
    const int next = from.city + 1;

    const int* h_next = heuristics ? (*heuristics)[next].data() : nullptr;
    int* g_next = best_g[next].data();
    const int nb = bivecs[next].size();
    const BiVec* to = bivecs[next].data();
    int sx, sy;                                    // rows of intervals of this source
    const Interval* rx = X.row(from.city, bv_from.ind_x, sx);   // (same access as the DPs)
    const Interval* ry = Y.row(from.city, bv_from.ind_y, sy);
    for (int j = 0; j < nb; j++) {
      // cheap tests first (tau >= 0), before the intervals are read
      if (from.distance >= g_next[j]) continue;                       // cannot improve
      if (h_next && from.distance + h_next[j] > UB) continue;         // cannot be optimal
      int tau = min_value(rx[to[j].ind_x * sx], ry[to[j].ind_y * sy]);

      int g = from.distance + tau;                 // g(next, j) through this path
      if (g >= best_g[next][j]) continue;          // not an improvement
      int h = h_of(next, j);
      if (g + h > UB) continue;                    // cannot be optimal

      bool was_closed = (best_g[next][j] != INF_COST) && (handles[next][j] == NOHANDLE);
      best_g[next][j] = g;
      heapdata data(Conff{next, j, g + h, g, h});

      if (handles[next][j] == NOHANDLE) {          // new node, or closed node to re-open
        handles[next][j] = heap.push(data);
        st.pushed++;
        if (was_closed) st.reopened++;
      } else {                                      // already in the heap with a larger key
        heap.update(handles[next][j], data);
        st.decreased++;
      }
    }
  }
  return INF_COST;
}


inline int solvedijkstra(const UniProfile& X, const UniProfile& Y,
                  const std::vector<std::vector<BiVec>>& bivecs,
                  SearchStats* stats = nullptr) {
  return astar_core(X, Y, bivecs, nullptr, stats);
}


inline int solveastarfh(const UniProfile& X, const UniProfile& Y,
                 const std::vector<std::vector<BiVec>>& bivecs,
                 const std::vector<std::vector<int>>& heuristics,
                 SearchStats* stats = nullptr, int UB = INF_COST) {
  return astar_core(X, Y, bivecs, &heuristics, stats, UB);
}


/**
 * Kernel of both DPs: relax all the targets b of the next layer from ONE source
 * (cost ca, intervals rows rx / ry of its x / y uni-vectors):
 *     out[b] = min(out[b], ca + tau(source, b)).
 * Plain pointers in the arguments: no aliasing, the compiler keeps everything in registers.
 */
inline void relax_from_source(int ca, const BiVec* __restrict to, int nb,
                              const Interval* __restrict rx, int sx,
                              const Interval* __restrict ry, int sy, int* __restrict out) {
  for (int b = 0; b < nb; b++) {
    int total = ca + min_value(rx[to[b].ind_x * sx], ry[to[b].ind_y * sy]);
    if (total < out[b]) out[b] = total;
  }
}


/**
 * DP of the ALGOWIN paper (the same DP as Bekos et al.): costs of layer i+1 from
 * costs of layer i, all pairs (a, b). Sequential.
 * (univecs, xs, ys are not used; they are kept so the call sites do not change.)
 */
inline int solve(const UniProfile& X, const UniProfile& Y,
          const std::vector<std::vector<BiVec>>& bivecs, const std::vector<UniVec>& univecs,
          const std::vector<int>& xs, const std::vector<int>& ys) {
  const int n = X.n;
  std::vector<int> cost(bivecs[0].size(), 0);    // layer 0: start, cost 0

  for (int i = 0; i < n - 1; i++) {
    const int na = bivecs[i].size(), nb = bivecs[i + 1].size();
    const BiVec* to = bivecs[i + 1].data();
    std::vector<int> next_cost(nb, INF_COST);
    for (int a = 0; a < na; a++) {
      const int ca = cost[a];
      if (ca >= INF_COST) continue;                // source not reachable
      const BiVec bv_from = bivecs[i][a];
      int sx, sy;
      const Interval* rx = X.row(i, bv_from.ind_x, sx);
      const Interval* ry = Y.row(i, bv_from.ind_y, sy);
      relax_from_source(ca, to, nb, rx, sx, ry, sy, next_cost.data());
    }
    cost.swap(next_cost);
  }
  return cost[0];   // last layer holds one bi-vector: stopped in the last city
}


/**
 * Parallel version of the same DP (same result as solve).
 * For every layer i -> i+1:
 *   1. the sources a of layer i are split between the threads (dynamic schedule,
 *      the unreachable sources are skipped); every thread computes, in its own
 *      array best[t][b], the min over ITS sources of cost[a] + tau(a, b);
 *   2. min-reduction: the targets b are split between the threads and
 *      next_cost[b] = min over t of best[t][b].
 * Each thread runs the same kernel as the sequential DP (relax_from_source), so the
 * only differences are the split of the work and the reduction. Compile with -fopenmp, else it is sequential.
 */
inline int dp_parallel(const UniProfile& X, const UniProfile& Y,
                       const std::vector<std::vector<BiVec>>& bivecs, int threads) {
#ifdef _OPENMP
  if (threads <= 0) threads = omp_get_num_procs();
#else
  threads = 1;
#endif
  const int n = X.n;
  std::vector<int> cost(bivecs[0].size(), 0);
  std::vector<std::vector<int>> best(threads);   // one array per thread, reused

  for (int i = 0; i < n - 1; i++) {
    const int na = bivecs[i].size(), nb = bivecs[i + 1].size();
    const BiVec* from = bivecs[i].data();
    const BiVec* to = bivecs[i + 1].data();
    std::vector<int> next_cost(nb);

    #pragma omp parallel num_threads(threads)
    {
#ifdef _OPENMP
      const int t = omp_get_thread_num(), team = omp_get_num_threads();
#else
      const int t = 0, team = 1;
#endif
      best[t].assign(nb, INF_COST);
      int* const mine = best[t].data();

      // 1. every thread: its share of the sources (same kernel as the sequential DP)
      #pragma omp for schedule(dynamic, 8)
      for (int a = 0; a < na; a++) {
        const int ca = cost[a];
        if (ca >= INF_COST) continue;
        const BiVec bv_from = from[a];
        int sx, sy;
        const Interval* rx = X.row(i, bv_from.ind_x, sx);
        const Interval* ry = Y.row(i, bv_from.ind_y, sy);
        relax_from_source(ca, to, nb, rx, sx, ry, sy, mine);
      }
      // (implicit barrier: all the arrays are complete)

      // 2. min-reduction, the targets are split between the threads
      #pragma omp for schedule(static)
      for (int b = 0; b < nb; b++) {
        int m = INF_COST;
        for (int u = 0; u < team; u++) m = std::min(m, best[u][b]);   // team <= threads
        next_cost[b] = m;
      }
    }
    cost.swap(next_cost);
  }
  return cost[0];
}
