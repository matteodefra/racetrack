/**
 * Self checks. Usage: ./tests [trials=200] [seed=1]
 *
 *  1. Fibonacci heap     : random push / pop / update against a std::multiset.
 *  2. compute_interval   : against a brute-force search on one axis.
 *  3. Heuristics         : admissible (h <= true time left) on EVERY node of random
 *                          instances, against the exact time left computed by backward DP;
 *                          consistent (h(u) <= w(u,v) + h(v)) is required for the poor one only.
 *  4. Solvers            : Dijkstra, A* (both heuristics) and DP against a reference
 *                          Dijkstra written with std::priority_queue.
 *
 * Exit code 0 = everything passed.
 */
#include <boost/heap/fibonacci_heap.hpp>
#include "utils/solvers.h"
#include "utils/heuristic.h"
#include "utils/pruning.h"
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <tuple>

using namespace std;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { failures++; if (failures <= 20) cout << "  FAIL: " << msg << endl; } } while (0)


/* ---------------- 1. heap ---------------- */
void test_heap(int trials, mt19937& rng) {
  cout << "[1] Fibonacci heap vs std::multiset ... " << flush;
  int before = failures;
  for (int t = 0; t < trials; t++) {
    FibHeap heap;
    vector<Handle> handles;                   // handle of each id
    vector<int> key;                          // current key of each id (-1 = not in heap)
    set<pair<int,int>> ref;                   // (key, id) pairs in the heap
    for (int op = 0; op < 400; op++) {
      int r = rng() % 10;
      if (r < 4 || ref.empty()) {             // push
        int id = key.size(), k = rng() % 50;
        key.push_back(k);
        handles.push_back(heap.push(heapdata(Conff{id, 0, k, k, 0})));
        ref.insert({k, id});
      } else if (r < 7) {                     // decrease a random key
        auto it = ref.begin(); advance(it, rng() % ref.size());
        int id = it->second, k = it->first;
        if (k == 0) continue;
        int nk = rng() % k;                    // strictly smaller
        ref.erase(it); ref.insert({nk, id});
        key[id] = nk;
        heap.update(handles[id], heapdata(Conff{id, 0, nk, nk, 0}));
      } else {                                // pop
        int k = heap.top().config.cost;
        CHECK(k == ref.begin()->first, "top key " << k << " != reference " << ref.begin()->first);
        int id = heap.top().config.city;
        CHECK(key[id] == k, "stored key of node differs");
        heap.pop();
        ref.erase({k, id});
        key[id] = -1;
      }
      CHECK(heap.size() == ref.size(), "size mismatch");
    }
    while (!heap.empty()) {                   // drain: must come out sorted
      int k = heap.top().config.cost, id = heap.top().config.city;
      CHECK(k == ref.begin()->first, "drain order");
      heap.pop(); ref.erase({k, id});
    }
  }
  cout << (failures == before ? "ok" : "FAILED") << endl;
}


/* ---------------- 2. intervals ---------------- */
void test_intervals() {
  cout << "[2] compute_interval vs brute force ... " << flush;
  int before = failures;
  const int T = 40;
  for (int s = -5; s <= 5; s++) {
    vector<set<pair<int,int>>> R(T + 1);       // reachable (position, speed) after t moves, from (0,s)
    R[0].insert({0, s});
    for (int t = 0; t < T; t++)
      for (auto [p, v] : R[t]) for (int d = -1; d <= 1; d++) R[t + 1].insert({p + v + d, v + d});
    for (int sp = -5; sp <= 5; sp++) for (int xp = -30; xp <= 30; xp++) {
      bool ever = false;
      for (int t = 0; t <= T; t++) if (R[t].count({xp, sp})) ever = true;
      if (!ever) continue;                      // too far for the brute force horizon
      Interval iv = compute_interval(0, s, xp, sp);
      for (int t = 0; t <= T - 8; t++) {
        bool pred = (iv.t1 != -1 && t >= iv.t1 && t <= iv.t2) || t >= iv.t3;
        CHECK(pred == (bool)R[t].count({xp, sp}),
              "interval (0," << s << ")->(" << xp << "," << sp << ") wrong at t=" << t);
      }
    }
  }
  cout << (failures == before ? "ok" : "FAILED") << endl;
}


/* ---------------- instance generation ---------------- */
struct Instance {
  vector<int> xs, ys;
  vector<UniVec> univecs;
  vector<vector<BiVec>> bivecs;
};

Instance random_instance(mt19937& rng, int n, int L, int smax, bool filter) {
  Instance in;
  in.xs.resize(n); in.ys.resize(n);
  for (int i = 0; i < n; i++) { in.xs[i] = rng() % L; in.ys[i] = rng() % L; }
  in.xs[n-1] = in.xs[0]; in.ys[n-1] = in.ys[0];
  in.univecs = uni_vectors(smax);
  auto all = bi_vectors(in.univecs);
  in.bivecs.assign(n, {});
  in.bivecs[0].push_back({0,0});
  in.bivecs[n-1].push_back({0,0});

  int S = 0;
  if (filter) {   // S = DP with |speed| <= 1
    auto u1 = uni_vectors(1);
    UniProfile X1(in.xs, u1), Y1(in.ys, u1);
    vector<vector<BiVec>> b1(n);
    auto all1 = bi_vectors(u1);
    b1[0].push_back({0,0}); b1[n-1].push_back({0,0});
    for (int i = 1; i < n - 1; i++) b1[i] = all1;
    S = solve(X1, Y1, b1, u1, in.xs, in.ys);
  }
  for (int i = 1; i < n - 1; i++)
    for (const auto& b : all)
      if (!filter || filteringS(b, in.univecs, in.xs, in.ys, i, S)) in.bivecs[i].push_back(b);
  return in;
}

int edge_weight(const UniProfile& X, const UniProfile& Y, const Instance& in, int i, int a, int b) {
  return min_value(X.get(i, in.bivecs[i][a].ind_x, in.bivecs[i+1][b].ind_x),
                   Y.get(i, in.bivecs[i][a].ind_y, in.bivecs[i+1][b].ind_y));
}


/* ---------------- 3. heuristics ---------------- */
struct HStats { long nodes = 0, inadmissible = 0, edges = 0, inconsistent = 0; };

void check_heuristic(const char* name, const Instance& in, const UniProfile& X, const UniProfile& Y,
                     const vector<vector<int>>& h, const vector<vector<int>>& ctg, HStats& st,
                     bool report_admissible, bool report_consistent) {
  int n = in.xs.size();
  for (int i = 0; i < n - 1; i++)
    for (size_t a = 0; a < in.bivecs[i].size(); a++) {
      st.nodes++;
      if (h[i][a] > ctg[i][a]) {
        st.inadmissible++;
        if (report_admissible) CHECK(false, name << " not admissible: city " << i << " h=" << h[i][a] << " true=" << ctg[i][a]);
      }
      for (size_t b = 0; b < in.bivecs[i+1].size(); b++) {
        st.edges++;
        if (h[i][a] > edge_weight(X, Y, in, i, a, b) + h[i+1][b]) {
          st.inconsistent++;
          if (report_consistent) CHECK(false, name << " not consistent: city " << i);
        }
      }
    }
}

void test_heuristics(int trials, mt19937& rng) {
  cout << "[3] heuristics: admissibility (and consistency) on every node" << endl;
  // (n, L, smax): small boxes give repeated coordinates, large smax gives long vectors
  const int cfg[][3] = { {6,6,2}, {7,10,3}, {8,12,3}, {8,8,5}, {9,15,4}, {6,20,6}, {10,10,4}, {7,5,7} };
  HStats poor, vtsp, paper;
  int before = failures;
  for (auto& c : cfg) {
    for (int t = 0; t < trials; t++) {
      Instance in = random_instance(rng, c[0], c[1], c[2], false);
      int n = c[0];
      UniProfile X(in.xs, in.univecs), Y(in.ys, in.univecs);

      vector<vector<int>> ctg(n);                // exact time left
      for (int i = 0; i < n; i++) ctg[i].assign(in.bivecs[i].size(), INF_COST);
      ctg[n-1][0] = 0;
      for (int i = n - 2; i >= 0; i--)
        for (size_t a = 0; a < in.bivecs[i].size(); a++)
          for (size_t b = 0; b < in.bivecs[i+1].size(); b++)
            ctg[i][a] = min(ctg[i][a], edge_weight(X, Y, in, i, a, b) + ctg[i+1][b]);

      auto hp = compute_heuristics(Heuristic::Poor,      in.bivecs, in.univecs, in.xs, in.ys);
      auto hv = compute_heuristics(Heuristic::VectorTSP, in.bivecs, in.univecs, in.xs, in.ys);
      check_heuristic("poor",      in, X, Y, hp, ctg, poor, true, true);
      // vectortsp must be admissible; it is not always consistent (A* re-opens nodes then)
      check_heuristic("vectortsp", in, X, Y, hv, ctg, vtsp, true, false);

      // the bound of the paper (theta = 0), only to show what happens (not asserted)
      VectorTSPHeuristic raw(in.xs, in.ys, 0);
      vector<vector<int>> hr(n);
      for (int i = 0; i < n; i++) for (auto& b : in.bivecs[i]) hr[i].push_back(raw.estimate(b, in.univecs, i));
      check_heuristic("paper", in, X, Y, hr, ctg, paper, false, false);
    }
  }
  auto line = [](const char* name, const HStats& s) {
    cout << "    " << name << ": nodes=" << s.nodes << " overestimated=" << s.inadmissible
         << " edges=" << s.edges << " inconsistent=" << s.inconsistent << endl;
  };
  line("poor                        ", poor);
  line("vectortsp (theta = smax)    ", vtsp);
  line("vectortsp (theta = 0, paper)", paper);
  cout << "    -> " << (failures == before ? "ok" : "FAILED") << endl;
}


/* ---------------- 4. solvers ---------------- */
// Reference: textbook Dijkstra with std::priority_queue and lazy deletion
int reference_dijkstra(const UniProfile& X, const UniProfile& Y, const Instance& in) {
  int n = in.xs.size();
  vector<vector<int>> dist(n);
  for (int i = 0; i < n; i++) dist[i].assign(in.bivecs[i].size(), INF_COST);
  using T = tuple<int,int,int>;
  priority_queue<T, vector<T>, greater<T>> pq;
  dist[0][0] = 0; pq.push({0, 0, 0});
  while (!pq.empty()) {
    auto [d, i, a] = pq.top(); pq.pop();
    if (d > dist[i][a]) continue;
    if (i == n - 1) return d;
    for (size_t b = 0; b < in.bivecs[i+1].size(); b++) {
      int nd = d + edge_weight(X, Y, in, i, a, b);
      if (nd < dist[i+1][b]) { dist[i+1][b] = nd; pq.push({nd, i+1, (int)b}); }
    }
  }
  return INF_COST;
}

void test_solvers(int trials, mt19937& rng) {
  cout << "[4] solvers vs reference Dijkstra ... " << flush;
  int before = failures;
  const int cfg[][3] = { {5,6,2}, {6,8,3}, {8,10,4}, {9,14,5}, {7,5,6}, {12,12,4} };
  long reopened = 0, instances = 0;
  for (auto& c : cfg)
    for (int t = 0; t < trials; t++) {
      for (int filter = 0; filter < 2; filter++) {
        Instance in = random_instance(rng, c[0], c[1], c[2], filter);
        UniProfile X(in.xs, in.univecs), Y(in.ys, in.univecs);
        int ref = reference_dijkstra(X, Y, in);
        auto hp = compute_heuristics(Heuristic::Poor,      in.bivecs, in.univecs, in.xs, in.ys);
        auto hv = compute_heuristics(Heuristic::VectorTSP, in.bivecs, in.univecs, in.xs, in.ys);
        SearchStats s1, s2, s3;
        int dj = solvedijkstra(X, Y, in.bivecs, &s1);
        int ap = solveastarfh(X, Y, in.bivecs, hp, &s2);
        int av = solveastarfh(X, Y, in.bivecs, hv, &s3);
        int dp = solve(X, Y, in.bivecs, in.univecs, in.xs, in.ys);
        CHECK(dj == ref && ap == ref && av == ref && dp == ref,
              "n=" << c[0] << " L=" << c[1] << " smax=" << c[2] << " filter=" << filter
              << ": ref=" << ref << " dijkstra=" << dj << " astar-poor=" << ap
              << " astar-vectortsp=" << av << " dp=" << dp);
        reopened += s1.reopened + s2.reopened + s3.reopened;
        instances++;
      }
    }
  cout << (failures == before ? "ok" : "FAILED") << "  (" << instances << " instances, "
       << reopened << " re-opened nodes)" << endl;
}


int main(int argc, const char* argv[]) {
  int trials = (argc > 1) ? atoi(argv[1]) : 200;
  unsigned seed = (argc > 2) ? atoi(argv[2]) : 1;
  mt19937 rng(seed);

  test_heap(trials, rng);
  test_intervals();
  test_heuristics(max(1, trials / 10), rng);
  test_solvers(max(1, trials / 10), rng);

  cout << (failures ? "SOME TESTS FAILED" : "ALL TESTS PASSED") << endl;
  return failures ? 1 : 0;
}
