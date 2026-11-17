#pragma once
/**
 * Heuristics h(city, bivec) for A*: lower bounds on the time left after the
 * vehicle crossed `city` with the given bi-vector, until it stops in the last city.
 *
 *   Heuristic::Zero      h = 0                      (Dijkstra)
 *   Heuristic::Poor      time to the LAST city only (insights.pdf, Bekos et al.):
 *                        the intermediate cities are ignored.
 *   Heuristic::VectorTSP per-axis bound that respects the visiting order
 *                        (VectorTSP paper, arXiv 2006.03666, Sec. 5.2).
 *
 * Notation (1D, one axis).
 *   After crossing city i with the uni-vector (from,to) the vehicle is at
 *   position P = c[i] + to with speed v = to - from.  Time is the number of moves.
 *
 * IMPORTANT (admissibility): the bound of the paper assumes that, to go from a
 * city "on the right" to a city "on the left", the vehicle must stop and turn.
 * In this model a vector can cross several cities, so if the two cities are
 * closer than the longest vector (smax) one move can cover both and no turn is
 * needed.  VectorTSPHeuristic therefore only counts a turn between two extreme
 * cities that are MORE than `theta` apart (theta = smax, the longest vector of
 * the instance).  With theta = 0 it is the bound of the paper; it can then
 * overestimate the true time (A* may return a non-optimal value).
 */
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <vector>
#include "profiles.cpp"

enum class Heuristic { Zero, Poor, VectorTSP };

inline const char* heuristic_name(Heuristic kind) {
  switch (kind) {
    case Heuristic::Zero:      return "zero";
    case Heuristic::Poor:      return "poor";
    case Heuristic::VectorTSP: return "vectortsp";
  }
  return "?";
}


/* ------------------------------------------------------------------ */
/* 1D building blocks                                                  */
/* ------------------------------------------------------------------ */

// Smallest number of moves to go from (x, speed s) to (xp, speed sp) on one axis.
// compute_interval returns [t1,t2] U [t3,inf) (t1 = -1 if the first part is empty),
// so the SMALLEST feasible time is t1 if it exists, t3 otherwise.
// (Using t3 alone can overestimate the time: it is a bug for a heuristic.)
inline int min_time_1d(int x, int s, int xp, int sp) {
  Interval interv = compute_interval(x, s, xp, sp);
  return (interv.t1 != -1) ? interv.t1 : interv.t3;
}

// Smallest number of moves to go from rest to rest over a distance d >= 0.
// In t moves the largest distance is floor(t*t/4), so this is ceil(2*sqrt(d)).
inline int rest_to_rest_time(int d) {
  int t = 0;
  while ((long long)t * t / 4 < d) t++;
  return t;
}

// Longest vector of the instance: crossing a city with one move never uses more.
inline int longest_vector(const std::vector<UniVec>& univecs) {
  int smax = 0;
  for (const UniVec& u : univecs) smax = std::max(smax, std::abs(u.to - u.from));
  return smax;
}


/* ------------------------------------------------------------------ */
/* Poor heuristic: time to the last city, intermediate cities ignored  */
/* ------------------------------------------------------------------ */

// Exact smallest time (both axes together) to reach the last city and stop there,
// starting from the state of `bivec` in `city`.
inline int computeheur_poor(const BiVec& bivec, const std::vector<UniVec>& univecs,
                            const std::vector<int>& xs, const std::vector<int>& ys, int city) {
  const int n = xs.size();
  if (city == n - 1) return 0;

  const UniVec& ux = univecs[bivec.ind_x];
  const UniVec& uy = univecs[bivec.ind_y];
  Interval ix = compute_interval(xs[city] + ux.to, ux.to - ux.from, xs[n - 1], 0);
  Interval iy = compute_interval(ys[city] + uy.to, uy.to - uy.from, ys[n - 1], 0);
  return min_value(ix, iy);
}


/* ------------------------------------------------------------------ */
/* VectorTSP heuristic: per axis, respecting the visiting order        */
/* ------------------------------------------------------------------ */

/**
 * Plan of one axis for the cities after `city`.
 *
 * The remaining cities c[city+1..n-1] are reduced to the "forced turns":
 * extreme positions u_1, u_2, ..., u_k (alternately a maximum and a minimum,
 * consecutive ones more than theta apart) after which the vehicle must stop and
 * go back, plus the last city. Time left >= (time to stop at u_1 or beyond)
 *                                         + sum of rest-to-rest times between turns.
 * Everything that does not depend on the current speed is computed here once.
 */
struct AxisPlan {
  int dir = 0;         // +1: first turn is a maximum, -1: a minimum, 0: no forced turn
  int first_turn = 0;  // u_1
  int tail_cost = 0;   // rest-to-rest times: u_1 -> u_2 -> ... -> u_k -> last city
};

inline AxisPlan make_axis_plan(const std::vector<int>& c, int city, int theta) {
  const int n = c.size();
  std::vector<int> turns;   // u_1 ... u_k

  // Zig-zag filter on c[city+1 .. n-2]: a turn is confirmed only when the
  // position comes back by more than theta. (Skipping cities only weakens the bound.)
  int start = c[city];
  int dir = 0;          // direction of the current run
  int extreme = start;  // furthest position of the current run
  for (int j = city + 1; j < n - 1; j++) {
    int x = c[j];
    if (dir == 0) {
      if (std::abs(x - start) > theta) { dir = (x > start) ? 1 : -1; extreme = x; }
    } else if (dir > 0) {
      if (x >= extreme)               extreme = x;
      else if (extreme - x > theta) { turns.push_back(extreme); dir = -1; extreme = x; }
    } else {
      if (x <= extreme)               extreme = x;
      else if (x - extreme > theta) { turns.push_back(extreme); dir = 1; extreme = x; }
    }
  }
  // The last city is the end of the tour: the vehicle must stop there, so the
  // last run needs no distance condition.
  const int last = c[n - 1];
  if ((dir > 0 && last < extreme) || (dir < 0 && last > extreme)) turns.push_back(extreme);

  AxisPlan plan;
  if (turns.empty()) return plan;

  plan.first_turn = turns[0];
  plan.dir = (turns[0] > start) ? 1 : -1;
  for (size_t k = 0; k + 1 < turns.size(); k++)
    plan.tail_cost += rest_to_rest_time(std::abs(turns[k] - turns[k + 1]));
  plan.tail_cost += rest_to_rest_time(std::abs(turns.back() - last));
  return plan;
}

// Lower bound of one axis, from position P with speed v.
inline int axis_bound(const AxisPlan& plan, int P, int v, int last) {
  // 1) go straight to the last city (this alone is the "poor" bound)
  int bound = min_time_1d(P, v, last, 0);
  if (plan.dir == 0) return bound;

  // 2) first stop at the first turn (or beyond), then follow the plan.
  //    Braking at once from speed v needs |v| moves and ends at position q.
  int av = std::abs(v);
  int brake = av * (av - 1) / 2;
  int q = (v > 0) ? P + brake : P - brake;

  bool brake_is_enough = (plan.dir > 0) ? (q >= plan.first_turn) : (q <= plan.first_turn);
  // Beyond the braking point the time to stop at p grows with the distance,
  // so the nearest admissible stop (the turn itself) is the cheapest.
  int first = brake_is_enough ? av : min_time_1d(P, v, plan.first_turn, 0);

  return std::max(bound, first + plan.tail_cost);
}


class VectorTSPHeuristic {
 public:
  // theta: two extreme cities need a turn between them only if they are more than
  //        theta apart. Use longest_vector(univecs) for a bound that is always valid.
  VectorTSPHeuristic(const std::vector<int>& xs, const std::vector<int>& ys, int theta)
      : xs_(xs), ys_(ys) {
    const int n = xs.size();
    plan_x_.resize(n);
    plan_y_.resize(n);
    for (int i = 0; i < n - 1; i++) {
      plan_x_[i] = make_axis_plan(xs, i, theta);
      plan_y_[i] = make_axis_plan(ys, i, theta);
    }
  }

  int estimate(const BiVec& bivec, const std::vector<UniVec>& univecs, int city) const {
    const int n = xs_.size();
    if (city == n - 1) return 0;

    const UniVec& ux = univecs[bivec.ind_x];
    const UniVec& uy = univecs[bivec.ind_y];
    int hx = axis_bound(plan_x_[city], xs_[city] + ux.to, ux.to - ux.from, xs_[n - 1]);
    int hy = axis_bound(plan_y_[city], ys_[city] + uy.to, uy.to - uy.from, ys_[n - 1]);
    return std::max(hx, hy);
  }

 private:
  std::vector<int> xs_, ys_;
  std::vector<AxisPlan> plan_x_, plan_y_;
};


/* ------------------------------------------------------------------ */
/* Table heuristics[city][bivec_idx], parallel to bivecs               */
/* ------------------------------------------------------------------ */

inline std::vector<std::vector<int>> compute_heuristics(
    Heuristic kind, const std::vector<std::vector<BiVec>>& bivecs,
    const std::vector<UniVec>& univecs, const std::vector<int>& xs, const std::vector<int>& ys) {
  const int n = xs.size();
  std::vector<std::vector<int>> table(n);
  for (int i = 0; i < n; i++) table[i].assign(bivecs[i].size(), 0);
  if (kind == Heuristic::Zero) return table;

  if (kind == Heuristic::Poor) {
    for (int i = 0; i < n - 1; i++)
      for (size_t j = 0; j < bivecs[i].size(); j++)
        table[i][j] = computeheur_poor(bivecs[i][j], univecs, xs, ys, i);
  } else {
    VectorTSPHeuristic vtsp(xs, ys, longest_vector(univecs));
    for (int i = 0; i < n - 1; i++)
      for (size_t j = 0; j < bivecs[i].size(); j++)
        table[i][j] = vtsp.estimate(bivecs[i][j], univecs, i);
  }
  return table;
}
