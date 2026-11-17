#pragma once
#include <iostream>
#include "intervals.cpp"
#include "vectors.cpp"

using namespace std;


struct TruePoint {
  int x;
  int y;
  int xdx;
  int ydy;

  friend std::ostream& operator <<(std::ostream& os, TruePoint const& p) {
    return os << "Passing with " << p.x << "," << p.y << " " << p.xdx << "," << p.ydy;
  }
};

struct UniProfile{
  int n;
  int nb_vecs;
  vector<Interval> intervals;
  // recover the corresponding real coordinates
 
  UniProfile(vector<int> cities, const vector<UniVec>& vecs){
    n = cities.size();
    nb_vecs = vecs.size();
    intervals.reserve(2*nb_vecs + nb_vecs*nb_vecs*(n-3) + 2);
 
    for (auto vec : vecs){
      Interval interv = compute_interval(cities[0], 0, cities[1]+vec.to, vec.to-vec.from);
      intervals.push_back(interv);
    }
    for (int i=1; i<n-2; i++){
      int x1 = cities[i];
      int x2 = cities[i+1];
      for (auto v1 : vecs){
        for (auto v2 : vecs){
          Interval interv = compute_interval(x1+v1.to, v1.to-v1.from, x2+v2.to, v2.to-v2.from);
          intervals.push_back(interv);
        }
      }
    }
    for (auto v : vecs){
      Interval interv = compute_interval(cities[n-2]+v.to, v.to-v.from, cities[n-1], 0);
      intervals.push_back(interv);
    }
  }

  Interval get(int index) const{
    return intervals[index];
  }
  Interval get(int city_from, int from_ind, int to_ind) const{ // vec_to crosses next city (implicit)
    if (city_from == 0){
      return intervals[to_ind];
    } else if (city_from == n-2){
      return intervals[nb_vecs + (city_from-1)*nb_vecs*nb_vecs + from_ind];
    } else {
      return intervals[nb_vecs + (city_from-1)*nb_vecs*nb_vecs + from_ind*nb_vecs + to_ind];
    }
  }
  // All the intervals of (city_from, from_ind) at once:
  //   get(city_from, from_ind, to_ind) == row[to_ind * stride]
  // (stride = 0 for the last pair of cities: the target state is fixed, stop in the last city).
  const Interval* row(int city_from, int from_ind, int& stride) const{
    if (city_from == 0)         { stride = 1; return intervals.data(); }
    else if (city_from == n-2)  { stride = 0; return intervals.data() + nb_vecs + (size_t)(city_from-1)*nb_vecs*nb_vecs + from_ind; }
    else                        { stride = 1; return intervals.data() + nb_vecs + (size_t)(city_from-1)*nb_vecs*nb_vecs + (size_t)from_ind*nb_vecs; }
  }
};


struct LInterv {
  int L,R;
};

// Compute intersection of two intervals.
// Returns empty vector if they don't intersect.
bool intersect(LInterv x, LInterv y, LInterv &out) {
  int L = max(x.L, y.L);
  int R = min(x.R, y.R);

  if (L <= R) {
      out = {L, R};
      return true;
  }
  return false;
}

int minimumIntersection(
    vector<LInterv>& A,
    vector<LInterv>& B) {
    int best = __INT_MAX__;

    for (const auto& a : A) {
        for (const auto& b : B) {
            LInterv opt;
            if (intersect(a, b, opt)) {
                best = min(best, opt.L);  // take minimum starting point
            }
        }
    }
    return best;
}


// Smallest t in ([x.t1,x.t2] U [x.t3,inf)) n ([y.t1,y.t2] U [y.t3,inf)).
// (t1 == -1: the first part is empty.) No heap allocation: this is called once per edge.
inline int min_value(Interval x, Interval y){
  int best = max(x.t3, y.t3);                                     // tail n tail
  if (x.t1 != -1) {
    int L = max(x.t1, y.t3);              if (L <= x.t2 && L < best) best = L;   // [x1,x2] n tail_y
    if (y.t1 != -1) { L = max(x.t1, y.t1); if (L <= min(x.t2, y.t2) && L < best) best = L; }
  }
  if (y.t1 != -1) { int L = max(y.t1, x.t3); if (L <= y.t2 && L < best) best = L; } // tail_x n [y1,y2]
  return best;
}
