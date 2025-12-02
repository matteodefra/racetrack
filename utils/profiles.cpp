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


int min_value(Interval x, Interval y){
  if (x.t1 == -1 && y.t1 == -1){
    return max(x.t3, y.t3);
  } else {
    bool len4_x = (x.t1 != -1);
    bool len4_y = (y.t1 != -1);

    int tau = -1;

    if (len4_x && len4_y) {
      vector<LInterv> A =  { {x.t1, x.t2}, {x.t3, __INT_MAX__} };
      vector<LInterv> B =  { {y.t1, y.t2}, {y.t3, __INT_MAX__} };
      tau = minimumIntersection(A, B);
    } else if (len4_x && !len4_y) {
      vector<LInterv> A =  { {x.t1, x.t2}, {x.t3, __INT_MAX__} };
      vector<LInterv> B =  { {y.t3, __INT_MAX__} };
      tau = minimumIntersection(A, B);
    } else if (!len4_x && len4_y) {
      vector<LInterv> A =  { {x.t3, __INT_MAX__} };
      vector<LInterv> B =  { {y.t1, y.t2}, {y.t3, __INT_MAX__} };
      tau = minimumIntersection(A, B);
    } else {
      tau = max( x.t3, y.t3 ); 
    }
    return tau;
  }
}

