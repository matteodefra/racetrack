#include <boost/heap/fibonacci_heap.hpp>
#include "utils/solvers.h"
#include "utils/pruning.h"
#include "utils/utimer.hpp"
#include "utils/heuristic.h"
#include <queue>
#include <set>
#include <iostream>
#include <fstream>


// Tangent-based
bool is_valid_inbox(const BiVec bivec, const vector<UniVec>& univecs, const vector<int>& xs, const vector<int>& ys, int city, vector<int> &boxlimits){
  // city location
  Point b = {(double)xs[city], (double)ys[city]};

  int x_from = univecs[bivec.ind_x].from;
  int y_from = univecs[bivec.ind_y].from;
  // vector tail
  Point d = {(double)(x_from + b.x), (double)(y_from + b.y)};
  
  int x_to = univecs[bivec.ind_x].to;
  int y_to = univecs[bivec.ind_y].to;
  // vector point
  Point e = {(double)(x_to + b.x), (double)(y_to + b.y)};
  // cout << d << " " << e << endl;

  return checklanding(b, d, e, boxlimits);

}



// All the bi-vectors in every city (no pruning in this scenario).
auto bi_vectors_per_city(const vector<UniVec> &univecs, int n){
  auto bivecs = bi_vectors(univecs);
  vector<vector<BiVec>> result(n, vector<BiVec>());

  result[0].push_back({0,0});      // start at rest
  result[n-1].push_back({0,0});    // end at rest
  for (int i=1; i<n-1; i++){
    for (auto bivec: bivecs){
      result[i].push_back(bivec);
    }
  }

  return result;
}



int main(int argc, const char * argv[]) {
  vector<int> xs, ys;
  int n;


  if (argc!=3) {
    cout << "Usage ./executable delta n" << endl;
    return 0;
  } 

  int delta = std::atoi(argv[1]); // first input is the number of experiments to run
  int numcities = std::atoi(argv[2]); // second input is the number of cities to visit

  xs.push_back(0); ys.push_back(0);

  for (int i=1; i<numcities; i++) {
    xs.push_back(i*delta);
    ys.push_back(-i);
  }
  n = xs.size();

  for (int i=0; i<n; i++) {
    std::cout << "(x,y)=(" << xs[i] << "," << ys[i] << ")" << std::endl;
  }

  int D = floor(sqrt(distance(xs[0], ys[0], xs[n-1], ys[n-1])));

  // int smax = compute_smax( xs, ys );
  int smax = floor(D/delta)*delta;
  cout << smax << endl;
  // int sqrtD = smax;
  vector<int> boxlimits = boundingbox_limits( xs, ys );
  for (auto x: boxlimits) {
    cout << x << endl;
  }
  // smax from safe solutions
  auto univecs = uni_vectors(smax);
  UniProfile X(xs, univecs);
  UniProfile Y(ys, univecs);
  vector<vector<BiVec>> bivecs = bi_vectors_per_city(univecs, n);
  // cout << "profiles OK" << endl;

  vector<vector<int>> h_poor = compute_heuristics(Heuristic::Poor,      bivecs, univecs, xs, ys);
  vector<vector<int>> h_vtsp = compute_heuristics(Heuristic::VectorTSP, bivecs, univecs, xs, ys);

  int dijkstra_opt, poor_opt, vtsp_opt, dp_opt;
  long dijktime, poortime, vtsptime, dptime;
  SearchStats djk_stats, poor_stats, vtsp_stats;

  {
    Utimer timer("Dijkstra Time", &dijktime);
    dijkstra_opt = solvedijkstra( X, Y, bivecs, &djk_stats );
  }
  {
    Utimer timer("A* (poor) Time", &poortime);
    poor_opt = solveastarfh( X, Y, bivecs, h_poor, &poor_stats );
  }
  {
    Utimer timer("A* (VectorTSP) Time", &vtsptime);
    vtsp_opt = solveastarfh( X, Y, bivecs, h_vtsp, &vtsp_stats );
  }
  {
    Utimer timer("DP Time", &dptime);
    dp_opt = solve( X, Y, bivecs, univecs, xs, ys );
  }

  cout << "Final cost DJK: " << dijkstra_opt << endl;
  cout << "Final cost A* (poor): " << poor_opt << endl;
  cout << "Final cost A* (VectorTSP): " << vtsp_opt << endl;
  cout << "Final cost DP: " << dp_opt << endl;

  cout << "Time djk: " << dijktime << " (expanded " << djk_stats.popped << ")" << endl;
  cout << "Time A* (poor): " << poortime << " (expanded " << poor_stats.popped << ")" << endl;
  cout << "Time A* (VectorTSP): " << vtsptime << " (expanded " << vtsp_stats.popped << ")" << endl;
  cout << "Time DP: " << dptime << endl;

  if (!(dijkstra_opt==poor_opt && poor_opt==vtsp_opt && vtsp_opt==dp_opt)) {
    cout << "MISMATCH between the algorithms" << endl;
    return 1;
  }

  return 0;

}
