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



auto bi_vectors_per_city(const vector<int>& xs, const vector<int>& ys, vector<UniVec> &univecs, bool heur, vector<vector<int>>& heuristics){
  int n = xs.size();

  auto bivecs = bi_vectors(univecs);
  vector<vector<BiVec>> result(n, vector<BiVec>());

  bool cond;

  int counter = 0;

  result[0].push_back({0,0});
  heuristics[0].push_back( computeheur( {0,0}, univecs, xs, ys, 0 ) );
  result[n-1].push_back({0,0});
  heuristics[n-1].push_back( 0 );
  for (int i=1; i<n-1; i++){
    for (auto bivec: bivecs){     
      // cout << cond << endl;
      counter += 1;
      result[i].push_back(bivec);
      if (heur) heuristics[i].push_back(computeheur( bivec, univecs, xs, ys, i ));
    }
  }
  //cout << "Bivecs used: " << counter << endl;

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

  double total_astar = 0.0;
  double total_dijkstra = 0.0;
  double total_ourdp = 0.0;

  bool heur = true;


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
  vector<vector<int>> heuristics(n, vector<int>());

  // smax from safe solutions
  auto univecs = uni_vectors(smax);
  UniProfile X(xs, univecs);
  UniProfile Y(ys, univecs);
  vector<vector<BiVec>> bivecs = bi_vectors_per_city(xs, ys, univecs, heur, heuristics);
  // cout << "profiles OK" << endl;

  int dijkstra_opt;
  long dijktime;

  int astar_opt;
  long astartime; 

  int dp_opt;
  long dptime; 

  // #pragma omp parallel num_threads(3)
  for (int alg=0; alg<=2; alg++) {
    if (alg==0) {
      // LAUNCHING DIJKSTRA
      {
        Utimer timer("Dijkstra Time", &dijktime);

        dijkstra_opt = solvedijkstra( X, Y, bivecs );
      }
    } else if (alg==1) {
      // LAUNCHING ASTAR
      {
        Utimer timer("A* Time", &astartime);

        astar_opt = solveastarfh( X, Y, bivecs, heuristics );
      }
    } else {
      // LAUNCHING OUR DP
      {
        Utimer timer("DP Time", &dptime);

        dp_opt = solve( X, Y, bivecs, univecs, xs, ys );
      } 
    }
  }
      
  // #pragma omp barrier
  if (dijkstra_opt==astar_opt && astar_opt==dp_opt ) {
    cout << "Final cost DJK: " << dijkstra_opt << endl;
    cout << "Final cost A*: " << astar_opt << endl;  
    cout << "Final cost DP: " << dp_opt << endl;
    
    cout << "Time djk: " << dijktime << endl;
    cout << "Time A*: " << astartime << endl;
    cout << "Time DP: " << dptime << endl;
  }

  return 0;

}
