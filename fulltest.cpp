#include <boost/heap/fibonacci_heap.hpp>
#include "utils/solvers.h"
#include "utils/utimer.hpp"
#include "utils/heuristic.h"
#include "utils/pruning.h"
// #include <queue>
// #include <set>
// #include <map>
// #include <numeric>
// #include <iostream>
// #include <fstream>
#include <mutex>
#include <omp.h>
#include <sstream>
// #include <string>

using namespace std;

auto bi_vectors_per_city(vector<int>& xs, vector<int>& ys, vector<UniVec> &univecs, bool heur, vector<vector<int>>& heuristics, bool filter, int S){
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
      if (filter) {
        cond = filteringS(bivec, univecs, xs, ys, i, S);
      }
      if (!cond) {
        counter += 1;
        result[i].push_back(bivec);
        if (heur) heuristics[i].push_back(computeheur( bivec, univecs, xs, ys, i ));
      }
      
    }
  }
  //cout << "Bivecs used: " << counter << endl;

  return result;
}


int trivialsolution(vector<int> &xs, vector<int> &ys ) {
  
  int n = xs.size();

  vector<vector<int>> heuristics(n, vector<int>());
  
  auto univecs = uni_vectors(1);
  UniProfile X(xs, univecs);
  UniProfile Y(ys, univecs);
  vector<vector<BiVec>> bivecs = bi_vectors_per_city(xs, ys, univecs, false, heuristics, false, 0);

  int S = solve(X, Y, bivecs, univecs, xs, ys);
  //cout << "First trivial solution: " << S << endl;
  return S;
}


int main(int argc, const char * argv[]) {

  if (argc!=4) {
    cout << "Usage ./executable #trials n L" << endl;
    return 0;
  } 

  int trials = std::atoi(argv[1]); // first input is the number of experiments to run
  int n = std::atoi(argv[2]); // second input is the number of cities to visit
  int L = std::atoi(argv[3]); // third input is the dimension of the box

  mutex update;
  double total_astar = 0.0;
  double total_dijkstra = 0.0;
  double total_ourdp = 0.0;

  bool heur = true;
  
  vector<vector<int>> allxs(trials);
  vector<vector<int>> allys(trials);

  for (int t=0; t<trials; t++) {
    // vector<int> xs, ys;
    for (int i = 0; i<n; i++){
      int x = int(rand()%L);
      int y = int(rand()%L);
      allxs[t].push_back(x); allys[t].push_back(y);
    }
    allxs[t][n-1] = allxs[t][0]; allys[t][n-1] = allys[t][0]; // setting last equal to first
  }

  // omp_set_max_active_levels(2);
  
  // #pragma omp parallel num_threads(trials)
  for (int t=0; t<trials; t++) {

    vector<int> xs = allxs[t];
    vector<int> ys = allys[t];

    vector<vector<int>> heuristics(n, vector<int>());

    // smax heuristic (sqrt(D))
    // int smax = compute_smax( xs, ys );
    
    // smax from safe solutions
    int S = trivialsolution( xs, ys );
    int smax = S/4;

    auto univecs = uni_vectors(smax);
    UniProfile X(xs, univecs);
    UniProfile Y(ys, univecs);
    vector<vector<BiVec>> bivecs = bi_vectors_per_city(xs, ys, univecs, heur, heuristics, true, S);
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
      {
        // lock_guard<mutex> lg(update);
        if (t==0) {
          total_dijkstra += dijktime;
          total_astar += astartime;
          total_ourdp += dptime;
        } else {
          total_dijkstra += (dijktime)/2;
          total_astar += (astartime)/2;
          total_ourdp += (dptime)/2;
        }
      }
    }


    cout << "Final cost DJK: " << dijkstra_opt << endl;

    

    cout << "Final cost A*: " << astar_opt << endl;

       
    
    cout << "Final cost DP: " << dp_opt << endl;


  }

  std::stringstream ss;
  ss << "n=" << n << "_t=" << trials << "_L=" << L;

  std::stringstream res;
  res << "Dijk " << total_dijkstra << endl;
  res << "Astar " << total_astar << endl;
  res << "DP " << total_ourdp << endl;

  return 0;
}
