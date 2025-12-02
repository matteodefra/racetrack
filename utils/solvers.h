#include "profiles.cpp"
#include "fibonacciheap.h"


int solvedijkstra(UniProfile& X, UniProfile& Y, vector<vector<BiVec>>& bivecs){
  int n = X.n;

  FibHeap heap;

  Handle h1 = heap.push(Conff {0, 0, {0,0}, 0, 0, 0});

  vector<vector<Handle>> handles(n);

  for (int i=0; i<n; i++) {
    handles[i].assign( bivecs[i].size(), NOHANDLE );
  }
  
  handles[0][0] = h1;
  (*h1).handle = h1;

  int mincost = __INT_MAX__;
  int heurcost = 0;

  while (!heap.empty()) {
    auto from = heap.top().config; 
    heap.pop();
    int currentcost = from.cost;

    handles[from.city][from.bivec_idx] = NOHANDLE;
    
    // cout << "Popping: " << from << endl;
    if (from.city == n-1){
      // cout << "Final cost: " << from.cost << endl;
      mincost = from.cost;
      return from.cost;
      // if (from.cost <= mincost) mincost = from.cost;
    } 
    else {
      for (int j=0; j<bivecs[from.city+1].size(); j++) {

        auto bivec = bivecs[from.city+1][j];
        Interval interv_x = X.get(from.city, from.bivec.ind_x, bivec.ind_x);
        Interval interv_y = Y.get(from.city, from.bivec.ind_y, bivec.ind_y);
        int tau = min_value(interv_x, interv_y);
        
        int distance = from.distance + tau;

        Conff to = {from.city+1, j, bivec, distance, distance, heurcost};
        
        if (handles[from.city+1][j] == NOHANDLE) {
          Handle h = heap.push(to);
          handles[from.city+1][j] = h;
          (*h).handle = h;
        } else {
          if ( to.cost < (*handles[from.city+1][j]).config.cost ) {
            (*handles[from.city+1][j]).config.cost = to.cost;
            (*handles[from.city+1][j]).config.distance = to.distance;
            heap.update( handles[from.city+1][j] );
          }
        }

      }
    }

  }
  return mincost;

}


int solveastarfh(UniProfile& X, UniProfile& Y, vector<vector<BiVec>>& bivecs, vector<vector<int>> &heuristics){
  int n = X.n;
  FibHeap heap;

  Handle h1 = heap.push(Conff {0, 0, {0,0}, 0 + heuristics[0][0], 0, heuristics[0][0]});

  vector<vector<Handle>> handles(n);

  for (int i=0; i<n; i++) {
    handles[i].assign( bivecs[i].size(), NOHANDLE );
  }
  
  handles[0][0] = h1;
  (*h1).handle = h1;

  int mincost = __INT_MAX__;
  // int heurcost = 0;

  while (!heap.empty()) {
    auto from = heap.top().config; 
    heap.pop();
    int currentcost = from.cost;

    handles[from.city][from.bivec_idx] = NOHANDLE;
    
    // cout << "Popping: " << from << endl;
    if (from.city == n-1){
      // cout << "Final cost: " << from.distance << endl;
      if (mincost <= from.distance) mincost = from.distance;
      return from.distance;
    } else {

      for (int j=0; j<bivecs[from.city+1].size(); j++) {

        auto bivec = bivecs[from.city+1][j];
        Interval interv_x = X.get(from.city, from.bivec.ind_x, bivec.ind_x);
        Interval interv_y = Y.get(from.city, from.bivec.ind_y, bivec.ind_y);
        int tau = min_value(interv_x, interv_y);

        // h(city+1 [j])
        int heurcost = heuristics[from.city+1][j];
        
        // g(city+1 [j])
        int distance = from.distance + tau;

        Conff to = {from.city+1, j, bivec, distance + heurcost, distance, heurcost};


        if (handles[from.city+1][j] == NOHANDLE) {
          Handle h = heap.push(to);
          handles[from.city+1][j] = h;
          (*h).handle = h;
        } else {
          if (to.distance < (*handles[from.city+1][j]).config.distance) {
            (*handles[from.city+1][j]).config.distance = to.distance;
            (*handles[from.city+1][j]).config.cost = to.cost;
            heap.update( handles[from.city+1][j] );
          }
        }

        
      }
    }
  }

  return mincost;

}


int solve(const UniProfile& X, const UniProfile& Y, 
  const vector<vector<BiVec>>& bivecs, std::vector<UniVec>& univecs, 
  const vector<int>& xs, const vector<int>& ys){

  int n = X.n;

  int number_bivecs = 0;
  for (int i=0; i<n; i++){
    number_bivecs += bivecs[i].size();
  }
  vector<int> costs(number_bivecs, __INT_MAX__);
  for (int i=0; i<bivecs[0].size();i++){
    costs[i] = 0;
  }
  // vector<TruePoint> points(n, {0,0,0,0});

  auto cost_from_it = costs.begin();
  for (int i = 0; i < n-1; i++){
    auto cost_to_it = cost_from_it + bivecs[i].size();
    for (auto bivec_from : bivecs[i]){
      for (auto bivec_to : bivecs[i+1]){
        // if (i+1==n-1) bivec_to.print(univecs);
        auto x_interv = X.get(i, bivec_from.ind_x, bivec_to.ind_x);
        auto y_interv = Y.get(i, bivec_from.ind_y, bivec_to.ind_y);
        int val = min_value(x_interv, y_interv);
        int total = *cost_from_it + val;
        if (total < *cost_to_it){
          *cost_to_it = total;
        }
        cost_to_it++;
      }
      cost_from_it++;
      cost_to_it -= bivecs[i+1].size();
    }

  // I know at the end of this for loop the costs in all potatoes
  }

  return costs[costs.size()-1];
}