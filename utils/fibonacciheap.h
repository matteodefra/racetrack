/**
 * Fibonacci Heap with Boost
 */
struct Conff{
  int city;
  int bivec_idx;
  BiVec bivec;
  int cost;
  int distance;
  int heurcost;
  

  friend std::ostream& operator <<(std::ostream& os, Conff const& in) {
    return os << "City " << in.city << " crossing with bivec (" << in.bivec.ind_x << "," << in.bivec.ind_y << ") and cost " << in.cost << " and true distance " << in.distance << " and heuristic " << in.heurcost;
  }

};


struct heapdata;
using FibHeap = boost::heap::fibonacci_heap<heapdata>;
using Handle  = FibHeap::handle_type;
static const Handle NOHANDLE{}; 

struct heapdata {

  Conff config;
  Handle handle;

  heapdata(Conff conf) : config(conf), handle() {};
 
  bool operator<(heapdata const &b) const {
    if (config.cost != b.config.cost) return config.cost > b.config.cost;        // smaller cost -> higher priority
    if (config.city != b.config.city) return config.city > b.config.city;        // tie-break on city
    return config.bivec_idx > b.config.bivec_idx;
  }

};
