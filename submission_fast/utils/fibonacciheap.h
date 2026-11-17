#pragma once
/**
 * Open list (Fibonacci heap, boost implementation) for Dijkstra / A* on the
 * layered graph "one layer per city, one node per bi-vector".
 *
 * A node of the search is identified by (city, bivec_idx).
 * The heap stores a small, trivially copyable record (Conff). Nothing else
 * (BiVec, handle, ...) is stored inside the heap: the BiVec is looked up in
 * bivecs[city][bivec_idx] and the handle of every node lives in a separate
 * table owned by the search (see solvers.h).
 */
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>

struct Conff {
  int city;       // layer of the node
  int bivec_idx;  // index of the node inside bivecs[city]
  int cost;       // f = g + h, the key of the heap
  int distance;   // g, best known time to reach the node
  int heurcost;   // h, estimate of the time left

  friend std::ostream& operator <<(std::ostream& os, Conff const& in) {
    return os << "City " << in.city << " with bivec #" << in.bivec_idx
              << ", cost " << in.cost << ", true distance " << in.distance
              << ", heuristic " << in.heurcost;
  }
};


struct heapdata {
  Conff config;

  heapdata(const Conff& conf) : config(conf) {}

  // boost::heap is a MAX-heap: "a < b" means "a has a LOWER priority than b".
  // We want the smallest f first, so the comparison is reversed.
  // Ties are broken in a fixed order (deeper city first, then smaller index),
  // so the search is deterministic.
  bool operator<(heapdata const &b) const {
    if (config.cost != b.config.cost)           return config.cost > b.config.cost;
    if (config.city != b.config.city)           return config.city < b.config.city;
    return config.bivec_idx > b.config.bivec_idx;
  }
};


using FibHeap = boost::heap::fibonacci_heap<heapdata>;
using Handle  = FibHeap::handle_type;

// A default-constructed handle: "this node is NOT in the heap".
static const Handle NOHANDLE{};
