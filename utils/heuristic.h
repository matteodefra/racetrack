auto computetime( BiVec bivec, const vector<UniVec>& univecs, const vector<int>& cities, int current_city, int isx ) {

  int n = cities.size();

  int city = cities[current_city];
  
  int j = current_city+1;

  int total = 0;

  int landingpoint = 0;

  int direction = -1;
  int last = city;
  int initial = city;
  int speed = 0;
  if (isx) speed = univecs[bivec.ind_x].to-univecs[bivec.ind_x].from;
  else speed = univecs[bivec.ind_y].to-univecs[bivec.ind_y].from;

  int headvector = 0;
  if (isx) headvector = univecs[bivec.ind_x].to;
  else headvector = univecs[bivec.ind_y].to;

  // std::cout << "Starting from city " << current_city <<  " with pos " << city << std::endl;
  // std::cout << "Speed vector " << speed << std::endl;

  bool firsttime = false;
  while (j<n-1) {
    int jx = cities[j];
    // std::cout << "next pos " << jx << std::endl;
    if (direction==-1) {
      if (jx >= city) { 
        // std::cout << "first city greater setting direction 1" << std::endl;
        last = jx;
        direction = 1;
      } else { 
        // std::cout << "first city smaller setting direction 0" << std::endl;
        direction=0;
        last = jx;
      }
      // j += 1;
    } else {
      if (jx==last) {
        last = jx;
      } else if (jx > last && direction==1) {
        // std::cout << "continuing good direction 1" << std::endl;
        last = jx;
      } else if (jx > last && direction==0) {
        // std::cout << "reversing direction" << std::endl;
        // accumulate cost
        direction=1;
        if (!firsttime) {
          firsttime = true;
          // compute the landing point
          int dist = (abs(speed)*(abs(speed)-1))/2;
          if (speed<0) landingpoint = (initial+headvector) - dist;
          else landingpoint = (initial+headvector) + dist;

          // the landing point is far away then the reversing point => 
          if (landingpoint < last) {
            Interval interv = compute_interval( initial+headvector, speed, landingpoint, 0 );
            total += interv.t3;
            initial = landingpoint;
          } else {
            Interval interv = compute_interval( initial+headvector, speed, last, 0 );
            total += interv.t3;
            initial = last;
          }
          // total += (interv.t1==-1) ? interv.t3 : interv.t1;
        } else {
          Interval interv = compute_interval( initial, 0, last, 0 );
          // total += (interv.t1==-1) ? interv.t3 : interv.t1;
          total += interv.t3;
          initial = last;
        }
        // std::cout << "Accumulated cost up to now " << total << std::endl; 
        last=jx;
      } else if (jx<last && direction==0) {
        // std::cout << "continuing good direction 0" << std::endl;
        last = jx;
      } else if (jx<last && direction==1) {
        // std::cout << "reversing direction" << std::endl;
        // accumulate cost
        direction=0;
        if (!firsttime) {
          firsttime = true;
          // compute the landing point
          int dist = (abs(speed)*(abs(speed)-1))/2;
          if (speed<0) landingpoint = (initial+headvector) - dist;
          else landingpoint = (initial+headvector) + dist;

          // the landing point is far away then the reversing point => 
          if (landingpoint >= last) {
            Interval interv = compute_interval( initial+headvector, speed, landingpoint, 0 );
            total += interv.t3;
            initial = landingpoint;
          } else {
            Interval interv = compute_interval( initial+headvector, speed, last, 0 );
            total += interv.t3;
            initial = last;
          }
        } else {
          Interval interv = compute_interval( initial, 0, last, 0 );
          // total += (interv.t1==-1) ? interv.t3 : interv.t1;
          total += interv.t3;
          initial = last;
        } 
        // std::cout << "Accumulated cost up to now " << total << std::endl; 
        last=jx;
      }

    }

    j+=1;
    
  }
  // last city
  int jx = cities[j];
  if (jx==last) {
    if (!firsttime) {
      firsttime = true;
      Interval interv = compute_interval( initial+headvector, speed, last, 0 );
      // total += (interv.t1==-1) ? interv.t3 : interv.t1;
      total += interv.t3;
    } else {
      Interval interv = compute_interval( initial, 0, last, 0 );
      // total += (interv.t1==-1) ? interv.t3 : interv.t1;
      total += interv.t3;
    } 
  }
  else if (jx > last && direction==1) {
    last = jx;
    if (!firsttime) {
      firsttime = true;
      Interval interv = compute_interval( initial+headvector, speed, last, 0 );
      // total += (interv.t1==-1) ? interv.t3 : interv.t1;
      total += interv.t3;
    } else {
      Interval interv = compute_interval( initial, 0, last, 0 );
      // total += (interv.t1==-1) ? interv.t3 : interv.t1;
      total += interv.t3;
    } 
    // std::cout << "Accumulated cost up to now " << total << std::endl; 
  } else if (jx > last && direction==0) {
    // accumulate cost
    direction=1;
    if (!firsttime) {
      firsttime = true;
      // compute the landing point
      int dist = (abs(speed)*(abs(speed)-1))/2;
      if (speed<0) landingpoint = (initial+headvector) - dist;
      else landingpoint = (initial+headvector) + dist;

      // the landing point is far away then the reversing point => 
      if (landingpoint < last) {
        Interval interv = compute_interval( initial+headvector, speed, landingpoint, 0 );
        total += interv.t3;
        initial = landingpoint;
      } else {
        Interval interv = compute_interval( initial+headvector, speed, last, 0 );
        total += interv.t3;
        initial = last;
      }
    } else {
      Interval interv = compute_interval( initial, 0, last, 0 );
      // total += (interv.t1==-1) ? interv.t3 : interv.t1;
      total += interv.t3;
      initial = last;
    } 
    last=jx;
    Interval interv = compute_interval( initial, 0, last, 0 );
    // total += (interv.t1==-1) ? interv.t3 : interv.t1;
    total += interv.t3;
    // std::cout << "Accumulated cost up to now " << total << std::endl; 
  } else if (jx<last && direction==0) {
    last = jx;
    if (!firsttime) {
      firsttime = true;
      Interval interv = compute_interval( initial+headvector, speed, last, 0 );
      // total += (interv.t1==-1) ? interv.t3 : interv.t1;
      total += interv.t3;
    } else {
      Interval interv = compute_interval( initial, 0, last, 0 );
      // total += (interv.t1==-1) ? interv.t3 : interv.t1;
      total += interv.t3;
    } 
    // std::cout << "Accumulated cost up to now " << total << std::endl; 
  } else if (jx<last && direction==1) {
    // accumulate cost
    direction=0;
    if (!firsttime) {
      firsttime = true;
      // compute the landing point
      int dist = (abs(speed)*(abs(speed)-1))/2;
      if (speed<0) landingpoint = (initial+headvector) - dist;
      else landingpoint = (initial+headvector) + dist;

      // the landing point is far away then the reversing point => 
      if (landingpoint >= last) {
        Interval interv = compute_interval( initial+headvector, speed, landingpoint, 0 );
        total += interv.t3;
        initial = landingpoint;
      } else {
        Interval interv = compute_interval( initial+headvector, speed, last, 0 );
        total += interv.t3;
        initial = last;
      }
    } else {
      Interval interv = compute_interval( initial, 0, last, 0 );
      // total += (interv.t1==-1) ? interv.t3 : interv.t1;
      total += interv.t3;
      initial = last;
    } 
    last=jx;    
    Interval interv = compute_interval( initial, 0, last, 0 );
    // total += (interv.t1==-1) ? interv.t3 : interv.t1;
    total += interv.t3;
  } else if (direction==-1) {
    last = jx;
    Interval interv = compute_interval( initial+headvector, speed, last, 0 );
    // total += (interv.t1==-1) ? interv.t3 : interv.t1;
    total += interv.t3;
  }
  // total += 1;
  // std::cout << "Total heuristic for this bivec " << total << std::endl;

  return total;
}


auto computeheur( BiVec bivec, const vector<UniVec>& univecs, const vector<int>& xs, const vector<int>& ys, int current_city) {

  int totx = computetime( bivec, univecs, xs, current_city, true );
  int toty = computetime( bivec, univecs, ys, current_city, false );

  if (totx > toty) return totx;
  else return toty;
}