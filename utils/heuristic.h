int computetime( BiVec bivec, const vector<UniVec>& univecs, const vector<int>& cities, int current_city, int isx ) {
    int n = (int)cities.size();
    if (current_city < 0 || current_city >= n) return 0; // defensive
    if (current_city == n-1) return 0; // no next cities

    int initial_city = cities[current_city];
    int last = initial_city;
    int total = 0;
    bool started = false;

    int speed = isx ? (univecs[bivec.ind_x].to - univecs[bivec.ind_x].from)
                    : (univecs[bivec.ind_y].to - univecs[bivec.ind_y].from);
    int head = isx ? univecs[bivec.ind_x].to : univecs[bivec.ind_y].to;
    int dir = (speed >= 0) ? +1 : -1;

    for (int j = current_city + 1; j < n; ++j) {
        int jx = cities[j];
        if (jx == last) {
            // nothing changes
            continue;
        }
        bool moving_positive = (jx > last);

        if ((moving_positive && dir > 0) || (!moving_positive && dir < 0)) {
            // continuing in same direction
            last = jx;
            continue;
        } else {
            // reversal detected
            // handle reversal via helper that accepts current origin, speed, head, target=last
            if (!started) {
                // compute landingpoint from origin initial_city and head
                int origin = initial_city + head; // <-- *verify this is correct*
                int dist = (abs(speed)*(abs(speed)-1))/2; // verify formula
                int landing = origin + dir * dist;
                int cutoff = (dir > 0) ? std::min(landing, last) : std::max(landing, last);
                Interval interv = compute_interval(origin, speed, cutoff, 0);
                total += (interv.t1 == -1) ? interv.t3 : interv.t1;
                initial_city = cutoff;
                started = true;
            } else {
                Interval interv = compute_interval(initial_city, 0, last, 0);
                total += (interv.t1 == -1) ? interv.t3 : interv.t1;
                initial_city = last;
            }
            // flip direction
            dir = -dir;
            last = jx;
        }
    }

    // process last->final_jx if needed (but loop already handled all j)
    return total;
}


auto computeheur( BiVec bivec, const vector<UniVec>& univecs, const vector<int>& xs, const vector<int>& ys, int current_city) {

  int totx = computetime( bivec, univecs, xs, current_city, true );
  int toty = computetime( bivec, univecs, ys, current_city, false );

  if (totx > toty) return totx;
  else return toty;
}
