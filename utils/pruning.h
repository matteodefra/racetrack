#include <math.h>
#include <iostream>

using namespace std;;

struct Point {
  double x;
  double y;


  friend std::ostream& operator <<(std::ostream& os, Point const& in) {
    return os << "(" << in.x << "," << in.y << ")";
  }
};

Point midpoint(const Point& a, const Point& b) {
  Point ret;
  ret.x = (a.x + b.x) / 2;
  ret.y = (a.y + b.y) / 2;
  return ret;
}
  

/**
 * Filtering based on S total cost
 */
bool filteringS(const BiVec bivec, const vector<UniVec>& univecs, const vector<int>& xs, const vector<int>& ys, int city, int tottime){
  int n = xs.size();

  int x_from = univecs[bivec.ind_x].from;
  int x_to = univecs[bivec.ind_x].to;
  int y_from = univecs[bivec.ind_y].from;
  int y_to = univecs[bivec.ind_y].to;
  
  Interval one_x = compute_interval( xs[0], 0, xs[city]+x_to, x_to-x_from );
  Interval one_y = compute_interval( ys[0], 0, ys[city]+y_to, y_to-y_from );

  int t1 = min_value( one_x, one_y);


  Interval two_x = compute_interval( xs[city]+x_to, x_to-x_from, xs[n-1], 0 );
  Interval two_y = compute_interval( ys[city]+y_to, y_to-y_from, ys[n-1], 0 );

  int t2 = min_value( two_x, two_y );

  return (t1+t2)<=tottime;
}

bool separating_hyperplane( double m, double q, double xtest, double ytest, bool equality ) {
  if (equality) return (ytest <= m*xtest + q);
  else return (ytest < m*xtest + q);
}


bool check_mediatrice( Point& a, Point& b, Point& c, Point& d, Point& e ) {

  Point mid = midpoint( a, c );
  // cout << "point mid: " << mid << endl;
  double m, q;
  bool aside, dside, cside, eside;

  bool ad, ce;

  if (mid.y == b.y) {
    // same yaxis
    m = numeric_limits<double>::infinity();
    q = numeric_limits<double>::infinity();

    aside = true ? (a.x <= mid.x) : false;
    cside = true ? (c.x <= mid.x) : false;

    if (aside==true) {
      dside = (d.x <= mid.x);
      if (dside==true) {
        ad = true;
      } else {
        ad = false;
      }
    } else {
      dside = (d.x < mid.x);
      if (dside==false) {
        ad = true;
      } else {
        ad = false;
      }
    }

    if (cside==true) {
      eside = (e.x <= mid.x);
      if (eside==true) {
        ce = true;
      } else {
        ce = false;
      }
    } else {
      eside = (e.x < mid.x);
      if (eside==false) {
        ce = true;
      } else {
        ce = false;
      }
    }

  } else if (mid.x == b.x) {
    // same yaxis
    m = 0.0;
    q = mid.y;

    aside = true ? (a.y <= q) : false;
    cside = true ? (c.y <= q) : false;

    if (aside==true) {
      dside = (d.y <= q);
      if (dside==true) {
        ad = true;
      } else {
        ad = false;
      }
    } else {
      dside = (d.y < q);
      if (dside==false) {
        ad = true;
      } else {
        ad = false;
      }
    }

    if (cside==true) {
      eside = (e.y <= q);
      if (eside==true) {
        ce = true;
      } else {
        ce = false;
      }
    } else {
      eside = (e.y < q);
      if (eside==false) {
        ce = true;
      } else {
        ce = false;
      }
    }

  } else {

    m = (mid.y - b.y)/(mid.x - b.x);
    q = mid.y - m*mid.x;
    
    aside = separating_hyperplane( m, q, a.x, a.y, true );
    cside = separating_hyperplane( m, q, c.x, c.y, true );

    if (aside==true) {
      dside = separating_hyperplane( m, q, d.x, d.y, true );
      if (dside==true) {
        ad = true;
      } else {
        ad = false;
      }
    } else {
      dside = separating_hyperplane( m, q, d.x, d.y, false );;
      if (dside==false) {
        ad = true;
      } else {
        ad = false;
      }
    }

    if (cside==true) {
      eside = separating_hyperplane( m, q, e.x, e.y, true );
      if (eside==true) {
        ce = true;
      } else {
        ce = false;
      }
    } else {
      eside = separating_hyperplane( m, q, e.x, e.y, false );
      if (eside==false) {
        ce = true;
      } else {
        ce = false;
      }
    }

  }

  return (ad && ce);

}


bool check_tangent( Point& a, Point& b, Point& c, Point& d, Point& e ) {

  double m, q;
  bool aside, dside, cside, eside;
  int sign;
  int sign_d, sign_e;

  bool res = false;

  if (a.x == c.x) {
    m = numeric_limits<double>::infinity();
    q = b.y;

    if (q >= a.y) sign=1;
    else sign=-1;

    if (q>=d.y) sign_d=1;
    else sign_d=-1;

    if (q>=e.y) sign_e=1;
    else sign_e=-1;

    if (sign_d && sign_e) res=true;
    else if (sign_e && !sign_d) res=true;

  } else if (a.y == c.y) {
    m = 0.0;
    q = b.x;

    if (q >= a.x) sign=1;
    else sign=-1;

    if (q>=d.x) sign_d=1;
    else sign_d=-1;

    if (q>=e.x) sign_e=1;
    else sign_e=-1;

    if (sign_d && sign_e) res=true;
    else if (sign_e && !sign_d) res=true;
  
  } else {

    m = (a.y - c.y)/(a.x - c.x);
    q = b.y;

    if (m*a.x+q >= a.y) sign=1;
    else sign=-1;

    if (m*d.x+q >= d.y) sign_d=1;
    else sign_d=-1;

    if (m*e.x+q >= e.y) sign_e=1;
    else sign_e=-1;

    if (sign_d && sign_e) res=true;
    else if (sign_e && !sign_d) res=true;
  
  }

  return res;
  
}



bool checklanding( Point city, Point tail, Point head, vector<int> &boxlimits) {
  int xspeed = head.x - tail.x;
  int yspeed = head.y - tail.y;



  if (xspeed > 0) {
    int xdist = (xspeed*(xspeed-1))/2;

    if ( (tail.x-xdist) < boxlimits[0]) return false;
    else if ( (head.x+xdist) > boxlimits[1] ) return false;
  } else if (xspeed < 0) {
    int xdist = (abs(xspeed)*(abs(xspeed)-1))/2;

    if ( (head.x-xdist) < boxlimits[0]) return false;
    else if ( (tail.x+xdist) > boxlimits[1] ) return false;
  }

  if (yspeed > 0) {
    int ydist = (yspeed*(yspeed-1))/2;

    if ( (tail.y-ydist) < boxlimits[2]) return false;
    else if ( (head.y+ydist) > boxlimits[3] ) return false;
  } else if (yspeed < 0) {
    int ydist = (abs(yspeed)*(abs(yspeed)-1))/2;
  
    if ( (head.y-ydist) < boxlimits[2]) return false;
    else if ( (tail.y+ydist) > boxlimits[3] ) return false;
  }

  return true;
}