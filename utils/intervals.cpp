#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <list>
#include <vector>
#include <utility>


using namespace std;


struct Interval {
	int t1;
	int t2;
	int t3;

	bool operator==(const Interval& rhs)
	{
		return (t1 == rhs.t1) && (t2 == rhs.t2) && (t3 == rhs.t3);
	}

  friend std::ostream& operator <<(std::ostream& os, Interval const& in) {
    return os << "[" << in.t1 << "," << in.t2 << "] U [" << in.t3 << "," << "inf" << ")";
  }
};

struct Couple {
  int x;
  int s;
  int xp;
  int sp;
  int shift;
  friend std::ostream& operator <<(std::ostream& os, Couple const& in) {
    return os << "[" << in.x << "," << in.s << "," << in.xp << "," << in.sp << "]";
  }
};

struct Result {
  double t1;
  double t2;
};

// Solve quadratic equation ax^2 + bx + c = 0
Result solve_quadratic(double a, double b, double c) {
	// std::cout << "a: " << a << ",b: " << b << ",c: " << c << std::endl;
	double discriminant = b * b - 4 * a * c;
	if (discriminant < 0) {
		return {NAN, NAN};
	}
	double sqrt_d = sqrt(discriminant);
	double t1 = (-b - sqrt_d) / (2 * a);
	double t2 = (-b + sqrt_d) / (2 * a);
  	Result res = {min( t1, t2 ), max(t1, t2)};
	return res; 
}


float dmin_integer( int s, int sp, int t ) {
	float dmin = 0.0;
	float alpha = 0.0;

	alpha = ((s-sp)*(s+sp-1))/2;
	dmin = (2*sp-(t-s+sp)/2)*((t-s+sp)/2);
	
	return dmin + alpha;
}


float dmax_integer( int s, int sp, int t ) {
	float dmax = 0.0;
	float alpha = 0.0;

	alpha = ((s-sp)*(s+sp-1))/2;
	dmax = (2*s+(t-s+sp)/2)*((t-s+sp)/2);
	
	return dmax + alpha;
}


bool check_cond( int s, int sp, int t, int delta ) {
	float dmin = dmin_integer( s, sp, t );
	float dmax = dmax_integer( s, sp, t );

	// cout << "dmin integer: " << dmin << endl;
	// cout << "dmax integer: " << dmax << endl;

	return dmin <= delta && delta <= dmax;
}



// Solve dmin equation
Result dmin_eqn(int s, int sp, int delta) {
	double alpha, a, b, c;
	
	alpha = ((sp+s+1)*(sp-s))/2.0;
	a = -0.25;
	b = s + 0.5*(sp-s);
	c = -0.25 * pow(sp-s, 2) - s*(sp-s) + alpha - delta;
	
	Result t = solve_quadratic( a, b, c );
	
	return t;
}


// Solve dmax equation
Result dmax_eqn(int s, int sp, int delta) {
	double alpha, a, b, c;
	
	alpha = ((sp+s+1)*(sp-s))/2.0;
	a = 0.25;
	b = sp - 0.5*(sp-s);
	c = 0.25*pow(sp-s,2) - sp*(sp-s) + alpha - delta;
	
	Result t = solve_quadratic( a, b, c );
	
	return t;
}



void algo1(Couple points, Interval &racetrack) {

	int delta = points.xp - points.x;

	Result dmin = dmin_eqn( points.s, points.sp, delta );
	Result dmax = dmax_eqn( points.s, points.sp, delta );

	int tmin = abs(points.sp-points.s);
	int tstar = points.s + points.sp;

	int t1 = ceil(dmax.t2);

	float dminssp = dmin_integer( points.s, points.sp, tstar );

	if (delta>=dminssp) {
		racetrack.t3 = t1 + points.shift;
	} else {

		int t3 = ceil(dmin.t2);
		float dmintmin = dmin_integer( points.s, points.sp, tmin );

		if (delta<dmintmin) {
			racetrack.t3 = t3 + points.shift;
		} else {
			int t2 = floor(dmin.t1);
			if (t2>=t1) {
				racetrack.t1 = t1 + points.shift;
				racetrack.t2 = t2 + points.shift;
				racetrack.t3 = t3 + points.shift;
			} else {
				racetrack.t3 = t3 + points.shift;
			}
		}
	}

}

Couple lemma1( Couple points ) {
	// 	For all $(x,s)$ and $(x',s')$,  $\I((x,s),(x',s')) = \I((0,s),(x'-x,s'))$ 
	Couple res = { 0, points.s, points.xp-points.x, points.sp, points.shift };

	return res;
}

Couple lemma2( Couple points ) {
	// 	 For all $(x,s)$ and $(x',s')$,  $\I((x,s),(x',s')) = \I((-x,-s),(-x',-s'))$ 
	Couple res = { -points.x, -points.s, -points.xp, -points.sp, points.shift };

	return res;
}

//  wlog $x=0$ and $x'\ge 0$

// (s<0 & s'<=0) || (s<=0 & s'<0)
Couple case1( Couple points ) {

	int p1 = -(abs(points.s)*(abs(points.s)-1))/2;
	int p2 = points.xp + (abs(points.sp)*(abs(points.sp)+1))/2;
	
	Couple res = { 0, 0, p2-p1, 0, abs(points.s)+abs(points.sp) };

	return res;
}

// s<0 & s'>0
Couple case2( Couple points ) {

	int p1 = -(abs(points.s)*(abs(points.s)-1))/2;
	int p2 = points.xp - (points.sp*(points.sp+1))/2;

	Couple res;

	if (p1 <= p2) {
		res = { 0, 0, points.xp-p1, points.sp, abs(points.s) };
	} else {
		res = { 0, -points.s, -p2, 0, abs(points.sp) };
	}
	
	return res;
}

// s>0 & s'<0
Couple case3( Couple points ) {

	int p1 = points.x + (points.s*(points.s-1))/2;
	int p2 = points.xp + (abs(points.sp)*(abs(points.sp)+1))/2;

	Couple res;

	if (p1 >= p2) {
		res = { 0, 0, -points.xp+p1, -points.sp, abs(points.s) };
	} else {
		res = { 0, points.s, p2-points.x, 0, abs(points.sp) };
	}
	
	return res;
}

// s > 0 and s'> 0
Couple case4( Couple points ) {

	// find_interval
	return points;
}

// Applyin Lemma 1 and Lemma 2 to reduce x=0 & x'>0
Couple lemmas( Couple points ) {
	Couple res;
	if (points.xp >= 0 && points.xp >= points.x) {
		res = lemma1(points);
	} else if (points.xp >= 0 && points.xp < points.x) {
		res = lemma1(lemma2(points));
	} else if (points.xp < 0 && points.xp >= points.x) {
		res = lemma1(points);
	} else if (points.xp < 0 && points.xp < points.x) {
		res = lemma1(lemma2(points));
	}
	// cout << res << endl;

	return res;
}
// This works fine ./


Interval compute_interval(Couple points) {
	Couple converted = lemmas( points );
	Couple newcouple;

	if (converted.s==0 && converted.sp==0) return {-1, -1, (int)ceil(2*sqrt(converted.xp))};
	else if ((converted.s<=0 && converted.sp<0) || (converted.s<0 && converted.sp<=0)) {
		newcouple = case1( converted );
		return {-1, -1, (int)ceil(2*sqrt(newcouple.xp))+newcouple.shift};
	}
	else if (converted.s<0 && converted.sp>0) newcouple = case2(converted);
	else if (converted.s>0 && converted.sp<0) newcouple = case3(converted);
	else newcouple = converted;

	Interval racetrack = {-1,-1,-1}; 
	algo1(newcouple, racetrack);

	return racetrack;

}


Interval compute_interval(int x, int sx, int xp, int sp) {
	if (x==xp && sx==sp && sx==0) return {0,0,0};
	else return compute_interval({x, sx, xp, sp});
}
