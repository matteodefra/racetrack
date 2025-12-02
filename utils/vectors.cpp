#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

struct UniVec{
    int from;
    int to;

    friend std::ostream& operator <<(std::ostream& os, UniVec const& in) {
        return os << "(from:" << in.from << ", to:" << in.to << ")";
    }
};

struct BiVec{
    int ind_x;
    int ind_y;

    void print(vector<UniVec>& vecs){
        UniVec vx = vecs[ind_x];
        UniVec vy = vecs[ind_y];
        cout << "from:(" << vx.from << "," << vy.from << "),\t to:(" << vx.to << "," << vy.to << ")";
    }
};

inline bool operator==(const BiVec a, const BiVec b){
    return (a.ind_x == b.ind_x) && (a.ind_y == b.ind_y);
}

inline int distoperator( const vector<BiVec> a, const vector<BiVec> b ) {
    int counter = a.size();
    for (int i=0; i<a.size(); i++) {
        if (a[i]==b[i]) counter-=1;
    }
    return counter;
}

auto uni_vectors(int smax){
    vector<UniVec> uni_vecs;
    uni_vecs.push_back({0,0});
    for (int a = -smax; a <= smax; a++){
        for (int b = -smax; b <= smax; b++){
            if (abs(b - a) <= smax){
                if (a*b < 0 || (a != 0 && b == 0)){
                    uni_vecs.push_back({a,b});
                }
            }
        }
    }
    return uni_vecs;
}

auto uni_vectors_discrete(int smax){
     vector<UniVec> uni_vecs;
     uni_vecs.push_back({0,0});
     for (int a = -smax; a <= smax; a++){
       if (a != 0){
           uni_vecs.push_back({a,0});
       }
     }
     return uni_vecs;
}

double distance(int ax, int ay, int bx, int by){
    return sqrt((ax - bx)*(ax - bx) + (ay - by)*(ay - by));
}
  
bool is_between(int ax, int ay, int bx, int by, int cx, int cy){
    double indirect_distance = distance(ax, ay, cx, cy) + distance(bx,by,cx,cy);
    double direct_distance = distance(ax,ay,bx,by);
    return abs(indirect_distance - direct_distance) < 0.001;
}
  
vector<BiVec> bi_vectors(const vector<UniVec>& univecs){
    vector<BiVec> vecs;
    for (int ix = 0; ix < univecs.size(); ix++){
        for (int iy = 0; iy < univecs.size(); iy++){
            int ax = univecs[ix].from;
            int ay = univecs[iy].from;
            int bx = univecs[ix].to;
            int by = univecs[iy].to;
            if (is_between(ax,ay,bx,by,0,0)){
                vecs.push_back({ix, iy});
            }
        }
    }
    return vecs;
}

vector<int> boundingbox_limits( vector<int> &x, vector<int> &y ) {

    vector<int> bounds = {0, 0, 0, 0};
    bounds[0] = std::numeric_limits<int>::max();
    bounds[1] = std::numeric_limits<int>::min();
    bounds[2] = std::numeric_limits<int>::max();
    bounds[3] = std::numeric_limits<int>::min();
    int n =x.size();
    for (int i=0; i<n; i++) {
        if (x[i] <= bounds[0]) bounds[0] = x[i];
        if (x[i] >= bounds[1]) bounds[1] = x[i];


        if (y[i] <= bounds[2]) bounds[2] = y[i];
        if (y[i] >= bounds[3]) bounds[3] = y[i];
    }

    return bounds;
}

int compute_smax( vector<int> &x, vector<int> &y ) {
	int xi, xj, yi, yj;
	int n = x.size();
	double dist = 0.0;
	for (int i=0; i<n; i++) {
		for (int j=i+1; j<n; j++) {
			xi = x[i];
			xj = x[j];
			yi = y[i];
			yj = y[j];
			// cout << "distance between (" << xi << "," << yi << ") and (" << xj << "," << yj << ")" << endl;
			
			dist = max( dist, sqrt( pow(xi-xj, 2) + pow(yi-yj,2) )  );
		}
	}
	// std::cout << "distance " << dist << std::endl;
	return int(floor( sqrt(dist) ));
}


void test_vectors(){
    int smax = 10;
    auto uni_vecs = uni_vectors(smax);
    for (auto v : uni_vecs){
        cout << v.from << "," << v.to << endl;
    }
    cout << endl;

    auto vecs = bi_vectors(uni_vecs);
    for (auto vec : vecs){
        auto vx = uni_vecs[vec.ind_x];
        auto vy = uni_vecs[vec.ind_y];
        cout << "(" << vx.from << ", " << vy.from << "), (" << vx.to << ", " << vy.to << ")" << endl;
    }
    
    cout << uni_vecs.size() << " | " << vecs.size() << endl;
}