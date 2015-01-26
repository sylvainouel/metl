#include "candidate_lists.hh"
#include "tsp_prob.hh"
#include <iostream>

using namespace std;

struct cmp_f {
  cmp_f(const vector<int>& _d):d(_d) {}
  bool operator()(unsigned v1, unsigned v2) { return d[v1]<d[v2]; }
private:
  const vector<int>& d;
};

void candidate_lists::init(unsigned k)
{
  cout << "building candidate lists .." << flush;
  _empty=false;
  c.reserve(tsp.size());
  for (unsigned i=0; i<tsp.size(); i++) {
    c.push_back(vector<unsigned>());
    vector<unsigned>& x =c.back();
    x.reserve(k+1);
    vector<int> d(tsp.size());   // this is a vector of distance from that city
    cmp_f functio(d);  // make a comparison fucntion;

    for (unsigned j=0; j<tsp.size(); j++) {
      if (i==j) continue;
      d[j] = tsp.dist(i,j);
      if (x.size()==k && d[j]>d[x[0]]) continue;

      x.push_back(j);
      
      push_heap(x.begin(), x.end(), functio);
      if (x.size()>k) {
	pop_heap(x.begin(), x.end(), functio);
	x.pop_back();
      }
    }
    sort_heap(x.begin(), x.end(), functio);
  }
  cout << "done" << endl;
}


