#ifndef TOUR_HH
#define TOUR_HH

#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>

class tsp_prob;


struct tour {
  //  tour(const tsp_prob& tsp); 
  tour();
  tour(const std::vector<unsigned>& x);

  inline unsigned next(unsigned a) const { 
    unsigned i=B[a]+1;
    if (i==size()) i=0;
    return A[i]; 
  }

  inline unsigned prev(unsigned a) const { 
    unsigned i=B[a];
    if (i==0) i=size();
    --i;
    return A[i]; 
  }
  inline bool between(unsigned a, unsigned b, unsigned c) const {
    if ((B[a] < B[b] && B[b] < B[c]) ||
	(B[a] < B[b] && B[c] < B[a]) ||
	(B[b] < B[c] && B[c] < B[a])) return true;
    return false;
  }

  void flip(unsigned a, unsigned b, unsigned c, unsigned d);
  
  unsigned operator[](unsigned i) const { return A[i]; }

  tour& operator=(const std::vector<unsigned>& _A) {
    A = _A;
    updateB();
    return *this;
  }

  bool operator==(const tour& rhs) const {
    if (A==rhs.A) return true;
    return false;
  }
  const std::vector<unsigned>& get_tour() const { return A; }
  unsigned size() const { return A.size(); }
  
  bool checkB() const;
  
  void print();
  void double_bridge();

  template<class Archive>
  void serialize(Archive & ar, const unsigned int /* file_version */){
    ar & A & B;
  }


private:
  std::vector<unsigned> A;    // hold the cities
  std::vector<unsigned> B;    // hold the positions of the cities in vector A

  void updateB();
};

#endif
