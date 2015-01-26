#ifndef CANDIDATE_LISTS_HH
#define CANDIDATE_LISTS_HH
class tsp_prob;

#include <vector>
struct candidate_lists {
  candidate_lists(const tsp_prob& my_tsp)
    : _empty(true), tsp(my_tsp) , c()
  {};
  void init(unsigned k);

  const std::vector<unsigned>& operator[](unsigned v) const { return c[v]; }
  bool empty() const { return _empty; }

private:
  bool _empty;
  const tsp_prob& tsp;
  std::vector<std::vector<unsigned> > c;
};

#endif
