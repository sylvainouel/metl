#ifndef TSP_EAX_HH
#define TSP_EAX_HH

#include "abstract_crossover.hh"
#include "tsp_prob.hh"

#include <vector>
#include <deque>

class tour;

struct tsp_eax: public metl::abstract_crossover<tsp_prob> {
  void operator()(const tsp_prob::soleval_t& Ap,
		  const tsp_prob::soleval_t& Bp,
		  tsp_prob::soleval_t& Cp) const;

private:
  void gready_subtours_recombine(std::vector<std::deque<int> >& subtours) const;

};

#endif
