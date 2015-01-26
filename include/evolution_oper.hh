#ifndef EVOLUTION_OPER_HH
#define EVOLUTION_OPER_HH

/*
metl: A generic framework for sequential and parallel metaheuristics
Copyright (c) 2005-2015, Sylvain Ouellet


Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

*/


// This file define generic evolution operators for selection,
// replacement and crossover.



#include <vector>
#include <algorithm>
#include <utility>

#include "exchange_oper.hh"

#include "abstract_selection.hh"
#include "abstract_replace.hh"
#include "abstract_crossover.hh"

namespace metl {


////////////////////////////////////////////////
// generic selection operations
////////////////////////////////////////////////


template <class prob_t>
struct select_random: public abstract_selection<prob_t> {
  // select two random individuals in the population
  // p1!=p2
  void operator()(const typename std::vector<typename prob_t::soleval_t>& pop, 
		  typename std::vector<typename prob_t::soleval_t>::const_iterator& p1, 
		  typename std::vector<typename prob_t::soleval_t>::const_iterator& p2) const
  {
    unsigned n1 = rng(pop.size());
    unsigned n2 = rng(pop.size()-1);

    if (n2>=n1) ++n2;

    p1 = pop.begin()+n1;
    p2 = pop.begin()+n2;
  }
};



template <class prob_t>
struct select_rank: public abstract_selection<prob_t> {
  // select two random individuals in the population
  // p1!=p2
  // use rank selection
  void operator()(const typename std::vector<typename prob_t::soleval_t>& pop, 
		  typename std::vector<typename prob_t::soleval_t>::const_iterator& p1, 
		  typename std::vector<typename prob_t::soleval_t>::const_iterator& p2) const
  {
    // FIXME: test
    unsigned n1 = ranking(pop.size());
    unsigned n2 = ranking(pop.size());
    while(n1==n2) n2=ranking(pop.size());   // FIXME: inefficient

    p1 = pop.begin()+n1;
    p2 = pop.begin()+n2;
  }
private:
  inline unsigned ranking(unsigned size) const {
    return unsigned((2*size-sqrt(double(4*size*size - 8*rng(size*size/2))))/2);
  }
};


template <class prob_t>
struct select_tournament: public abstract_selection<prob_t> {

  void operator()(const typename std::vector<typename prob_t::soleval_t>& pop, 
		  typename std::vector<typename prob_t::soleval_t>::const_iterator& p1, 
		  typename std::vector<typename prob_t::soleval_t>::const_iterator& p2) const
  {
    // FIXME: test
    unsigned n1 = rng(pop.size());
    unsigned n2 = rng(pop.size());

    unsigned p_1 = std::min(n1,n2);   // they are sorted, so keep lower rank

    unsigned n3 = rng(pop.size()-1);
    unsigned n4 = rng(pop.size()-1);
    if (n3>=p_1) ++n3;
    if (n4>=p_1) ++n4;

    unsigned p_2 = std::min(n3,n4);

    p1 = pop.begin()+p_1;
    p2 = pop.begin()+p_2;
    
  }
};


template <class prob_t>
struct select_roulette: public abstract_selection<prob_t> {
  void operator()(const typename std::vector<typename prob_t::soleval_t>& pop, 
		  typename std::vector<typename prob_t::soleval_t>::const_iterator& p1, 
		  typename std::vector<typename prob_t::soleval_t>::const_iterator& p2) const
  {
    typename prob_t::eval_t sum=0;
    std::vector<typename prob_t::eval_t> s(pop.size());
    
    for (typename std::vector<typename prob_t::soleval_t>::const_iterator i=pop.begin(); i!=pop.end(); ++i)
      {
	sum+=i->second;
	s[i-pop.begin()] = (pop.back().second - i->second);

	if (i!=pop.begin) {
	  s[i-pop.begin()] += s[i-1-pop.begin()];
	}
      }

    // choose a random number between 0 and sum-pop.back().second
    typename prob_t::eval_t choice = rng(sum-pop.back().second);
    
    p1 = std::lower_bound(s.begin(), s.end(), choice);
    p2 = std::lower_bound(s.begin(), s.end(), choice);
    
  }
};




///////////////////////////////////////////////////////
// generic population replacement operations
///////////////////////////////////////////////////////


// I'm not sure but I think that this version might apply too much selection
// presure.
template <class prob_t>
struct replace_worst: public abstract_replace<prob_t> {

  // this version keeps the best pop.size() individuals from (pop union childs)
  void operator()(std::vector<typename prob_t::soleval_t>& pop, 
		  std::vector<typename prob_t::soleval_t>& childs,
		  const typename std::vector<typename prob_t::soleval_t>::iterator& wp) const 

  {
    // sort childrens
    sort(childs.begin(), childs.end(), individus_compare<typename prob_t::soleval_t>);

    typename std::vector<typename prob_t::soleval_t>::const_iterator ci = childs.begin();
    typename std::vector<typename prob_t::soleval_t>::iterator i = std::lower_bound(pop.begin(), pop.end(), *ci, individus_compare<typename prob_t::soleval_t>);

    if (i==pop.end()) return;

    while(1) {
      pop.insert(i,*ci);
      pop.pop_back();
      if (++ci==childs.end()) return;
      while (individus_compare<typename prob_t::soleval_t>(*i,*ci)) 
	{
	  if (++i==pop.end()) return;
	}
    }
  }
};




template <class prob_t>
struct replace_worst_from_pop: public abstract_replace<prob_t> {
  // this version replaces the worst individuals from pop with the childs
  void operator()(std::vector<typename prob_t::soleval_t>& pop, 
		  std::vector<typename prob_t::soleval_t>& childs,
		  const typename std::vector<typename prob_t::soleval_t>::iterator& wp) const 
  {
    const std::vector<typename prob_t::soleval_t>& const_childs = childs;

    typename std::vector<typename prob_t::soleval_t>::const_iterator ci = const_childs.begin();

    for (typename std::vector<typename prob_t::soleval_t>::iterator i=pop.end()-const_childs.size(); 
	 i!=pop.end(); ++i, ++ci)
      {
	*i = *ci;
      }
  }
};


template <class prob_t>
struct replace_worst_parent: public abstract_replace<prob_t> {
  void operator()(std::vector<typename prob_t::soleval_t>& pop, 
		  std::vector<typename prob_t::soleval_t>& childs,
		  const typename std::vector<typename prob_t::soleval_t>::iterator& wp) const 
  {
    const std::vector<typename prob_t::soleval_t>& const_childs = childs;
    // keeps the best child and replace the worst parent with it.
    *wp = *std::min_element(const_childs.begin(), const_childs.end(),individus_compare<typename prob_t::soleval_t>);
  }
};




/////////////////////////////////////////////////////
// Generic blind crossover operations
////////////////////////////////////////////////////


// no crossover, only copy parent 1 into child
// looks pointless, but it is usefull for testing
template <class prob_t>
struct no_xover: public abstract_crossover<prob_t> {

  void operator()(const typename prob_t::soleval_t& parent1, 
		  const typename prob_t::soleval_t& parent2, 
		  typename prob_t::soleval_t& child) const 
  {
    child=parent1;
  }
};


  // perform uniform crossover
template <class prob_t>
struct uniform_xover: public abstract_crossover<prob_t> {

  void operator()(const typename prob_t::soleval_t& parent1, 
		  const typename prob_t::soleval_t& parent2, 
		  typename prob_t::soleval_t& child) const 
  {
    child=parent1; // FIXME: a complete copy should be avoided because some of the values will be overwritten anyway
    typename prob_t::sol_t::iterator s_it = child.first.begin();

    for (typename prob_t::sol_t::const_iterator it = parent2.first.begin();
	 it!=parent2.first.end();
	 ++it, ++s_it)
      {
	if (rng()<0.5)
	  *s_it= *it;
      }
    child.second = prob_t::instance().evaluation(child.first);
  }
};



  // perform single point crossover
template <class prob_t>
struct single_point_xover: public abstract_crossover<prob_t> {

  void operator()(const typename prob_t::soleval_t& parent1, 
		  const typename prob_t::soleval_t& parent2, 
		  typename prob_t::sol_t& child) const 
  {
    child=parent1; // FIXME: a complete copy should be avoided because some of the values will be overwritten anyway

    const unsigned crosspoint = rng(parent2.first.size());
    typename prob_t::sol_t::iterator s_it = child.first.begin()+crosspoint;
    
    for (typename prob_t::sol_t::const_iterator it = parent2.first.begin()+crosspoint;
	 it!=parent2.first.end();
	 ++it, ++s_it)
      {
	  *s_it= *it;
      }
    child.second = prob_t::instance().evaluation(child.first);
  }
};




  // perform 2 point crossover
template <class prob_t>
struct two_point_xover: public abstract_crossover<prob_t> {

  void operator()(const typename prob_t::soleval_t& parent1, 
		  const typename prob_t::soleval_t& parent2, 
		  typename prob_t::sol_t& child) const 
  {
    child=parent1; // FIXME: a complete copy should be avoided because some of the values will be overwritten anyway

    const unsigned point1 = rng(parent2.first.size());
    const unsigned point2 = rng(parent2.first.size());

    const typename prob_t::sol_t::const_iterator stop_it = 
      parent2.first.begin()+point2;
    
    typename prob_t::sol_t::const_iterator it = parent2.first.begin()+point1;
    typename prob_t::sol_t::iterator s_it = child.first.begin()+point1;
    while (it!=stop_it) {
      
      *s_it = *it;
      
      ++it;
      ++s_it;
      if (it==parent2.first.end()) {
	it = parent2.first.begin();
	s_it = child.first.begin();
      }
    }

    child.second = prob_t::instance().evaluation(child.first);
  }
};



}
#endif
