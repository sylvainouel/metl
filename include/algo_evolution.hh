#ifndef ALGO_EVOLUTION_HH
#define ALGO_EVOLUTION_HH

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

#include <assert.h>
#include <utility>

#include "dummy_oper.hh"
#include "evolution_oper.hh"
#include "generator.hh"

#include <iostream>


namespace metl {

template <class prob_t, 
	  class generator_type, 
	  class crossover_op, 
	  class mutation_op=no_mutation, 
	  class localsearch=no_localsearch<prob_t>, 
	  class select_op=select_random<prob_t>, 
	  class reduce_op=replace_worst_from_pop<prob_t> >
class evolution: public meta_gen<prob_t, generator_type> {

public:
  // this allow user to configure internal objects.
  localsearch& local_search() { return ls; } 
  select_op& selection() { return _selection; }
  crossover_op& crossover() { return x_over; }   
  reduce_op& reduction() { return reduce_population; }


  void set_popsize(unsigned popsize) { _popsize = popsize; }
  void set_generations(unsigned generations) { _generations = generations; }
  void set_mutation_rate(float mutation_rate) { mut_rate = mutation_rate; }
  void set_childs_per_gen(unsigned childs_per_gen) {_childs = childs_per_gen; }


  evolution(unsigned popsize=10, unsigned generations=200, float mutation_rate=0.0, unsigned childs=1)
    : _popsize(popsize),
      _generations(generations),
      mut_rate(mutation_rate),
      _childs(childs),
      ls(),
      _selection(),
      x_over(),
      reduce_population()
  {
    assert(_popsize>1);
    assert(_childs<=_popsize);
  }


  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se) {
    dummy_op dummy;
    return operator()(se,dummy);
  }


  const std::string name() const { return "Evolutionary algorithm with localsearch=("+this->ls.name()+")"; }

#ifdef XLC_WORKAROUND
  // IBM XLC does not let the derived classes to access the protected member.
  // this is stupid and I could not find out why.
public:
#else
protected:
  using meta_gen<prob_t, generator_type>::operator();
#endif


  template <class PE>
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in, PE& periodic_exchange) {
    typename prob_t::soleval_t se(se_in);
    CHECK_EVAL(se);

    std::vector<typename prob_t::soleval_t> pop(_popsize);

    pop[0]=se;
//     std::cout << pop[0].second << " " << std::flush;
    
    int signed_pops = static_cast<int>(_popsize);
    int i;

    for (i=1; i<signed_pops; ++i) {
      pop[i]=this->gen();
      CHECK_EVAL(pop[i]);
    }

    std::vector<typename prob_t::soleval_t> childrens(_childs);

    mutation_op mutation; 
    mutation.set_rate(mut_rate);

    const int signed_nchilds = static_cast<int>(_childs);

    unsigned wp_i=0;
    int n;
    unsigned gen;
    const prob_t& instance = prob_t::instance();
    
    // sort initial pop
    std::sort(pop.begin(), pop.end(), individus_compare<typename prob_t::soleval_t>);
    
    for (gen=0; 
	 gen<_generations && pop.front().second > instance.optimum();    // FIXME: méthode de terminaison parallèle.
	 ++gen) {
      

      for (n=0; n<signed_nchilds; ++n) {
	typename std::vector<typename prob_t::soleval_t>::const_iterator p1, p2;
	
	_selection(pop, p1, p2);
	x_over(*p1, *p2, childrens[n]);
	
	if (mutation.rate()>0 && rng()<mutation.rate()) {
	  mutation(childrens[n].first);
	  // re-evaluate after mutation
	  childrens[n].second = instance.evaluation(childrens[n].first);
	}
	
	childrens[n]=ls(childrens[n]);
	CHECK_EVAL(childrens[n]);
	
	unsigned lwp = std::max(p1-pop.begin(),p2-pop.begin());
	if (lwp > wp_i) wp_i=lwp;
	
	reduce_population(pop, childrens, pop.begin()+wp_i);      // does not keep population sorted
	periodic_exchange(*(pop.begin()+rng(pop.size())));
	//	  periodic_exchange.recv(pop.front());
	std::sort(pop.begin(), pop.end(), individus_compare<typename prob_t::soleval_t>);
	  
#ifdef USE_MPI
	if (MPI::COMM_WORLD.Get_rank()==0)
#endif
 	  if(gen%32==0) 
	    {
	      std::cout << "\tgeneration:" << gen << "  best eval: " << pop.front().second << std::endl;
	    }
      }
    }

    const typename prob_t::soleval_t& x = pop.front();
    CHECK_EVAL(x);
    return x;
  }


private:
  unsigned _popsize;     // population size
  unsigned _generations;  // number of generations
  float mut_rate;   // mutation rate
  unsigned _childs;  // this is the number of childs per generation

  localsearch ls;
  select_op _selection;
  crossover_op x_over;  
  reduce_op reduce_population;

};



}

#endif
