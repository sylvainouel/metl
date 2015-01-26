#ifndef SA_HH
#define SA_HH

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

// generic simulated annealing


#include "sprng_rand.hh"
#include <math.h>

#include "exchange_oper.hh"
#include "meta_gen.hh"
#include "dummy_oper.hh"
#include "sa_oper.hh"

namespace metl {

// generic Simulated Annealing kernel
template<class prob_t, 
	 class _move, 
	 class _neighborhood, 
	 class generator_type=no_generator<prob_t>, 
	 class _accept_scheme=metropolis<prob_t,_move>, 
	 class _cooling_scheme=cooling_geometric_steps<_neighborhood> >
struct simulated_annealing: public meta_gen<prob_t, generator_type> {

  using meta_gen<prob_t, generator_type>::operator();
  using meta_gen<prob_t, generator_type>::generator;

  // functions that set algorithm parameters
  void set_init_temp(double init_temp) {
    init_T=init_temp;
  }

  void set_final_temp(double final_temp) {
    final_T=final_temp;
  }

  void set_max_rejects(unsigned m_rejects) {
    max_rejects=m_rejects;
  }

  // this is to allow configuring the cooling scheme in a generic maner.
  _cooling_scheme& cooling_scheme() { return cooler; }
  
  _accept_scheme& accept_scheme() { return accept; }

 
  // constructor
  simulated_annealing(double init_temp=10, double final_temp=0.1) 
    :init_T(init_temp),
     current_temp(0),
     final_T(final_temp),
     cooler(&current_temp),
     accept(current_temp),
     max_rejects(std::numeric_limits<unsigned>::max())
  {  };


  simulated_annealing(const simulated_annealing& r) 
    : init_T(r.init_T),
      current_temp(r.current_temp),
      final_T(r.final_T),
      cooler(&current_temp),
      accept(current_temp),
      max_rejects(r.max_rejects)
  {  
    cooler = r.cooler;
    accept = r.accept;
  }


  virtual ~simulated_annealing() {}

  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se) {
    dummy_op dummy;
    return operator()(se,dummy);
  }


  const std::string name() const { return "Simulated annealing"; }
  


#ifdef XLC_WORKAROUND
  // IBM XLC does not let the derived classes to access the protected member.
  // this is stupid and I could not find out why.
public:
#else
protected:
#endif

  // do the search
  template <class PE> 
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in, PE& periodic_exchange) {
    typename prob_t::soleval_t se(se_in);
    CHECK_EVAL(se);
    
    accept.set_sol_and_eval(&se.first, &se.second);

    current_temp = init_T;

    _neighborhood n;
    cooler.set_neighborhood(n);
    const typename prob_t::eval_t optimum = prob_t::instance().optimum();
    accept.reset_counters();

#ifdef SA_TRACE
    typename prob_t::eval_t best=se.second;
    unsigned it=0;
#endif 

    while (current_temp>final_T && 
	   se.second > optimum && 
	   accept.get_rejects()<max_rejects) {

      n(accept, se.first);
      cooler();

      // exchange solution. This is for the Multiple Markov chains 
      //   parallel algorithms
      periodic_exchange(se);
#ifdef SA_TRACE
      if (se.second<best) {
	best = se.second;
	std::cout << 
#ifdef USE_MPI
	  MPI::COMM_WORLD.Get_rank() << " " <<
#else 
	  "0 "<<
#endif
	  it << " "<< current_temp <<" "<<se.second<<std::endl;
      }
      ++it;
#endif

    }
    CHECK_EVAL(se);
    return se;
  }


private:
  double init_T;
  double current_temp;
  double final_T;
  _cooling_scheme cooler;
  _accept_scheme accept;
  
  unsigned max_rejects;
};

}

#endif
