#ifndef TABU_HH
#define TABU_HH

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



// generic tabu engine
// variants are: tabu, tabu_gain, tabu_ns_omp, tabu_ns_mpi


#include <limits>
#include <iostream>
#include <algorithm>

#include "neighborhood_oper.hh"
#include "move_reduction.hh"
#include "meta_gen.hh"
#include "dummy_oper.hh"
#include <time.h>

namespace metl {

// this is an internal class used by tabu algorithms
template<class prob_t, class _move, class _tabu_list>
class tabu_kernel : public move_reduction<prob_t, _move> {
  typedef move_reduction<prob_t, _move> base;
public:

  tabu_kernel(const typename prob_t::sol_t& sol, 
	      _move& best_move,
	      const typename prob_t::eval_t& current_eval, 
	      const typename prob_t::eval_t& best_eval, 
	      const _tabu_list& tabu_list,
	      const unsigned& cycle,
	      bool using_mpi=false)
    : move_reduction<prob_t, _move>(best_move, using_mpi),
      
      s(sol), 
      ce(current_eval), 
      be(best_eval), 
      _tl(tabu_list),
      _cycle(cycle)
  {}


  // tabu kernel. Called from the neighborhood evaluation (defined by user)
  // compute move evaluation
  inline bool operator()(const _move m) {
    operator()(m,m.internal_cost(s));
    return false;  // always return false because move not directly applied
  }

  // evaluation is available. This is called either from previous version or
  //    from tabu_gain because the gain structure already has the cost of the
  //    move pre-computed
  inline bool operator()(const _move m, const typename prob_t::eval_t& e) {
    if ((e<base::thread_cost()) && (!_tl.is_tabu(m, s, _cycle) || (e+ce < be))) {
      base::update_thread_move_and_cost(m,e);
    }
    return false;  // always return false because move not directly applied
  }
  
private:
  const typename prob_t::sol_t& s;
  const typename prob_t::eval_t& ce;
  const typename prob_t::eval_t& be;

  const _tabu_list& _tl;
  const unsigned& _cycle;
};




// this is a generic implementation of the tabu algorithm
// this is not usable directly. Use one of its derived class.
template <class prob_t, 
	  class _move, 
	  class _neighborhood, 
	  class _tabu_list, 
	  class generator_type=no_generator<prob_t> >
struct tabu_base: public meta_gen<prob_t, generator_type> {

  void set_tenur_out(unsigned tabu_tenur) {
    // tabu tenur
    tenur_out=tabu_tenur;
  }

  void set_tenur_in(unsigned tabu_tenur) {
    // tabu tenur
    tenur_in=tabu_tenur;
  }

  void set_tenur(unsigned tabu_tenur) {
    // tabu tenur
    tenur_out=tenur_out=tabu_tenur;
  }

  void set_n_iter(unsigned n_iter) {
    // number of iteration of tabu search
    _n_iter = n_iter;
  }

  void set_return_current(float _p) {
    // probability of returning with current solution rather than the
    // best solution if the algorithm could not improve the best
    // solution
    return_current = _p;
  }

#ifdef XLC_WORKAROUND
  // IBM XLC does not let the derived classes to access the protected member.
  // this is stupid and I could not find out why.
public:
#else
protected:
#endif

#ifndef XLC_WORKAROUND
  using meta_gen<prob_t, generator_type>::operator();
#endif

  using meta_gen<prob_t, generator_type>::generator;
  
  tabu_base(unsigned tabu_tenur=8, unsigned n_iter=200, float r_current=0)
    : tenur_in(tabu_tenur),
      tenur_out(tabu_tenur),
      _n_iter(n_iter),
      return_current(r_current)
  {}

  virtual ~tabu_base() {}

  template <class _nh_eval, class _gain_t, class PE, class RNG>
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in, 
					_nh_eval& nh_eval, _gain_t& gain, 
					PE& periodic_exchange, RNG& srng, 
					bool mpi=false) const {

    typename prob_t::soleval_t se(se_in);

    CHECK_EVAL(se);

    _move next_move;  // best non tabu move or aspirated
    _tabu_list tl;

    typename prob_t::soleval_t best=se;
    unsigned cycle=1; // need to start at 1 because tabu list is filled with 0s. (every moves are tabu if cycle=0)

    typedef tabu_kernel<prob_t, _move, _tabu_list> tabu_k_t;

    tabu_k_t tabu_k(se.first, next_move, se.second, best.second, tl, cycle, mpi);

    const typename prob_t::eval_t optimum = prob_t::instance().optimum();

    while(cycle<=_n_iter) {
      nh_eval(tabu_k, se.first);
      tabu_k.reduce();

      if (tabu_k.cost() == std::numeric_limits<typename prob_t::eval_t>::max())
	{
	  std::cerr << "No move left. (tabu tenur probably too high)"<< std::endl;
	  abort();
	}

      gain.update_before(next_move, se.first);
      //********* APPLY SELECTED MOVE
      next_move(se.first);           // apply the move to sol
      //*********
      // one of the two update functions usually does nothing and will
      // be optimised away.
      gain.update_after(next_move, se.first);
      
      // make the move tabu
      tl.make_tabu(next_move, se.first, cycle,  
		   srng(9*tenur_in/10, 11*tenur_in/10), 
		   srng(9*tenur_out/10, 11*tenur_out/10));
      
      se.second += tabu_k.cost();

      if (se.second < best.second) {
	// save new best solution
	best = se;
#ifdef TABU_TRACE
	std::cout <<
#ifdef USE_MPI
	  MPI::COMM_WORLD.Get_rank() << " " <<
#else 
	  "0 "<<
#endif
	  cycle << " xxxx "<< se.second << std::endl;
#endif
	if (se.second <= optimum) break;
      }

      if (periodic_exchange(se)) {
	gain.init(se.first);
      }

      tabu_k.reset();
      ++cycle;
    }

    if ((se_in.second > best.second) || (srng()>return_current)) {
      se = best;
    }

    CHECK_EVAL(se);

    return se;
  }


protected:
  unsigned tenur_in, tenur_out;
  unsigned _n_iter;
private:
  float return_current;
};


// this is a basic version
template<class prob_t, 
	 class _move, 
	 class _neighborhood, 
	 class _tabu_list, 
	 class generator_type=no_generator<prob_t> >
class tabu : public tabu_base<prob_t, _move, _neighborhood, _tabu_list, generator_type> {

  typedef tabu_base<prob_t, _move, _neighborhood, _tabu_list, generator_type> base;

public:
    using base::operator();
    using base::generator;


  tabu(unsigned tabu_tenur=8, unsigned n_iter=1000)
    : base(tabu_tenur, n_iter)
  {}
  const std::string name() const { return "Tabu Search, base version"; }

  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se) {
    dummy_op dummy;
    return operator()(se, dummy);
  }

protected:
  template <class PE> 
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se, PE& periodic_exchange) {
    _neighborhood n;
    basic_nh_eval<_neighborhood, _move, typename prob_t::sol_t> nh_eval(n);
    dummy_op dummy;

    return base::operator()(se, nh_eval, dummy, periodic_exchange, rng);
  }
};


// this version uses OpenMP for parallel neighborhood evaluation
template<class prob_t, 
	 class _move, 
	 class _neighborhood, 
	 class _tabu_list, 
	 class generator_type=no_generator<prob_t> >
class tabu_ns_omp : public tabu_base<prob_t, _move,_neighborhood, _tabu_list, generator_type> {

  typedef tabu_base<prob_t, _move, _neighborhood, _tabu_list, generator_type> base;

public:
   using base::operator();
   using base::generator;


  tabu_ns_omp(unsigned tabu_tenur=8, unsigned n_iter=1000)
    : base(tabu_tenur, n_iter)
  {}

  const std::string name() const { return "Parallel tabu search based on neighborhood separation using OpenMP"; }

  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se) {
    dummy_op dummy;
    return operator()(se, dummy);
  }

protected:
  template<class PE>
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se, PE& periodic_exchange) {
    _neighborhood n;
    ns_nh_omp<_neighborhood, _move, typename prob_t::sol_t> nh_eval(n);
    dummy_op dummy;
    sync_rng srng(time(0));

    return base::operator()(se, nh_eval, dummy, periodic_exchange, srng);
  }
};




#ifdef USE_MPI

// this version use MPI to seperate neighborhood evaluation
template<class prob_t, 
	 class _move, 
	 class _neighborhood, 
	 class _tabu_list, 
	 class generator_type=no_generator<prob_t> >
class tabu_ns_mpi : public tabu_base<prob_t, _move, _neighborhood, _tabu_list, generator_type> {

  typedef tabu_base<prob_t, _move, _neighborhood, _tabu_list, generator_type> base;

public:
   using base::operator();
   using base::generator;


  tabu_ns_mpi(unsigned tabu_tenur=8, unsigned n_iter=1000)
    : base(tabu_tenur, n_iter)
  {}
  
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se) {
    dummy_op dummy;
    return operator()(se, dummy);
  }

  const std::string name() const { return "Parallel tabu search based on neighborhood separation using MPI"; }

protected:
  template<class PE>
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in, PE& periodic_exchange) {
    typename prob_t::soleval_t se(se_in);
    sol_reducter<prob_t>::Allreduce(se);

    _neighborhood n;
    ns_nh_mpi<_neighborhood, _move, typename prob_t::sol_t> nh_eval(n);
    dummy_op dummy;
    sync_rng srng(time(0));

    return base::operator()(se, nh_eval, dummy, periodic_exchange, srng,true);
  }
};
#endif




// this version uses a gain structure to accelerate the computation of the next move
template<class prob_t,
	 class _move, 
	 class _gain_struct_t, 
	 class _tabu_list, 
	 class generator_type=no_generator<prob_t> >
class tabu_gain : public tabu_base<prob_t, _move, dummy_neighborhood<prob_t>, _tabu_list,generator_type> {
  typedef tabu_base<prob_t, _move, dummy_neighborhood<prob_t>, _tabu_list, generator_type> base;

public:
   using base::operator();
   using base::generator;


  tabu_gain(unsigned tabu_tenur=8, unsigned n_iter=1000)
    : base(tabu_tenur, n_iter)
  {}

  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se) {
    dummy_op dummy;
    return operator()(se, dummy);
  }
  
  const std::string name() const { return "Tabu Search using a gain structure"; }  

#ifdef XLC_WORKAROUND
  // IBM XLC does not let the derived classes to access the protected member.
  // this is stupid and I could not find out why.
public:
#else
protected:
#endif

  template <class PE>
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se, PE& periodic_exchange) {

    _gain_struct_t _gain;
    _gain.init(se.first);
    gain_nh_eval<_gain_struct_t, _move, typename prob_t::sol_t> nh_eval(_gain);

    return base::operator()(se, nh_eval, _gain, periodic_exchange, rng);
  }
};


}

#endif
