#ifndef DESCENT_HH
#define DESCENT_HH

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

#include "metl_def.hh"

#ifdef USE_MPI
#ifdef HAVE_MPIPP_H     
#include <mpi++.h>
#else
#include <mpi.h>
#endif  // HAVE_MPIPP_H
#endif  // USE_MPI


#include "move_reduction.hh"

#include "neighborhood_oper.hh"
#include "meta_gen.hh"
#include "dummy_oper.hh"

namespace metl {

template<class prob_t, class _move>
struct keep_best : public move_reduction<prob_t, _move> {
  typedef move_reduction<prob_t, _move> base;

  keep_best(const typename prob_t::sol_t& sol, 
	    _move& best_move, 
	    bool using_mpi=false)
    : move_reduction<prob_t, _move>(best_move, using_mpi),
      s(sol)
  {  }

  inline bool operator()(const _move m) {
    const typename prob_t::eval_t e=m.internal_cost(s);
    operator()(m,e);
    return false;
  }

  inline bool operator()(const _move m, const typename prob_t::eval_t& e) {
    if (e<0 && e<base::thread_cost()) { 
      base::update_thread_move_and_cost(m,e);
    }
    return false;
  }

private:
  const typename prob_t::sol_t& s;
};




// this policy applies the first improving move
template <class prob_t, class _move>
struct first_improve {
  first_improve(typename prob_t::sol_t &sol, typename prob_t::eval_t &sol_eval)
    : s(sol),
      ev(sol_eval),
      _found(false)
  { }

  inline bool operator()(const _move m)
  {
    const typename prob_t::eval_t e=m.internal_cost(s);

    if (e<0) {
      m(s);    // apply improving move
      ev+=e;
      _found=true;
      return true;  // tell caller we applied the move
    }
    return false;
  }
  
  bool found() const { return _found; }
  void reset() { _found=false; }

private:
  typename prob_t::sol_t& s; 
  typename prob_t::eval_t& ev;
  bool _found;  // remember if we have found an improving move during
		// this neighborhood evaluation
};


// descent that applies first improving move
template<class prob_t, class _move, class _neighborhood, class generator_type=no_generator<prob_t> >
struct descent_fm : public meta_gen<prob_t, generator_type> {
public:
    using meta_gen<prob_t, generator_type>::operator();
    using meta_gen<prob_t, generator_type>::generator;


  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in) {
    typename prob_t::soleval_t se(se_in);
    CHECK_EVAL(se);

    _neighborhood n;
    first_improve<prob_t, _move> fm(se.first, se.second);

    do {
      fm.reset();
      n(fm, se.first);
    } while (fm.found());

    CHECK_EVAL(se);
    return se;
  }

  const std::string name() const { return "Descent that accept first improving move"; }

};



// generic DESCENT
template<class prob_t, class _move, class _neighborhood, class generator_type=no_generator<prob_t> >
class descent_base: public meta_gen<prob_t, generator_type> {

#ifdef XLC_WORKAROUND
  // IBM XLC does not let the derived classes to access the protected member.
  // this is stupid and I could not find out why.
public:
#else
protected:
#endif


  using meta_gen<prob_t, generator_type>::operator();
  using meta_gen<prob_t, generator_type>::generator;


  descent_base() { }


  template <class _nh_eval, class _gain_t>
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in, _nh_eval& nh_eval, _gain_t& gain, bool mpi=false) const {
    typename prob_t::soleval_t se(se_in);
    CHECK_EVAL(se);
 
    _move best_move;

    keep_best<prob_t, _move> keeper(se.first, best_move, mpi);
    while (1) {
      nh_eval(keeper, se.first);
      keeper.reduce();

      if (keeper.cost()>=0) break;

//       assert(gain.update_before != abstract_gain<prob_t>
//       std::cout <<" applying move: " << best_move << std::endl;
      gain.update_before(best_move, se.first);
      best_move(se.first);           // applique le mouvement
      gain.update_after(best_move, se.first);

      se.second+=keeper.cost();
  
//       keeper.print_cost();
      keeper.reset();
    }
    CHECK_EVAL(se);
    return se;
  }

};


template<class prob_t, class _move, class _neighborhood, class generator_type=no_generator<prob_t> >
class descent : public descent_base<prob_t, _move, _neighborhood, generator_type> {
  typedef descent_base<prob_t, _move, _neighborhood, generator_type> base;

public:
  using meta_gen<prob_t, generator_type>::operator();
  using meta_gen<prob_t, generator_type>::generator;

  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se) {

    _neighborhood n;
    basic_nh_eval<_neighborhood, _move, typename prob_t::sol_t> nh_eval(n);
    dummy_op dummy;

    return base::operator()(se, nh_eval, dummy);
  }
  const std::string name() const { return "Descent"; }
};




template<class prob_t, class _move, class _neighborhood, class generator_type=no_generator<prob_t> >
class descent_ns_omp : public descent_base<prob_t, _move, _neighborhood, generator_type> {
  typedef descent_base<prob_t, _move, _neighborhood, generator_type> base;
public:
  using base::operator();
  using base::generator;


  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se) {
    _neighborhood n;
    ns_nh_omp<_neighborhood, _move, typename prob_t::sol_t> nh_eval(n);
    dummy_op dummy;

    return base::operator()(se, nh_eval, dummy);
  }
  const std::string name() const { return "Parallel descent based on neighborhood seperation using OpenMP"; }
};




#ifdef USE_MPI
template<class prob_t, class _move, class _neighborhood, class generator_type=no_generator<prob_t> >
class descent_ns_mpi : public descent_base<prob_t, _move, _neighborhood, generator_type> {
  typedef descent_base<prob_t, _move, _neighborhood, generator_type> base;

public:
  using base::operator();
  using base::generator;


  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in) {

    typename prob_t::soleval_t se(se_in);
    sol_reducter<prob_t>::Allreduce(se); // make sure everyone start form same solution

    _neighborhood n;
    ns_nh_mpi<_neighborhood, _move, typename prob_t::sol_t> nh_eval(n);
    dummy_op dummy;

    return base::operator()(se, nh_eval, dummy, true);
  }
  const std::string name() const { return "Parallel descent based on neighborhood seperation using MPI"; }
};
#endif


// this version uses a gain structure to accelerate the computation of the next move
template<class prob_t, class _move, class _gain_struct_t, class generator_type=no_generator<prob_t> >
class descent_gain : public descent_base<prob_t, _move, dummy_neighborhood<prob_t>, generator_type> {
  typedef descent_base<prob_t, _move, dummy_neighborhood<prob_t>, generator_type> base;

public:
  using base::operator();
  using base::generator;


  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se) {
    
    _gain_struct_t _gain;
    _gain.init(se.first);
    gain_nh_eval<_gain_struct_t, _move, typename prob_t::sol_t> nh_eval(_gain);
    
    return base::operator()(se, nh_eval, _gain);
  }
  const std::string name() const { return "Descent using a gain structure"; }
};

}

#endif
