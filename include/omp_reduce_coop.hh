#ifndef OMP_REDUCE_COOP_HH
#define OMP_REDUCE_COOP_HH

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

#ifdef _OPENMP
#ifndef INCLUDED_OMP_H   
#define INCLUDED_OMP_H   
#include <omp.h>         
#endif                   
#else                    
#include "omp_stub.h"    
#endif


#include "exchange_oper.hh"

#include <limits>

namespace metl {


template <class _metaheuristic>
class omp_reduce_coop : public _metaheuristic {
  typedef _metaheuristic base;
  typedef typename _metaheuristic::problem_type prob_t;

public:
    using base::operator();
    using base::generator;

  omp_reduce_coop(unsigned exchange_period=10, bool eval_recv=false)
    : shared_individu(typename prob_t::sol_t(), std::numeric_limits<typename prob_t::eval_t>::max()),
      xchange_op(exchange_period, &shared_individu, eval_recv)
  {}

  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in) {
    typename prob_t::soleval_t se(se_in);
    CHECK_EVAL(se);

#pragma omp parallel 
    {
      typename prob_t::soleval_t tsol = se;  // REQUIS: sol must be copy constructible
      omp_reduce_coop<_metaheuristic> thread_copy(*this);
      tsol = thread_copy.runme(tsol);
#pragma omp critical
      {
	if (tsol.second < se.second) {
	  se = tsol;
	}
      }
    }

    CHECK_EVAL(se);
    return se;
  }

  const std::string name() const { return "OpenMP reduce cooperative algorithme., base algorithm=(" + base::name()+")" ; }
  

private:
  typename prob_t::soleval_t shared_individu;
  
  omp_reduce_coop_op<prob_t> xchange_op;

  typename prob_t::soleval_t runme(const typename prob_t::soleval_t& sol) {
    return base::operator()(sol, xchange_op);
  }


  // it is not possible to use 2 cooperative algorithm
  template <class PE>
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in, PE& periodic_exchange);
 
};

}
#endif
