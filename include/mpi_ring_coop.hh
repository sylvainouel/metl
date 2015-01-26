#ifndef MPI_RING_COOP_HH
#define MPI_RING_COOP_HH
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


#ifdef USE_MPI

#include "exchange_oper.hh"

namespace metl {


template <class _metaheuristic>
class mpi_ring_coop : public _metaheuristic {
  typedef _metaheuristic base;
  typedef typename _metaheuristic::problem_type prob_t;

public:
  using base::operator();
  using base::generator;

  mpi_ring_coop(unsigned exchange_period=10, bool eval_recv=false)
    : xchange_op(exchange_period, eval_recv)
  {}


  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in) {
    typename prob_t::soleval_t se(se_in);
    CHECK_EVAL(se);

    se = base::operator()(se, xchange_op);
    sol_reducter<prob_t>::Allreduce(se);

    CHECK_EVAL(se);
    return se;
  }

  const std::string name() const 
  { 
    return "MPI ring Synchronous exchange coop. algo., base algorithm=(" + base::name()+")" ; 
  }
  

private:
  mpi_ring_coop_op<prob_t> xchange_op;

  // it is not possible to use 2 cooperative algorithm
  template <class PE>
  typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se_in, PE& periodic_exchange);
  
};


}

#endif   // USE_MPI

#endif
