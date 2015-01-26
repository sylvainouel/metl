#ifndef NEIGHBORHOOD_OPER_HH
#define NEIGHBORHOOD_OPER_HH

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

// this defines some operations that are used by the neighborhood based searches

#include "metl_def.hh"

#ifdef USE_MPI
#ifdef HAVE_MPIPP_H     
#include <mpi++.h>
#else
#include <mpi.h>
#endif  // HAVE_MPIPP_H
#endif  // USE_MPI

#ifdef _OPENMP
#ifndef INCLUDED_OMP_H   
#define INCLUDED_OMP_H   
#include <omp.h>         
#endif                   
#else                    
#include "omp_stub.h"    
#endif


namespace metl {


// ===========================================
// operations that search a given neighborhood
// ===========================================

// this is a basic neighborhood evaluation. Evaluate whole neighborhood
template <class _nh_t, class _move, class sol_t>
struct basic_nh_eval {
  basic_nh_eval(_nh_t& neighborhood) : nh(neighborhood) {};
  template<class _op>
  inline void operator()(_op& op, const sol_t& sol) {
    nh(op, sol);        // evaluate whole neighborhood
  }
private:
  _nh_t& nh;
};


// Evaluate whole neighborhood using OMP to split the work on multiple processors using shared memory
template <class _nh_t, class _move, class sol_t>
struct ns_nh_omp {
  ns_nh_omp(_nh_t& neighborhood) : nh(neighborhood) {};
  template<class _op>
  inline void operator()(_op& op, const sol_t& sol) {
    const int size = static_cast<int>(nh.size());
    int i;

#pragma omp parallel for schedule(guided)
    for (i=0; i<size; ++i) {
      nh.iteration(op, sol, i);        // evaluate partial neighborhoods
    }
  }
private:
  _nh_t& nh;
};


#ifdef USE_MPI
// Evaluate whole neighborhood using MPI to split the work on multiple processors using message passing
template <class _nh_t, class _move, class sol_t>
struct ns_nh_mpi {
  ns_nh_mpi(_nh_t& neighborhood) : nh(neighborhood) {};
  template<class _op>
  inline void operator()(_op& op, const sol_t& sol) {
    const unsigned rank = MPI::COMM_WORLD.Get_rank();
    const unsigned size = MPI::COMM_WORLD.Get_size();
    const unsigned ns = nh.size();

    // this is "alternated" interlaced work distribution.

    // size=2
    //p0: 0     3 4     7 8      11
    //p1:   1 2     5 6     9 10

    //size=3
    //p0: 0         5 6          11
    //p1:   1     4     7     10
    //p2:     2 3         8 9

    unsigned i=rank;
    unsigned j=1;
    while (i<ns) {
      //      std::cerr << rank << "   " << i << std::endl;
      nh.iteration(op, sol, i);        // evaluate partial neighborhoods
      if (j & 1) {
	i=(j+1)*size-1-rank;
      } else {
	i=j*size+rank;
      }
      ++j;
    }

    // this is interlaced work distribution
//     for (unsigned i=rank; i<ns; i+=size) {
//       nh.iteration(op, sol, i);        // evaluate partial neighborhoods
//     }

  }
private:
  _nh_t& nh;
};
#endif


// this version uses the gain structure to evaluate the cost of every moves
template <class _gain_t, class _move, class sol_t>
struct gain_nh_eval {
  gain_nh_eval(_gain_t& gain) : g(gain) {};

  template<class _op>
  void operator()(_op& op, const sol_t& sol) {
    // search the gain structure for the best possible move
    const typename _gain_t::iterator gend = g.end();
    for (typename _gain_t::iterator i = g.begin(); i!=gend; ++i) {
      op(_move(i), *i);   // cost of the move is already available

#ifndef NDEBUG
      if (fabs(*i - _move(i).internal_cost(sol))>0.01) {
	std::cout << "bad gain : " << " gain:" << *i << " expected: " <<_move(i).cost(sol) << std::endl;
      }
#endif

    }
  }
private:
  _gain_t& g;
  gain_nh_eval(const gain_nh_eval&);
  gain_nh_eval& operator=(const gain_nh_eval&);
};

}
#endif

