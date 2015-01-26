#ifndef MOVE_REDUCTION_HH
#define MOVE_REDUCTION_HH

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

#include <limits>
#include <utility>
#include "metl_def.hh"

#ifdef _OPENMP
#ifndef INCLUDED_OMP_H   
#define INCLUDED_OMP_H   
#include <omp.h>         
#endif                   
#else                    
#include "omp_stub.h"    
#endif

#include "mpi_reduce_op.hh"

namespace metl {


// this is a base class. Each move selection policies that needs a reduction at the end should inherit from this class.
template<class prob_t, class _move>
class move_reduction {
private:
  std::vector<_move> best_moves;
  std::vector<typename prob_t::eval_t> _bm_costs;
  const bool _using_mpi;

#ifdef USE_MPI
  MPI::Op pair_reduce_op;
  MPI::Datatype type_pair;
#endif  
  _move& _bm;

protected:
  typename prob_t::eval_t _bm_cost;


  move_reduction(_move& best_move, bool using_mpi)
    : best_moves(omp_get_max_threads(), _move()),  //REQUIS: move must be default constructible
      _bm_costs(omp_get_max_threads(), std::numeric_limits<typename prob_t::eval_t>::max()),
      _using_mpi(using_mpi),
#ifdef USE_MPI
      pair_reduce_op(),
      type_pair(),
#endif
      _bm(best_move),
      _bm_cost(std::numeric_limits<typename prob_t::eval_t>::max())
  {
#ifdef USE_MPI
    pair_reduce_op.Init(pair_reduce<std::pair<_move,typename prob_t::eval_t> >, 0);
    type_pair = MPI::INT.Create_contiguous(sizeof(std::pair<_move,typename prob_t::eval_t>)/sizeof(int));
    type_pair.Commit();
#endif
  }

  virtual ~move_reduction() {
#ifdef USE_MPI
    pair_reduce_op.Free();
    type_pair.Free();
#endif    
  }

  inline void update_thread_move_and_cost(const _move& m, const typename prob_t::eval_t e) {
    if (omp_get_max_threads()>1) {
      _bm_costs[omp_get_thread_num()] = e;
      best_moves[omp_get_thread_num()] = m;
    } else {
      _bm_cost = e;
      _bm = m;
    }
  }

  inline typename prob_t::eval_t thread_cost() const {
    return omp_get_max_threads()==1 ? _bm_cost : _bm_costs[omp_get_thread_num()];
  }

public:
  // should always call reduce() before accessing the move or calling cost
  inline void reduce() {
    if (omp_get_max_threads()>1) {

      _bm_cost=_bm_costs[0];
      unsigned min_i=0;
      
      for (int i=1; i<omp_get_max_threads(); ++i) {
	if (_bm_costs[i] < _bm_cost) {
	  _bm_cost = _bm_costs[i];
	  min_i = i;
	}
      }
      _bm = best_moves[min_i];
    }

#ifdef USE_MPI
    if (_using_mpi) {
      std::pair<_move,typename prob_t::eval_t> x;
      std::pair<_move,typename prob_t::eval_t> y;
      x.first = _bm;
      x.second=_bm_cost;

      // FIXME: should use a serialization

      //    cerr << "type_pair size:" << type_pair.Get_size() << endl;
      //      std::cout << MPI::COMM_WORLD.Get_rank() << " before: " << x.second << " " << x.first << "    ";
      
      MPI::COMM_WORLD.Allreduce(&x,&y, 1, type_pair, pair_reduce_op);
      //      std::cout << MPI::COMM_WORLD.Get_rank() << " after: " << y.second << " " << y.first << std::endl;
      _bm = y.first;
      _bm_cost = y.second;

    }
#endif
  }
  
  // forget the current best move
  void reset() {
    _bm_cost = std::numeric_limits<typename prob_t::eval_t>::max();
    if (omp_get_max_threads()>1) {
      for (unsigned i=0; i<_bm_costs.size(); ++i) { 
	_bm_costs[i] = std::numeric_limits<typename prob_t::eval_t>::max();
      }
    }
  }

  // return the cost of the selected move
  inline typename prob_t::eval_t cost() const { return _bm_cost; } 

  virtual bool operator()(const _move m)=0;
};

}
#endif
