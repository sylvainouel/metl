#ifndef ABSTRACT_MOVE_HH
#define ABSTRACT_MOVE_HH

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

#include <iostream>
#include <limits>

namespace metl {

template <class prob_t>
struct abstract_move {
  virtual ~abstract_move() {}; 

  // this is the function that is used internally. In debug mode, It
  // checks that the returned value of user defined cost function is
  // consistant with the evaluation cost.
  inline typename prob_t::eval_t internal_cost(const typename prob_t::sol_t &sol) const {
    const typename prob_t::eval_t c = cost(sol);
#ifndef NDEBUG
    const typename prob_t::eval_t evc = evaluation_cost(sol);
    if (!(fabs(evc - c)<0.1)) {
      std::cerr << "move evaluation returns wrong result: "<< c << " expected: " << evc << std::endl;
      assert(0);
    }
#endif
    return c;
  }
  

  // apply move on sol
  virtual void operator()(typename prob_t::sol_t &sol) const =0;

private:

  // evaluate move using the evaluation function, this is very slow
  typename prob_t::eval_t evaluation_cost(const typename prob_t::sol_t &sol) const {
    typename prob_t::sol_t tmp_sol = sol;
    (*this)(tmp_sol);
    // moves that does nothing should never be considered.
    if (sol == tmp_sol) 
      return std::numeric_limits<typename prob_t::eval_t>::max();

    const prob_t& instance = prob_t::instance();

    return instance.evaluation(tmp_sol) - instance.evaluation(sol);
  }



  // this is the user defined cost. By default, it uses the evaluation cost. User should override this function as soon as possible.
  virtual typename prob_t::eval_t cost(const typename prob_t::sol_t &sol) const {
    return evaluation_cost(sol);
  }

};

}
#endif
