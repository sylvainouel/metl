#ifndef ABSTRACT_PROBLEM_HH
#define ABSTRACT_PROBLEM_HH

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

namespace metl {

// asol_t is the data type holing the solution
// aeval_t is a scalar type returned by the evaluation function
template<class asol_t, class aeval_t>
struct abstract_problem {
  typedef asol_t sol_t;
  typedef aeval_t eval_t;
  typedef std::pair<sol_t, eval_t> soleval_t;

  virtual ~abstract_problem() {}
  virtual eval_t evaluation(const sol_t& sol) const =0;
  virtual eval_t optimum() const { return std::numeric_limits<eval_t>::min(); }
};

}
#endif
