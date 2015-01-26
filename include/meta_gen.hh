#ifndef META_GEN_HH
#define META_GEN_HH

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

#include "metaheuristic.hh"
#include "generator.hh"
#include <utility>

namespace metl {

// this is base class for a metaheuristic that is also a generator
template <class prob_t, class generator_type>
class meta_gen : public metaheuristic<prob_t> {
protected:
  typedef prob_t problem_type;

  meta_gen() 
    : gen() 
  {}
  

public:
  generator_type& generator() { return gen; }

   // don't know initial evaluation
  virtual typename  prob_t::soleval_t operator()(const typename prob_t::sol_t& s) {
    return operator()(std::make_pair(s, prob_t::instance().evaluation(s)));
  }
  
  // initial solution evaluation is known
  virtual typename prob_t::soleval_t operator()(const typename prob_t::soleval_t& se) =0;
  
  virtual typename prob_t::soleval_t operator()() 
  {
    // applique la recherche sur la solution generee
    return operator()(gen());
  }

protected:
  generator_type gen;
};

}
#endif
