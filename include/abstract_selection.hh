#ifndef ABSTRACT_SELECTION_HH
#define ABSTRACT_SELECTION_HH

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

#include <vector>
#include <utility>

namespace metl {

template <class prob_t>
struct abstract_selection {
  virtual ~abstract_selection() {}

  // pop hold the popluation. It is sorted in ascending evaluation
  // order so that pop.front() is the best individual.
  //
  // the function has to set p1 the the first parent in pop and p2 to
  // the second parent.
  virtual void operator()(const typename std::vector<typename prob_t::soleval_t>& pop, 
			  typename std::vector<typename prob_t::soleval_t>::const_iterator& p1, 
			  typename std::vector<typename prob_t::soleval_t>::const_iterator& p2) const=0;
  
};
}

#endif

