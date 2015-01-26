#ifndef ABSTRACT_REPLACE_HH
#define ABSTRACT_REPLACE_HH

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
struct abstract_replace {
  virtual ~abstract_replace() {}

  // pop is the population vector (sorted)
  // childs is the childrens vector (not sorted)
  // wp is an iterator that points to the worst parent. Only usefull
  //   if you want to replace the worst parent otherwise ignore it.
  virtual void operator()(std::vector<typename prob_t::soleval_t>& pop, 
			  std::vector<typename prob_t::soleval_t>& childs,
			  const typename std::vector<typename prob_t::soleval_t>::iterator& wp) 
    const=0;
};

}
#endif

