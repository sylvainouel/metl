#ifndef PERMUTATION_NEIGHBORHOOD_HH
#define PERMUTATION_NEIGHBORHOOD_HH

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


#include "permutation_move.hh"
#include "separable_neighborhood.hh"

namespace metl {

  // this is a pre-defined neighborhood for permutation problems
  // to use this neighborhood, your problem must define a size function.
template <class prob_t, class move=permutation_move<prob_t> >
struct permutation_neighborhood: separable_neighborhood {
  permutation_neighborhood()
    : _size(prob_t::instance().size())
 {}

  template <class _oper>
  void operator()(_oper& op, const typename prob_t::sol_t& s) {
    for (unsigned i=0; i<_size-1; ++i)
      iteration(op,s,i);
  }

  template <class _oper>
  void iteration(_oper& op, const typename prob_t::sol_t& s, unsigned i) {
    for (unsigned j=i+1; j<_size; ++j)
      op(move(i,j));
  }

  unsigned size() const { return _size; };

private:
  const unsigned _size;
};


}

#endif
