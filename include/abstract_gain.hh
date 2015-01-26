#ifndef ABSTRACT_GAIN_HH
#define ABSTRACT_GAIN_HH

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

namespace metl {

template<class prob_t, class _move, class _iterator_type>
struct abstract_gain {
  typedef _iterator_type iterator;

  virtual ~abstract_gain() {}

  virtual void init(const typename prob_t::sol_t& sol) {
    const iterator gend = end();
    for (iterator i=begin(); i!=gend; ++i)
      *i = _move(i).internal_cost(sol);
  }


  // this one is called before the move is made
  // it is so that all information is available because sometime a move is not reversible
  // be updating the gain structure before making the move, all information is always available
  virtual void update_before(const _move& m, const typename prob_t::sol_t& sol) {};
  // if the move is reversible, the implementer could override this function instead. This one is called after the move is made.
  virtual void update_after(const _move& m, const typename prob_t::sol_t& sol) {};

  virtual iterator begin() =0;
  virtual iterator end() =0;
};
}

#endif
