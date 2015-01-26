#ifndef PERMUTATION_MOVE_HH
#define PERMUTATION_MOVE_HH

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


#include "abstract_move.hh"
#include <algorithm>

namespace metl {

// this is a special move that is a permutation of 2 elements in a vector.  
// for permutation problems, user can inherit from this class and only
// define the move cost function (cost())
template <class prob_t>
struct permutation_move : public abstract_move<prob_t> {
  permutation_move(unsigned _i=0, unsigned _j=1) : i(_i), j(_j) {}

  // move application function
  inline void operator()(typename prob_t::sol_t& sol) const {
    //    std::cout << "move " << i << " " <<j << std::endl;
    std::swap(sol[i], sol[j]);
  }
 
  inline unsigned get_i() const { return i; }
  inline unsigned get_j() const { return j; }

private:
  unsigned i,j;
};  


}

#endif
