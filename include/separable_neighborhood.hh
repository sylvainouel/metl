#ifndef SEPARABLE_NEIGHBORHOOD_HH
#define SEPARABLE_NEIGHBORHOOD_HH

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

struct separable_neighborhood {
  virtual ~separable_neighborhood() {}
  // It is not possible to have template virtual functions. Users must
  // still implement these functions in their derivated classes.

  // template <class _oper>
  // virtual void operator()(_oper& op, const sol_t& sol) const=0;

  //  template <class _oper>
  //  virtual void iteration(_oper& op, const sol_t& sol, unsigned i) const=0;
  

  virtual unsigned size() const =0;
};

}

#endif
