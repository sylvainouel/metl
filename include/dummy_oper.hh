#ifndef DUMMY_OPER_HH
#define DUMMY_OPER_HH

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

namespace metl {

  // this is a dummy generator that cannot even be called.
  // I would have liked to make it a compile time error but I could not find a way.
template <class prob_t>
struct no_generator: public generator<prob_t> 
{
  typename prob_t::soleval_t operator()()
  {
    std::cerr << "invalid call to no_generator::operator()" << std::endl;
    abort();
    return typename prob_t::soleval_t();
  }
}; 



// this is a dummy operation that should fit as a place holder for
// most optional operations
template <class useless>
struct dummy_op_t {
  template <class a1>
  inline bool operator()(const a1& a) const { return false; }

  template <class a1, class a2>
  inline bool operator()(const a1& a, const a2& b) const { return false; }

  template <class a1, class a2>
  inline void update_before(const a1& a, const a2& b) const { }

  template <class a1, class a2>
  inline void update_after(const a1& a, const a2& b) const { }

  template <class crap>
  inline void reset(const crap& c) const {}

  template <class crap>
  inline void set_rate(const crap& c) const {}

  inline float rate() const { return 0; }

  template <class a1>
  inline void send(const a1& a) const {}

  template <class a1, class b1>
  inline void send(const a1& a,const b1& b) const {}

  template <class a1>
  inline bool recv(const a1& a) const { return false; }

  template <class a1>
  inline void init(const a1& a) const {}

};

typedef dummy_op_t<int> dummy_op;
typedef dummy_op_t<int> no_mutation;
#define dummy_neighborhood dummy_op_t


// this is a dummy local search
// this one really is specific.
template <class prob_t>
struct no_localsearch {
  const typename prob_t::soleval_t& operator()(const typename prob_t::soleval_t& se) {
    return se;  // this should allow assignment to self to be optimized away
  }

  const std::string name() const { return "no local search"; }
};

}

#endif
