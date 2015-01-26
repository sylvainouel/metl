#ifndef PERMUTATION_GENERATOR_HH
#define PERMUTATION_GENERATOR_HH

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


// this is a generic generator for simple permutation problems

#include "generator.hh"

namespace metl {

template <class prob_t>
struct permutation_generator : generator<prob_t>
{
  typename prob_t::soleval_t operator()() {
    typename prob_t::sol_t s;

    const prob_t& instance=prob_t::instance();

    // FIXME: for this to work, sol_t has to be a vector.
    //        but why would it not be a vector anyway ?
    s.reserve(instance.size());
    for (unsigned i=0; i<instance.size(); ++i)
      s.push_back(i);

    std::random_shuffle(s.begin(), s.end(), rng);
    return make_pair(s, instance.evaluation(s));
  }
};

}
#endif
