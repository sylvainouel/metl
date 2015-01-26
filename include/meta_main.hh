#ifndef META_MAIN_HH
#define META_MAIN_HH

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

#include "metl_def.hh"

#ifdef USE_MPI
#ifdef HAVE_MPIPP_H     
#include <mpi++.h>
#else
#include <mpi.h>
#endif  // HAVE_MPIPP_H
#endif  // USE_MPI


#include <iostream>

#ifdef BOOST_NO_EXCEPTIONS
namespace boost {
  void throw_exception(std::exception const & e) {}
}
#endif


namespace metl {
sprng_rand rng;
}

extern int _main(int argc, char* argv[]);

int main(int argc, char* argv[]) 
{
#ifndef NDEBUG
  std::cout << "Compiled without NDEBUG, will be slow." << std::endl;
#endif
#ifdef USE_MPI
  MPI::Init(argc, argv);
#endif
  metl::rng.init();
#ifdef _OPENMP
  std::cout << "Using OpenMP. omp_get_max_threads=" << omp_get_max_threads() << std::endl;
#endif

  _main(argc, argv);

#ifdef USE_MPI
  MPI::Finalize();
#endif
  return 0;
}
#endif
