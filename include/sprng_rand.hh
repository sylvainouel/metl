#ifndef SPRNG_RAND
#define SPRNG_RAND

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

#ifdef _OPENMP
#ifndef INCLUDED_OMP_H   
#define INCLUDED_OMP_H   
#include <omp.h>         
#endif                   
#else                    
#include "omp_stub.h"    
#endif


#include <sprng.h>
#include "random_generator.hh"

namespace metl {

class sprng_rand: public random_generator {
 public:
  sprng_rand(int gen_type=0)
    : _gen_type(gen_type)
  {
    for (int i=0; i<omp_get_max_threads(); ++i) {
      stream[i]=0;
    }
  }

  void init() {
#ifdef USE_MPI
    int streamnum = MPI::COMM_WORLD.Get_rank();
    int nstreams = MPI::COMM_WORLD.Get_size();
#else
    int streamnum=0;
    int nstreams=1;
#endif

    for (int i=0; i<omp_get_max_threads(); ++i) {
      int seed = make_sprng_seed();
      stream[i] = init_sprng(_gen_type,streamnum,nstreams,seed,SPRNG_DEFAULT);
    }
  }

  ~sprng_rand() {
    for (int i=0; i<omp_get_max_threads(); ++i) {
      if (stream[i]!=0)
	free_sprng(stream[i]);
    }
  }

  // return a random number [0,1)
  inline double operator()() {
    return sprng(stream[omp_get_thread_num()]);
  }

  // return a random number between 0 and N
  inline unsigned operator()(unsigned N) {
    return static_cast<unsigned>(N*sprng(stream[omp_get_thread_num()]));
  }

  inline unsigned operator()(unsigned _min, unsigned _max) {
    return operator()(_max-_min)+_min;
  }

 private:
  int* stream[MAX_THREADS];
  const int _gen_type;

  sprng_rand& operator=(const sprng_rand&);
  sprng_rand(const sprng_rand&);

};


extern sprng_rand rng;

}

#endif

