#ifndef SYNC_RNG
#define SYNC_RNG

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


#include <sprng.h>
#include "metl_def.hh"

#ifdef _OPENMP
#ifndef INCLUDED_OMP_H   
#define INCLUDED_OMP_H   
#include <omp.h>         
#endif                   
#else                    
#include "omp_stub.h"    
#endif


#include "random_generator.hh"
namespace metl {


  // this class defines a "synchronized" rng. That is a random number generator that gives the same sequence on every threads/processes. It is a lecuyer RNG.

class sync_rng : public random_generator {
public:
  sync_rng(unsigned _x10=12345, unsigned _x11=67890, unsigned _x12=13579, unsigned _x20=24680, unsigned _x21=98765, unsigned _x22=43210) 
    : x10(_x10), x11(_x11), x12(_x12), x20(_x20), x21(_x21), x22(_x22)
  {
    sync();
  }

  double operator() () {
    /*************** L'Ecuyer random number generator ***************/
    const long m = 2147483647; const long m2 = 2145483479; 
    const long a12= 63308; const long q12=33921; const long r12=12979; 
    const long a13=-183326; const long q13=11714; const long r13=2883; 
    const long a21= 86098; const long q21=24919; const long r21= 7417; 
    const long a23=-539608; const long q23= 3976; const long r23=2071;
    const double invm = 4.656612873077393e-10;
    
    long h, p12, p13, p21, p23;
    h = x10/q13; 
    p13 = -a13*(x10-h*q13)-h*r13;
    h = x11/q12; 
    p12 = a12*(x11-h*q12)-h*r12;
    if (p13 < 0) 
      p13 = p13 + m; 
    if (p12 < 0) 
      p12 = p12 + m;
    x10 = x11; 
    x11 = x12; 
    x12 = p12-p13; 
    if (x12 < 0) 
      x12 = x12 + m;
    
    h = x20/q23; 
    p23 = -a23*(x20-h*q23)-h*r23;
    h = x22/q21; 
    p21 = a21*(x22-h*q21)-h*r21;
    
    if (p23 < 0) 
      p23 = p23 + m2; 
    if (p21 < 0) 
      p21 = p21 + m2;
    
    x20 = x21; 
    x21 = x22; 
    x22 = p21-p23; 
    if(x22 < 0) 
      x22 = x22 + m2;
    if (x12 < x22) 
      h = x12 - x22 + m; 
    else h = x12 - x22;
    
    if (h == 0) 
      return(1.0); 
    else 
      return(h*invm);
  }

  unsigned operator()(unsigned max) {
    unsigned x = unsigned((*this)()*std::numeric_limits<unsigned>::max()) % max;
    //  assert(x<max);
    return x;
  }

  unsigned operator()(unsigned low, unsigned high) { 
    return low + unsigned(double(high-low)*operator()()); 
  }

  void sync()
  {
#ifdef USE_MPI
    MPI::COMM_WORLD.Bcast(&x10, 1, MPI::LONG, 0);
    MPI::COMM_WORLD.Bcast(&x11, 1, MPI::LONG, 0);
    MPI::COMM_WORLD.Bcast(&x12, 1, MPI::LONG, 0);
    MPI::COMM_WORLD.Bcast(&x20, 1, MPI::LONG, 0);
    MPI::COMM_WORLD.Bcast(&x21, 1, MPI::LONG, 0);
    MPI::COMM_WORLD.Bcast(&x22, 1, MPI::LONG, 0);
#endif
  }

private:
  long x10, x11, x12, x20, x21, x22;
};

}

#endif


