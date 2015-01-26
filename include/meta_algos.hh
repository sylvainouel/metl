#ifndef META_ALGOS_HH
#define META_ALGOS_HH

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


#include "meta_base.hh"
#include "meta_utility.hh"
#include "abstract_move.hh"


// basic algorithms
#include "algo_sa.hh"
#include "algo_descent.hh"
#include "algo_tabu.hh"
#include "algo_evolution.hh"

// cooperative communication schemes
#include "mpi_blackboard_coop.hh"
#include "omp_blackboard_coop.hh"
#include "mpi_ring_coop.hh"
#include "omp_ring_coop.hh"
#include "mpi_reduce_coop.hh"
#include "omp_reduce_coop.hh"


#endif
