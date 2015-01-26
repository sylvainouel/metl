// Some of this code is based on taboo_qap by Eric Taillard.
// The following is the original header of taboo_qap.cpp

/*****************************************************************/
// Implementation of the robust taboo search of: E. Taillard
// "Robust taboo search for the quadratic assignment problem", 
// Parallel Computing 17, 1991, 443-455.
//
// Data file format: 
//  n,
// (nxn) flow matrix,
// (nxn) distance matrix
//
// Copyright : E. Taillard, 1990-2004
// This code can be freely used for non-commercial purpose.
// Any use of this implementation or a modification of the code
// must acknowledge the work of E. Taillard
/****************************************************************/


#include <iostream>
#include <algorithm>
#include "qap_prob.hh"
#include "meta_algos.hh"
#include "meta_main.hh"
#include "meta_permutation.hh"
#include "../test_functions.h"

using namespace metl;

// define a move for the quadratic assignment problem
struct move: public permutation_move<qap_prob> {
  move(unsigned i=0, unsigned j=1) : 
    permutation_move<qap_prob>(i,j) {}
  // constructor from an iterator of the gain structure.
  move(const utrig_matrix<qap_prob::eval_t>::iterator& it) :
    permutation_move<qap_prob>(it.get_i(),it.get_j()) {}

  // move cost evaluation function
  inline qap_prob::eval_t cost(const qap_prob::sol_t &p) const {
    return qap_prob::instance().compute_delta(p, get_i(), get_j());
  }
};  

typedef permutation_neighborhood<qap_prob, move> neighborhood;
typedef permutation_generator<qap_prob> qap_gen;


// define a tabu list for QAP
struct tabu_list: public abstract_tabu_list<qap_prob, move> {
  tabu_list(): t_list(qap_prob::instance().size(), qap_prob::instance().size()) {}

  inline bool is_tabu(const move& m, const qap_prob::sol_t& sol, unsigned current_cycle) const {
    if (t_list(m.get_i(), sol[m.get_j()]) < current_cycle ||
	t_list(m.get_j(), sol[m.get_i()]) < current_cycle) 
      return false;
    return true;
  }

  inline void make_tabu(const move& m, const qap_prob::sol_t& sol, unsigned cycle, unsigned tenur_in, unsigned tenur_out) {
    t_list(m.get_i(), sol[m.get_j()]) = cycle+tenur_in;
    t_list(m.get_j(), sol[m.get_i()]) = cycle+tenur_out;
  }


private:
  Matrix<unsigned> t_list;
};


// define a gain structure for QAP
struct gain : public abstract_gain<qap_prob, move, utrig_matrix<qap_prob::eval_t>::iterator> {

  gain()
    : G(qap_prob::instance().size(), qap_prob::instance().size()) {}
  
  void update_after(const move& m, const qap_prob::sol_t& p) {
    // move m was selected and has been performed on solution p.
    // update the gain structure
    const unsigned r = m.get_i();
    const unsigned s = m.get_j();

    // update matrix of the move costs
    const unsigned size=qap_prob::instance().size();
    const qap_prob& instance = qap_prob::instance();

    for (unsigned i = 0; i < size-1; ++i) {
      for (unsigned j = i+1; j < size; ++j)
	if (i != r && i != s && j != r && j != s)
	    G(i,j) += instance.compute_delta_part(p,i,j,r,s);
	else
	    G(i,j) = move(i,j).cost(p);
    }
  }

  inline iterator begin() {
    return G.begin();
  }

  inline iterator end() {
    return G.end();
  }

private:
  utrig_matrix<qap_prob::eval_t> G;
};




struct mutation : abstract_mutation<qap_prob> {
  void operator()(qap_prob::sol_t& s) const {
    unsigned size = qap_prob::instance().size();

    for (int i=0; i<4; ++i)
      {
	unsigned a = rng(size);
	unsigned b = rng(size);
	
	if (a!=b) 
	  std::swap(s[a],s[b]);
      }
  }
};



int _main(int argc, char* argv[])
{
  qap_prob::instance().load(argv[1]);

  typedef descent_fm<qap_prob, permutation_move<qap_prob>, permutation_neighborhood<qap_prob> ,qap_gen> ls_slow;

  typedef descent_fm<qap_prob, move, neighborhood ,qap_gen> ls;
  typedef tabu_gain<qap_prob,  move, gain, tabu_list,qap_gen> tabu_qap;
  
  typedef evolution<qap_prob, ls_slow, path_xover_insert<qap_prob>, no_mutation, ls_slow, select_random<qap_prob>, replace_worst_parent<qap_prob> > evo1_type;

  typedef evolution<qap_prob, descent_fm<qap_prob, move, neighborhood ,qap_gen>, path_xover_insert<qap_prob>, no_mutation, tabu_qap , select_random<qap_prob>, replace_worst_parent<qap_prob> > evo2_type;

  typedef simulated_annealing<qap_prob, move, neighborhood, qap_gen, metropolis<qap_prob, move> > qap_sa_t;



#ifdef USE_PAR
#ifdef USE_MPI
  typedef tabu_ns_mpi<qap_prob,  move, permutation_neighborhood<qap_prob, move>, 
    tabu_list, qap_gen> tabu_ns_mpi;
  
  tabu_ns_mpi qap_ts(qap_prob::instance().size(), 40000);
  test_generator(&qap_ts);

#else
  typedef tabu_ns_omp<qap_prob,  move, permutation_neighborhood<qap_prob, move>, 
    tabu_list, qap_gen> tabu_ns_omp;
  
  tabu_ns_omp qap_ts(qap_prob::instance().size(), 40000);
  test_generator(&qap_ts);
#endif
#else
#ifndef USE_MPI
  qap_sa_t sa;
  sa.set_init_temp(atoi(argv[2]));
  sa.set_final_temp(atoi(argv[3]));
  sa.cooling_scheme().set_cooling_factor(atof(argv[4]));
  sa.cooling_scheme().set_step_length(atoi(argv[5]));
  test_generator(&sa);

  tabu_qap qap_ts(qap_prob::instance().size(), 40000);
  test_generator(&qap_ts);
  
  evo1_type evo1(200, 2000,0);
  test_generator(&evo1);
  
  evo2_type evo2(40, 2000, 0.2);
  evo2.local_search().set_tenur(qap_prob::instance().size());
  evo2.local_search().set_n_iter(500);
  evo2.local_search().set_return_current(0.5);
  test_generator(&evo2);
#endif
  
#ifdef USE_MPI
  unsigned n = qap_prob::instance().size();
  
  mpi_reduce_coop<qap_sa_t > sa_p(atoi(argv[5]));
  
  sa_p.set_init_temp(atoi(argv[2]));
  sa_p.set_final_temp(atoi(argv[3]));
  sa_p.cooling_scheme().set_cooling_factor(atof(argv[4]));
  sa_p.cooling_scheme().set_step_length(atoi(argv[5]));
  test_generator(&sa_p);  
  
  
  mpi_blackboard_coop<tabu_qap> ts_p(200, RANDOM, 20);
  ts_p.set_tenur(rng(0.5*n, 1.5*n));
  ts_p.set_n_iter(40000);
  test_generator(&ts_p);  
  
  
  mpi_ring_coop<evo2_type> evo_p(10);
  evo_p.set_generations(2000);
  evo_p.set_popsize(40);
  
  evo_p.local_search().set_tenur(rng(0.5*n,1.5*n));
  evo_p.local_search().set_n_iter(200);
  evo_p.local_search().set_return_current(0.5);  
  
  test_generator(&evo_p);  
#endif

#endif
  return 0;
}

