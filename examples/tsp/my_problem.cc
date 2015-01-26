#include <iostream>
#include <utility>

#include "tsp_prob.hh"
#include "Matrix.hh"

#include "meta_algos.hh"
#include "meta_main.hh"

#include "../test_functions.h"


#include "two_opt.hh"
#include "three_opt.hh"
#include "tsp_eax.hh"

//void boost::throw_exception(std::exception const &) {}

using namespace metl;

struct tsp_double_bridge : public abstract_mutation<tsp_prob>
{
  void operator()(tour& s) const {
    s.double_bridge();
  }
};



template <class neighborhood>
struct special_cooler: public cooling_geometric_steps<neighborhood> {
  special_cooler(double* temp)
    : cooling_geometric_steps<neighborhood>(temp)
  { }

  bool operator()() {
    if (cooling_geometric_steps<neighborhood>::operator()()) {
      // also reset dl_bits when temperature is modified
      this->nh->reset_dl_bits();
      return true;
    }
    return false;
  }
};



struct tsp_gen: public generator<tsp_prob> {
  tsp_prob::soleval_t operator()() {
    std::vector<unsigned> tmp_sol;
    const unsigned size = tsp_prob::instance().size();
    tmp_sol.reserve(size+1);
    unsigned i=0;

    for (i=0; i<size; i++)
      tmp_sol.push_back(i);
    
    //shuffle the cities in a random order
    random_shuffle(tmp_sol.begin(), tmp_sol.end(), rng);

    return std::make_pair(tour(tmp_sol),tsp_prob::instance().evaluation(tmp_sol));
  }
};



// ********** Main program *****************
int _main(int argc, char* argv[]) {
  tsp_prob::instance().load(argv[1]);

  tsp_prob::soleval_t solution=tsp_gen()();

  typedef simulated_annealing<tsp_prob, two_opt_move, two_opt_nh, tsp_gen, metropolis<tsp_prob, two_opt_move>, special_cooler<two_opt_nh> > ls2opt;

  typedef descent_fm<tsp_prob, three_opt_move, three_opt_nh, tsp_gen> descent3opt;

  typedef simulated_annealing<tsp_prob, three_opt_move, three_opt_nh, descent3opt, metropolis<tsp_prob, three_opt_move>, special_cooler<three_opt_nh> > ls3opt;



  typedef evolution<tsp_prob, descent3opt, tsp_eax, tsp_double_bridge, descent3opt, select_random<tsp_prob>, replace_worst_parent<tsp_prob> > tsp_evo;

#ifndef USE_MPI

#ifndef USE_PAR
  // using simulated annealing
//   simulated_annealing<tsp_prob, two_opt_move, two_opt_nh> tsp_sa(50, 0.1);
//   tsp_sa.cooling_scheme().set_step_length(10);
//   tsp_sa.cooling_scheme().set_cooling_factor(0.1);

//   test_metaheuristic<tsp_prob>(&tsp_sa, solution);

//   // using descent
//   //   descent<tsp_prob, two_opt_move, two_opt_nh> tsp_descent;
//   //   test_metaheuristic<tsp_prob>(&tsp_descent, solution);

//   // using descent that accept first improving move


  descent_fm<tsp_prob, two_opt_move, two_opt_nh> tsp_2opt;
  test_metaheuristic<tsp_prob>(&tsp_2opt,solution);


  /////////////

  descent_fm<tsp_prob, three_opt_move, three_opt_nh> tsp_3opt;
  test_metaheuristic<tsp_prob>(&tsp_3opt,solution);

  ////////////////


  ls2opt tsp_sa_2opt(100, 0.1);
  tsp_sa_2opt.cooling_scheme().set_step_length(1);
  tsp_sa_2opt.cooling_scheme().set_cooling_factor(0.01);
  test_generator(&tsp_sa_2opt);

  ////////////////////

  ls3opt tsp_sa_3opt(100, 0.1);
  tsp_sa_3opt.cooling_scheme().set_step_length(1);
  tsp_sa_3opt.cooling_scheme().set_cooling_factor(0.01);
  test_generator(&tsp_sa_3opt);



  /////////////////////
  //  typedef descent_fm<tsp_prob, two_opt_move, two_opt_nh> ls;
  //  typedef simulated_annealing<tsp_prob, two_opt_move, two_opt_nh> ls;

  tsp_evo tsp_evo1(10, 1000 , 0.2);
  tsp_evo1.set_childs_per_gen(4);
  
//   // configure the localsearch of the memetic algorithm
//   tsp_evo1.local_search().set_init_temp(10);
//   tsp_evo1.local_search().set_final_temp(0.1);
//   tsp_evo1.local_search().cooling_scheme().set_cooling_factor(0.1);

  test_metaheuristic<tsp_prob>(&tsp_evo1, solution);


#else


  /////////////////////////////////////////

  omp_blackboard_coop<ls3opt> tsp_ls_bb(10, RANDOM, 20);
  tsp_ls_bb.set_init_temp(100);
  tsp_ls_bb.set_final_temp(0.5);

  tsp_ls_bb.cooling_scheme().set_step_length(50);
  tsp_ls_bb.cooling_scheme().set_cooling_factor(0.05);

  test_metaheuristic<tsp_prob>(&tsp_ls_bb, solution);

  ////////////////////////////

  omp_ring_coop<ls3opt> tsp_ls_r(10);
  tsp_ls_r.set_init_temp(100);
  tsp_ls_r.set_final_temp(0.5);

  tsp_ls_r.cooling_scheme().set_step_length(50);
  tsp_ls_r.cooling_scheme().set_cooling_factor(0.05);

  test_metaheuristic<tsp_prob>(&tsp_ls_r, solution);

  ////////////////////////////

  omp_reduce_coop<ls3opt> tsp_ls_re(10);
  tsp_ls_re.set_init_temp(100);
  tsp_ls_re.set_final_temp(0.5);

  tsp_ls_re.cooling_scheme().set_step_length(50);
  tsp_ls_re.cooling_scheme().set_cooling_factor(0.05);

  test_metaheuristic<tsp_prob>(&tsp_ls_re, solution);


  ////////////////////////////////

  omp_blackboard_coop<tsp_evo> tsp_evo_bb(10, RANDOM, 20);

  tsp_evo_bb.set_popsize(30);
  tsp_evo_bb.set_generations(1000);
  tsp_evo_bb.set_mutation_rate(0.02);
  tsp_evo_bb.set_childs_per_gen(4);

  tsp_evo_bb.local_search().set_init_temp(10);
  tsp_evo_bb.local_search().set_final_temp(0.1);
  tsp_evo_bb.local_search().cooling_scheme().set_cooling_factor(0.20);

  test_metaheuristic<tsp_prob>(&tsp_evo_bb, solution);
  
  ///////////////////////////


  omp_ring_coop<tsp_evo> tsp_evo_r(10);

  tsp_evo_r.set_popsize(30);
  tsp_evo_r.set_generations(1000);
  tsp_evo_r.set_mutation_rate(0.02);
  tsp_evo_r.set_childs_per_gen(4);

  tsp_evo_r.local_search().set_init_temp(10);
  tsp_evo_r.local_search().set_final_temp(0.1);
  tsp_evo_r.local_search().cooling_scheme().set_cooling_factor(0.20);

  test_metaheuristic<tsp_prob>(&tsp_evo_r, solution);

  /////////////////////////

  omp_reduce_coop<tsp_evo> tsp_evo_re(10);

  tsp_evo_re.set_popsize(30);
  tsp_evo_re.set_generations(1000);
  tsp_evo_re.set_mutation_rate(0.02);
  tsp_evo_re.set_childs_per_gen(4);

  tsp_evo_re.local_search().set_init_temp(10);
  tsp_evo_re.local_search().set_final_temp(0.1);
  tsp_evo_re.local_search().cooling_scheme().set_cooling_factor(0.20);

  test_metaheuristic<tsp_prob>(&tsp_evo_re, solution);


#endif

#else   // USE_MPI
//   /////////////////////////////////////////

//   mpi_blackboard_coop<ls3opt> tsp_ls_bb(10, RANDOM, 20);
//   tsp_ls_bb.set_init_temp(100);
//   tsp_ls_bb.set_final_temp(0.5);

//   tsp_ls_bb.cooling_scheme().set_step_length(50);
//   tsp_ls_bb.cooling_scheme().set_cooling_factor(0.05);

//   test_metaheuristic<tsp_prob>(&tsp_ls_bb, solution);

//   ////////////////////////////

  mpi_ring_coop<tsp_evo> tsp_evo_r(1);
  tsp_evo_r.set_popsize(20);
  tsp_evo_r.set_generations(2000);
  tsp_evo_r.set_mutation_rate(0.2);

  tsp_evo_r.generator().set_init_temp(100);
  tsp_evo_r.generator().set_final_temp(1);
  tsp_evo_r.generator().cooling_scheme().set_step_length(30);
  tsp_evo_r.generator().cooling_scheme().set_cooling_factor(0.1);

  test_generator<tsp_prob>(&tsp_evo_r);

//   ////////////////////////////

  mpi_reduce_coop<ls3opt> tsp_ls_re(50);
  tsp_ls_re.set_init_temp(200);
  tsp_ls_re.set_final_temp(0.5);
  
  tsp_ls_re.cooling_scheme().set_step_length(60);
  tsp_ls_re.cooling_scheme().set_cooling_factor(0.05);

  for (int i=0; i<20; ++i)
    test_generator<tsp_prob>(&tsp_ls_re);


  mpi_blackboard_coop<ls3opt> tsp_ls_bb(50, RANDOM_BETTER);
  tsp_ls_bb.set_init_temp(200);
  tsp_ls_bb.set_final_temp(0.5);
  
  tsp_ls_bb.cooling_scheme().set_step_length(60);
  tsp_ls_bb.cooling_scheme().set_cooling_factor(0.05);

  for (int i=0; i<20; ++i)
    test_generator<tsp_prob>(&tsp_ls_bb);


//   ////////////////////////////////

  mpi_blackboard_coop<tsp_evo> tsp_evo_bb(1, RANDOM, 20);
  tsp_evo_bb.set_popsize(20);
  tsp_evo_bb.set_generations(2000);
  tsp_evo_bb.set_mutation_rate(0.2);

  tsp_evo_bb.generator().set_init_temp(100);
  tsp_evo_bb.generator().set_final_temp(1);
  tsp_evo_bb.generator().cooling_scheme().set_step_length(30);
  tsp_evo_bb.generator().cooling_scheme().set_cooling_factor(0.1);
  tsp_evo_bb.set_mutation_rate(0.2);

  test_generator<tsp_prob>(&tsp_evo_bb);
  
//   ///////////////////////////


//   mpi_ring_coop<tsp_evo> tsp_evo_r(10);

//   tsp_evo_r.set_popsize(30);
//   tsp_evo_r.set_generations(1000);
//   tsp_evo_r.set_mutation_rate(0.02);
//   tsp_evo_r.set_childs_per_gen(4);

//   tsp_evo_r.local_search().set_init_temp(10);
//   tsp_evo_r.local_search().set_final_temp(0.1);
//   tsp_evo_r.local_search().cooling_scheme().set_cooling_factor(0.20);

//   test_metaheuristic<tsp_prob>(&tsp_evo_r, solution);

//   /////////////////////////

//   mpi_reduce_coop<tsp_evo> tsp_evo_re(10);

//   tsp_evo_re.set_popsize(30);
//   tsp_evo_re.set_generations(1000);
//   tsp_evo_re.set_mutation_rate(0.02);
//   tsp_evo_re.set_childs_per_gen(4);

//   tsp_evo_re.local_search().set_init_temp(10);
//   tsp_evo_re.local_search().set_final_temp(0.1);
//   tsp_evo_re.local_search().cooling_scheme().set_cooling_factor(0.20);

//   test_metaheuristic<tsp_prob>(&tsp_evo_re, solution);


#endif  // USE_MPI
  return 0;
}
