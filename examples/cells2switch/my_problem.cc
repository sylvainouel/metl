#include "cell2switch.hh"
#include "meta_algos.hh"
#include "meta_main.hh"

#include <iostream>

#include "Matrix.hh"
#include <sys/time.h>
#include "../test_functions.h"


using namespace metl;

// incremental penality computation
float delta_penality(float cap_resi_cell, const std::vector<float>& cap_resi, unsigned cell, unsigned sw) {
  float p=0;
  const cell2switch& instance = cell2switch::instance();
  
  if (cap_resi_cell<0) {
    float p1 = cap_resi_cell + instance.get_load(cell);
    if (p1>0) {
      p = cap_resi_cell;
    }
    else {
      p = -instance.get_load(cell);
    }
  }
  
  float p2=cap_resi[sw] - instance.get_load(cell);
  if (p2 < 0) {
    if (cap_resi[sw] < 0)    
      p += instance.get_load(cell);
    else
      p+= -p2;
  }
  return 10*p;
}


// ******** definition of a move for the cell to switch assignment problem
struct move: public metl::abstract_move<cell2switch> {
  move(unsigned _c=0, int _s=0) : c(_c), s(_s)  {}
  move(const Matrix<cell2switch::eval_t>::iterator& it) 
    : c(it.get_i()), s(it.get_j())  {}

  
  cell2switch::eval_t cost(const cell2switch::sol_t &sol) const {
    if (sol[c]==s) return std::numeric_limits<cell2switch::eval_t>::max();
    // compute the cost of affecting cell c to switch s

    cell2switch::eval_t G=0;
    const cell2switch& instance = cell2switch::instance();

    for (unsigned i=0; i<instance.get_ncell(); ++i) {
      if (i==c) continue;
      if (sol[i]==sol[c])
  	G+=(instance.h_cost(c,i) + instance.h_cost(i,c));
      if (sol[i]==s)
 	G-=(instance.h_cost(c,i) + instance.h_cost(i,c));
    }
    G+=instance.c_cost(c,s);
    G-=instance.c_cost(c,sol[c]);

    // for improved performance, cap_resi could be stored together with
    // the assignation vector in a "solution structure"
    // However, it does not really matter because we are going to use
    //   a gain structure for incremental cost computation.
    std::vector<float> cap_resi;
    instance.compute_cap_resi(cap_resi, sol);

    G+=delta_penality(cap_resi[sol[c]], cap_resi, c, s);
    
    return G;
  }

  // effectue le mouvement sur la solution sol
  void operator()(cell2switch::sol_t &sol) const {
    assert(sol[c]!=s);
    sol[c]=s;
  }

  unsigned get_c() const { return c; }
  unsigned get_s() const { return s; }

private:
  //  friend class tabu_list;
  friend class gain;
//   friend std::ostream& operator<<(std::ostream& x , const move& m);
  unsigned c,s;
};


//  std::ostream& operator<<(std::ostream& x , const move& m)
//  {
//    x<<"c: " <<m.c<<" s: "<<m.s;
//    return x;
//  }


struct init_sol_gen: public metl::generator<cell2switch> {
  cell2switch::soleval_t operator()() {
    const cell2switch& instance = cell2switch::instance();
    cell2switch::sol_t s;

    for (unsigned i=0; i<instance.get_ncell(); ++i) {
      s.push_back(rng(instance.get_nswitch()));
    }
    return cell2switch::soleval_t(s, instance.evaluation(s));
  }
};



//******** Definition of a neighborhood for my problem ************
// comply with the definition of a partial neighborhood for 
// neighborhood separation parallel algorithms
class neighborhood : public metl::separable_neighborhood {
  const unsigned cells,switches;

public:
  neighborhood()
    : cells(cell2switch::instance().get_ncell()),
      switches(cell2switch::instance().get_nswitch())
  {}

  template <class _oper>
  void operator()(_oper& op, const cell2switch::sol_t& sol) {
    for (unsigned c=0; c<cells; ++c)
      iteration(op,sol,c);
  }

  template <class _oper>
  void iteration(_oper& op, const cell2switch::sol_t& sol, unsigned c) {
    for (unsigned s=0; s<switches; ++s) {
      if (sol[c]==s) continue;
      op(move(c, s));
    }
  }
  unsigned size() const { return cells; }
};



class tabu_list : public metl::abstract_tabu_list<cell2switch, move> {
public:
  tabu_list(): t_list(cell2switch::instance().get_ncell(), cell2switch::instance().get_nswitch()) {}

  inline bool is_tabu(const move& m, const cell2switch::sol_t& sol, unsigned current_cycle) const {
    return 
      current_cycle<t_list(m.get_c(), m.get_s()) ||  // interdit d'enlever le caractere
      current_cycle<t_list(m.get_c(), sol[m.get_c()]); // interdit de revenir à l'ancien 
   }

  inline void make_tabu(const move& m, const cell2switch::sol_t& sol, unsigned cycle, 
			unsigned tenur_in, unsigned tenur_out) {
    t_list(m.get_c(), m.get_s()) = cycle+tenur_in;
    t_list(m.get_c(), sol[m.get_c()]) = cycle+tenur_out;
  }
private:
   Matrix<unsigned> t_list;
};


struct gain : public 
metl::abstract_gain<cell2switch, move, Matrix<cell2switch::eval_t>::iterator> 
{
  typedef metl::abstract_gain<cell2switch, move, 
			      Matrix<cell2switch::eval_t>::iterator> base;

  gain()
    : G(cell2switch::instance().get_ncell(), 
	cell2switch::instance().get_nswitch()), 
      Gcap(G)
  {}

  void init(const cell2switch::sol_t& sol) {
    base::init(sol);
    // has to be specialized because we also have to compute the
    // vector of residual capacities and the gain without the penalities.
    cell2switch::instance().compute_cap_resi(cap_resi, sol);
    
    for (unsigned c=0; c<G.get_rows(); ++c)
      for (unsigned s=0; s<G.get_cols(); ++s) {
	  G(c,s) = Gcap(c,s)- delta_penality(cap_resi[sol[c]], cap_resi,c,s);
      }
  }

  void update_before(const move& m,const cell2switch::sol_t& sol) {
    const unsigned Acomm = sol[m.c];   // ancienne valeur du commutateur
    const unsigned comm = m.s;         // nouvelle valeur du commutateur
    const unsigned cell = m.c;    // cellule déplacée
    const unsigned nbr_cell = G.get_rows();
    const unsigned nbr_comm = G.get_cols();

    const cell2switch& instance = cell2switch::instance();

    /* mise a jour des colonnes des ancien et nouveau comm*/
    for (unsigned p=0;p<nbr_cell;p++) {
      if (p!=cell) {
	unsigned Cp = sol[p];
	if (p==m.c) Cp = m.s;

	if (Cp==Acomm) {
	  for (unsigned q=0; q<nbr_comm;q++) {
	    if (q!=Acomm && q!=comm) {
	      G(p,q)-=(instance.h_cost(p,cell)+instance.h_cost(cell,p));    // 2.23
	    }
	  }
	  G(p,comm)-=2*(instance.h_cost(p,cell)+instance.h_cost(cell,p));   // 2.22
	}
	else if (Cp==comm) {
	  for (unsigned q=0; q<nbr_comm;q++) {
	    if (q!=Acomm && q!= comm) {
	      G(p,q)+=(instance.h_cost(p,cell)+instance.h_cost(cell,p));     // 2.25
	    }
	  }
	  G(p,Acomm)+=2*(instance.h_cost(p,cell)+instance.h_cost(cell,p));  // 2.24
	}
	else { //Cp !=Acomm et !=comm
	  G(p,Acomm)+=(instance.h_cost(p,cell)+instance.h_cost(cell,p));
	  G(p,comm)-=(instance.h_cost(p,cell)+instance.h_cost(cell,p));
	}
      }
    }

    for (unsigned q=0; q<nbr_comm;q++) {
      if (q!=Acomm && q!= comm) {
	G(cell,q)-=G(cell,comm);           // 2.27
      }
    }
    G(cell,Acomm)=-G(cell,comm);
    G(cell,comm)=0;
    
    // mise à jour de cap_resi
    cap_resi[Acomm] += instance.get_load(cell);
    cap_resi[comm] -= instance.get_load(cell);

    // update Gain penality 
    
    for (unsigned c=0; c<nbr_cell; ++c) {
      const unsigned sw =          // switch apres le mouvement
	c==cell ? comm : sol[c];
      
      if (sw==Acomm || sw==comm) {
	// la cellule est affectee a Acomm ou a comm, sa penalite change
	for (unsigned s=0; s<nbr_comm; ++s)
	  {
	    if (s!=sw) {
	      if (c==cell) {
		Gcap(c,s)=G(c,s) + delta_penality(cap_resi[comm],cap_resi, c,s);
	      }
	      else 
		Gcap(c,s) = G(c,s) + delta_penality(cap_resi[sol[c]],cap_resi, c,s);
	    }
	    else
	      Gcap(c,s) = std::numeric_limits<cell2switch::eval_t>::max();
	  }
      } else {
	// tous les mouvements vers Acomm et comm changent de penalite
	if (Acomm!=sw)
	  Gcap(c,Acomm) = G(c,Acomm)+delta_penality(cap_resi[sol[c]],cap_resi,c,Acomm);
	else
	  Gcap(c,Acomm) = std::numeric_limits<cell2switch::eval_t>::max();
	if (comm!=sw)
	  Gcap(c,comm) = G(c,comm)+delta_penality(cap_resi[sol[c]],cap_resi,c,comm);
	else
	  Gcap(c,comm) = std::numeric_limits<cell2switch::eval_t>::max();
      }
    }
  }

  iterator begin() {
    return Gcap.begin();
  }

  iterator end() {
    return Gcap.end();
  }
  
private:
  Matrix<cell2switch::eval_t> G;
  Matrix<cell2switch::eval_t> Gcap;
  std::vector<float> cap_resi;  // capacite residuelle des commutateurs
};



class cell2switch_mutation: public metl::abstract_mutation<cell2switch> {
  const unsigned cells,switches;

public:
  cell2switch_mutation()
    : cells(cell2switch::instance().get_ncell()),
      switches(cell2switch::instance().get_nswitch())
  {}
  
  void operator()(solution& s) const {
    s[rng(cells)]=rng(switches);
    s[rng(cells)]=rng(switches);
  }
};



// ********** Main program *****************
int _main(int argc, char* argv[]) {
  cell2switch::instance().load(argv[1]); // load the problem from the command line

#ifndef USE_PAR
  // using descent with gain
  descent_gain<cell2switch, move, gain, init_sol_gen> 
    cell2switch_descent_g;
  test_generator(&cell2switch_descent_g);


//   // using descent that accept first improving move
  descent_fm<cell2switch, move, neighborhood, init_sol_gen> 
    cell2switch_descent_fm;
  test_generator(&cell2switch_descent_fm);

  // using tabu search  (with a gain structure)
  tabu_gain<cell2switch, move, gain, tabu_list,init_sol_gen> 
    cell2switch_ts_g(10, 1000);
  
  test_generator(&cell2switch_ts_g);

  evolution<cell2switch, descent_fm<cell2switch, move, neighborhood, init_sol_gen>, uniform_xover<cell2switch>, cell2switch_mutation, tabu_gain<cell2switch, move, gain, tabu_list> >
    cell2switch_evo2(10, 200,0.2);
  
  // configure the localsearch of the memetic algorithm
  cell2switch_evo2.local_search().set_tenur_out(10);
  cell2switch_evo2.local_search().set_tenur_in(5);
  cell2switch_evo2.local_search().set_n_iter(500);

  //  test_generator<cell2switch>(&cell2switch_evo1);
  test_generator<cell2switch>(&cell2switch_evo2);

#else

  // using parallel descent based on neighborhood separation
  descent_ns_omp<cell2switch, move, neighborhood, init_sol_gen> 
    cell2switch_descent2;
  test_generator<cell2switch>(&cell2switch_descent2);


  // using tabu search  (parallel version)
  tabu_ns_omp<cell2switch, move, neighborhood, tabu_list, init_sol_gen> 
    cell2switch_ts_omp(10, 10000);
  test_generator<cell2switch>(&cell2switch_ts_omp);


#endif

#ifdef USE_MPI

  // using descent
  descent_ns_mpi<cell2switch, move, neighborhood> 
    cell2switch_descent_mpi;
  test_generator<cell2switch>(&cell2switch_descent_mpi);


  tabu_ns_mpi<cell2switch, move, neighborhood, tabu_list> 
    cell2switch_ts_mpi(10, 10000);
  test_generator<cell2switch>(&cell2switch_ts_mpi);



  mpi_blackboard_coop<tabu_gain<cell2switch, move, gain, tabu_list> > async_exchange_ts(100,RANDOM_BETTER);


  async_exchange_ts.set_tenur_out(10);
  async_exchange_ts.set_tenur_in(5);
  async_exchange_ts.set_n_iter(10000);

  test_generator<cell2switch>(&async_exchange_ts);  


  /////////////
  typedef evolution<cell2switch, init_sol_gen, uniform_xover<cell2switch>, cell2switch_mutation, tabu_gain<cell2switch, move, gain, tabu_list> > evolution_algo;

  mpi_blackboard_coop<evolution_algo> async_exchange_evo(10,RANDOM);
  async_exchange_evo.set_popsize(30);
  async_exchange_evo.set_generations(129);
  async_exchange_evo.set_mutation_rate(0.02);
  async_exchange_evo.set_childs_per_gen(1);

  async_exchange_evo.local_search().set_tenur_out(10);
  async_exchange_evo.local_search().set_tenur_in(5);
  async_exchange_evo.local_search().set_n_iter(100);


  test_generator<cell2switch>(&async_exchange_evo);  


#endif // USE_MPI

  return 0;
}
