#ifndef TSP_PROB_HH
#define TSP_PROB_HH

#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <cmath>
#include <deque>
#include <limits>
#include "tour.hh"
#include "candidate_lists.hh"

#include "meta_base.hh"


#define MAX_MATRIX_SIZE -1

// this is incomplet, TSPLIB define other types and formats
typedef enum {TYPE_NONE, EUC_2D, CEIL_2D, EXPLICIT, ATT} t_weight_type;
typedef enum {FORMAT_NONE, FULL_MATRIX, LOWER_DIAG_ROW, UPPER_DIAG_ROW, UPPER_ROW} t_weight_format;


class tsp_prob: public metl::abstract_problem<tour, int> {
public:
  void load(const std::string& tspfile);

  // evaluate a solution, returns solution cost.
  int evaluation(const tour& sol) const;
  unsigned size() const { return prob_size; }
  
  // returns distance between city i and city j
  inline unsigned dist(unsigned i, unsigned j) const { 
    if (MAX_MATRIX_SIZE==-1) {
      return d[i][j];
    } // else {
//       if (i<maxrow) return d[i][j];
//       if (j<maxrow) return d[j][i];   // symetric
//       const double dx = cities[i].first  - cities[j].first;
//       const double dy = cities[i].second - cities[j].second;

//       switch (weight_type) {
//       case EUC_2D: {
// 	const double euc_d = sqrt(dx*dx+dy*dy);
// 	return (int)(euc_d+0.5);   // round to nearest integer
//       }
//       case CEIL_2D: {
// 	const double euc_d = sqrt(dx*dx+dy*dy);
// 	return (int)(ceil(euc_d));  // CEIL_2D
//       }
//       case ATT: {
// 	const double euc_d = sqrt((dx*dx+dy*dy)/10.0);
// 	const int tij = (int)(euc_d+0.5);
// 	if (tij<euc_d) 
// 	  return tij+1;
// 	else 
// 	  return tij;
//       }
//       default:
// 	return 0;  // should not happend
//       }
//     }
  }

//   // return the cost of inserting cityC between cityA and cityB
//   inline int insert_diff(int cityA, int cityB, int cityC) const {
//     return -dist(cityA,cityB) + dist(cityA,cityC) + dist(cityB,cityC);
//   }


  bool is_valid(const tour& sol) const;
  void plot_sol(const tour& sol, std::ostream& x) const;  // work only if weigth_type=EUC_2D | CEIL_2D
//   void print_sol(const vector<int>& sol, ostream& x) const;  

//   // generate a random solution for the problem.
  void canonical_sol(tour& sol) const;


  const std::vector<unsigned>& get_candidate_list(unsigned i) const {
    return candidate[i];
  }

//   int tsp_bestfirst(tour& sol) const;
//   int tsp_glouton_rand(tour& sol, float rand_param) const;

  static tsp_prob& instance() {
    static tsp_prob _instance;
    return _instance;
  }

private:
  tsp_prob();
  ~tsp_prob();

  int** d;                // distance matrix
  unsigned prob_size;          // problem size
  unsigned maxrow;

  std::vector<std::pair<double,double> > cities;
  t_weight_type weight_type;

  candidate_lists candidate;
  // forbid copy contruction
  tsp_prob(const tsp_prob&);

  void gready_subtours_recombine(std::vector<std::deque<int> > &subtours) const;
  //  double compute_diversity(const vector<Chromo>& pop) const;
  tsp_prob& operator=(const tsp_prob&);
};
#endif
