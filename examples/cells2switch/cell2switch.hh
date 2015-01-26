#ifndef CELL2SWITCH
#define CELL2SWITCH

#include <vector>
#include <string>

#include "meta_base.hh"
#include "Matrix.hh"

typedef std::vector<unsigned> solution;

class cell2switch: public metl::abstract_problem<solution, double> {
public:
  double evaluation(const solution& sol) const;

  static cell2switch& instance() {
    static cell2switch _instance;
    return _instance;
  }

  void load(const std::string& f);

  unsigned get_ncell() const { return nbr_cell; }
  unsigned get_nswitch() const { return nbr_comm; }
  float c_cost(unsigned c, unsigned s) const { return cable_cost(c,s); }
  float h_cost(unsigned c1, unsigned c2) const { return handover_cost(c1,c2); }
  float get_load(unsigned i) const { return cell_load[i]; }


  void compute_cap_resi(std::vector<float>& cap_resi, const solution& sol) const;
  bool is_valid(const solution& sol) const;


private:
  cell2switch() {}
  double penality(const solution& sol) const;

  metl::Matrix<float> cable_cost;
  metl::Matrix<float> handover_cost;

  std::vector<float> cell_load;
  std::vector<float> cap_switch;

  unsigned nbr_cell;
  unsigned nbr_comm;
};

#endif
