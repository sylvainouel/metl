#include "cell2switch.hh"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;


void cell2switch::load(const string& f)
{
  ifstream fichier;

  fichier.open((f+".don").c_str());
  if (fichier.fail()) {
    fichier.open((f+".dat").c_str());
    if (fichier.fail()) {
      cerr<<"Fichier de donnees "<< f+".don" << " non trouve\n";
      exit (1);
    }
  } 

  fichier>>nbr_cell>>nbr_comm;
  
  // FIXME: this cause an useless copy
  cable_cost = metl::Matrix<float>(nbr_cell, nbr_comm);
  handover_cost = metl::Matrix<float>(nbr_cell, nbr_cell);
  
  for (unsigned i=0; i<nbr_cell;i++) {
    for (unsigned j=0;j<nbr_comm;j++) {
      fichier>> cable_cost(i,j) ;
    }
  }
  for (unsigned i=0; i<nbr_cell;i++) {
    for (unsigned j=0;j<nbr_cell;j++) {
      fichier>>handover_cost(i,j) ;
    }
  }
  fichier.close();
  fichier.open((f+".cap").c_str());

  if (fichier.fail()) {
    cout<<"Fichier de capacite" << f+".cap" << " non trouve\n";
    exit(1);
  }
  /*lecture des donnees*/
  float cap;
  for (unsigned i=0; i<nbr_cell;i++) {    
    fichier>> cap;
    cell_load.push_back(cap);
  }
  for (unsigned i=0; i<nbr_comm;i++) {
    fichier>> cap;
    cap_switch.push_back(cap);
  }
  fichier.close();
}


cell2switch::eval_t cell2switch::evaluation(const sol_t& sol) const 
{
  cell2switch::eval_t acc=0;

  for (unsigned i=0; i<nbr_cell; ++i) {
    acc+=cable_cost(i,sol[i]);
    for (unsigned j=0; j<nbr_cell; ++j)
      if (sol[i] != sol[j])
	acc+=handover_cost(i,j);
  }

  acc += penality(sol);
  return acc;
}


cell2switch::eval_t cell2switch::penality(const sol_t& sol) const
{
  vector<float> cap_resi;
  compute_cap_resi(cap_resi, sol);

  float P=0;
  for (std::vector<float>::const_iterator i=cap_resi.begin(); i!=cap_resi.end(); ++i) {
    if (*i < 0)
      P+= -(*i);
  }

  return 10*P;
}


void cell2switch::compute_cap_resi(vector<float>& cap_resi, const sol_t& sol) const {
  cap_resi = cap_switch;

  for (sol_t::const_iterator i = sol.begin(); i!=sol.end(); ++i) {
    cap_resi[*i] -= cell_load[i-sol.begin()];
  }
}


bool cell2switch::is_valid(const sol_t& sol) const {
  return penality(sol)<0.000010;
}
