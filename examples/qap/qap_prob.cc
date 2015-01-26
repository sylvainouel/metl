#include <fstream>
#include <iostream>
#include "qap_prob.hh"
#include <meta_utility.hh>

using namespace std;
using namespace metl;


qap_prob::qap_prob()
  : a(), b(), _size(0)
 {}

void qap_prob::load(const string& qaplib_file) 
{

  /************** read file name and problem size ***************/
  cout << "Data file name : " << qaplib_file.c_str() << std::endl;

  ifstream data_file(qaplib_file.c_str());
  data_file >> _size;

  a = Matrix<long>(_size,_size);
  b = Matrix<long>(_size,_size);
  
  /************** read flows and distances matrices **************/
  for (unsigned i = 0; i < _size; ++i) 
    for (unsigned j = 0; j < _size; ++j)
      data_file >> a(i,j);
  for (unsigned i = 0; i < _size; ++i) 
    for (unsigned j = 0; j < _size; ++j)
      data_file >> b(i,j);
  data_file.close();
}


qap_prob::eval_t qap_prob::evaluation(const sol_t& sol) const
{
  eval_t current_cost = 0;
  for (unsigned i = 0; i < _size; ++i) 
    for (unsigned j = 0; j < _size; ++j)
      {
	current_cost += a(i,j) * b(sol[i],sol[j]);
      }
  return current_cost;
}

