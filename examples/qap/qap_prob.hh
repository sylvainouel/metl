#ifndef QAP_PROB_HH
#define QAP_PROB_HH

#include <vector>
#include <Matrix.hh>
#include "meta_base.hh"

struct qap_prob: public metl::abstract_problem<std::vector<int>, long> {
  void load(const std::string& qaplib_file);
  eval_t evaluation(const sol_t& sol) const;
  bool is_valid(const sol_t& sol) { return true; }

  inline long compute_delta(const std::vector<int>& p, unsigned i, unsigned j) const {
    eval_t d = 
      (a(i,i)-a(j,j))*
      (b(p[j],p[j])-b(p[i],p[i])) +
      (a(i,j)-a(j,i))*
      (b(p[j],p[i])-b(p[i],p[j]));

    for (unsigned k = 0; k < _size; ++k)
      if (k!=i && k!=j)
	d += (a(k,i) - a(k,j))*
	  (b(p[k],p[j]) - b(p[k],p[i])) +
	  (a(i,k) - a(j,k))*
	  (b(p[j],p[k])-b(p[i],p[k]));
    return(d);
  }

  inline long compute_delta_part(const std::vector<int>& p, unsigned i, unsigned j, unsigned r, unsigned s) const
  {
    return (a(r,i)-a(r,j)+
	    a(s,j)-a(s,i))*
      (b(p[s],p[i])-b(p[s],p[j])+
       b(p[r],p[j])-b(p[r],p[i]))+
      (a(i,r)-a(j,r)+
       a(j,s)-a(i,s))*
      (b(p[i],p[s])-b(p[j],p[s])+
       b(p[j],p[r])-b(p[i],p[r]));
  }

  inline unsigned size() const { return _size; }

  inline static qap_prob& instance() { 
    static qap_prob _instance;
    return _instance; 
  }


private:
  qap_prob();

  metl::Matrix<long> a;
  metl::Matrix<long> b;
  unsigned _size;
};


#endif
