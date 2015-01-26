#ifndef TWO_OPT_HH
#define TWO_OPT_HH

// ******** definition of a move for my problem ************
struct two_opt_move: public metl::abstract_move<tsp_prob> {
  two_opt_move(unsigned _a=0, unsigned _b=1, unsigned _c=1, unsigned _d=1) 
    : a(_a), b(_b), c(_c), d(_d) 
  {}
  
  tsp_prob::eval_t cost(const tsp_prob::sol_t &sol) const {
    const tsp_prob& problem = tsp_prob::instance();
    return 
      problem.dist(a,d)  + problem.dist(b,c) 
      -problem.dist(a, b) - problem.dist(c, d);
  }
  
  void operator()(tsp_prob::sol_t &sol) const {
    sol.flip(a,b,c,d);
  }

private:
  unsigned a,b,c,d;

  friend std::ostream& operator<<(std::ostream&x, const two_opt_move&m);
};


std::ostream& operator<<(std::ostream& x, const two_opt_move&m)
{
  x << m.a <<" " << m.b << " " << m.c << " " << m.d;
  return x;
}



//******** Definition of a neighborhood for my problem ************
struct two_opt_nh {
  two_opt_nh() 
    : dl_bits(tsp_prob::instance().size(),0)
    {}

  template <class _oper>
  void operator()(_oper& op, const tsp_prob::sol_t& s) {
    bool found;
    const tsp_prob& problem = tsp_prob::instance();

    for (unsigned a=0; a<s.size(); ++a) {
      found=false;
      if (dl_bits[a]) {
	continue;   
      }
      const unsigned b = s.next(a);
      const std::vector<unsigned>& clist = problem.get_candidate_list(b);
      const unsigned dab=problem.dist(a,b);
      for(std::vector<unsigned>::const_iterator j=clist.begin();
	  j!=clist.end() && problem.dist(b,*j)<dab;
	  ++j) {
	const unsigned& c = *j;
	unsigned d = s.prev(c);
	
	if (a==c || b==d) continue;
	if (op(two_opt_move(a,b,c,d))) {
	  found=true;
	  dl_bits[a]=0;
	  dl_bits[b]=0;
	  dl_bits[c]=0;
	  dl_bits[d]=0;
	  break;
	}
      }
      if (found==false) {
	const unsigned b2 = s.prev(a);
	const std::vector<unsigned>& clist2 = problem.get_candidate_list(b2);
	const unsigned dab2=problem.dist(a,b2);
	
	for(std::vector<unsigned>::const_iterator j=clist2.begin();
	    j!=clist2.end() && problem.dist(b2,*j)<dab2;
	    ++j) {
	  const unsigned& c = *j;
	  unsigned d = s.next(c);
	  if (a==c || b2==d) continue;
	  if (op(two_opt_move(b2,a,d,c))) {
	    dl_bits[a]=0;
	    dl_bits[b2]=0;
	    dl_bits[c]=0;
	    dl_bits[d]=0;
	    found=true;
	    break;
	  }
	}
      }
      if (found==false) {
	dl_bits[a]=true;
      }
    }
  }
  void reset_dl_bits() {
    fill(dl_bits.begin(), dl_bits.end(), 0);
  }
private:
  std::vector<bool> dl_bits;
};


#endif

