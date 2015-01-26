#ifndef THREE_OPT_HH
#define THREE_OPT_HH



struct three_opt_move: public metl::abstract_move<tsp_prob> {
  three_opt_move(unsigned _a=0, unsigned _b=0, unsigned _c=0, unsigned _d=0, unsigned _e=0, unsigned _f=0) 
    : A(_a), B(_b), C(_c), D(_d), E(_e), F(_f), which(0) 
  {}
  
  tsp_prob::eval_t cost(const tsp_prob::sol_t &sol) const {
    if (C==E || A==E) {
      which=4;
      if(B==C) return std::numeric_limits<tsp_prob::eval_t>::max();
      return two_opt_move(A,B,D,C).cost(sol);
    }
    
    unsigned c(C),d(D),e(E),f(F);
    if(!sol.between(A,C,E)) {
      std::swap(c,e);
      std::swap(d,f);
    }
    
    const tsp_prob& problem = tsp_prob::instance();

    int sub = problem.dist(A,B)+problem.dist(c,d)+problem.dist(e,f);

    which=0;
    int cost = problem.dist(A,d)+problem.dist(e,B)+problem.dist(c,f);
    int c2 = problem.dist(A,d)+problem.dist(e,c)+problem.dist(B,f);
    if (c2<cost) {
      which =1;
      cost = c2;
    }
    c2 = problem.dist(A,e)+problem.dist(B,d)+problem.dist(c,f);
    if (c2<cost) {
      which =2;
      cost = c2;
    }
    c2 = problem.dist(A,c)+problem.dist(B,e)+problem.dist(d,f);
    if (c2<cost) {
      which =3;
      cost = c2;
    }
    return cost-sub;
  }
  
  void operator()(tsp_prob::sol_t &sol) const {
    unsigned c(C),d(D),e(E),f(F);

#ifndef NDEBUG
    cost(sol);   // in debug mode, need to know which move to apply
#endif
    //    cout << " which: " << which << endl;

    if(!sol.between(A,C,E) && which!=4) {
      std::swap(c,e);
      std::swap(d,f);
      //      std::cout << "swaped" << endl;
    }
    switch (which) {   
    case 0:
      sol.flip(A,B,f,e);
      sol.flip(A,e,c,d);
      sol.flip(e,c,f,B);
      break;
    case 1:
      sol.flip(A,B,f,e);
      sol.flip(A,e,c,d);
      break;
    case 2:
      sol.flip(A,B,f,e);
      sol.flip(d,c,f,B);
      break;
    case 3:
      sol.flip(A,B,d,c);
      sol.flip(B,d,f,e);
      break;
    case 4:
      sol.flip(A,B,d,c);
    }
  }

private:
  unsigned A,B,C,D,E,F;
  mutable unsigned which;

  friend std::ostream& operator<<(std::ostream&x, const three_opt_move&m);
};


std::ostream& operator<<(std::ostream& x, const three_opt_move&m)
{
  x << m.A << " " << m.B << " " << m.C << " " << m.D << " " << m.E << " " << m.F << " " << m.which;
  return x;
}



struct three_opt_nh {
  three_opt_nh() 
    : dl_bits(tsp_prob::instance().size(), 0)
  { }

  template <class _oper>
  void operator()(_oper& op, const tsp_prob::sol_t& s) {
    const tsp_prob& problem = tsp_prob::instance();
    
    for (unsigned a=0; a<s.size(); ++a) {
      if (dl_bits[a]) {
	// 	continue;
      }
      
      int improve=0;
      const unsigned b = s.next(a);
      const std::vector<unsigned>& clist = problem.get_candidate_list(b);
      const unsigned dab=problem.dist(a,b);
      
      for(std::vector<unsigned>::const_iterator j=clist.begin();
	  j!=clist.end() && problem.dist(b,*j)<dab && improve==0;
	  ++j) {
	const unsigned& d = *j;
	const unsigned c = s.prev(d);
	
	const std::vector<unsigned>& elist = problem.get_candidate_list(d);
	const unsigned dab_cd = dab+problem.dist(c,d);
	const unsigned dbc = problem.dist(b,c);
	for(std::vector<unsigned>::const_iterator k=elist.begin();
	    k!=elist.end() && (dab_cd>dbc+problem.dist(d,*k));
	    ++k) {
	  const unsigned& f = *k;
	  const unsigned e = s.prev(f);

	  int same=0;
	  if (b==c) same++;
	  if (b==e) same++;
	  if (a==f) same++;
	  if (d==e) same++;
	  if (c==f) same++;
	  
	  if (same>1) continue; // don't consider degenerated moves that does nothing

	  if(op(three_opt_move(a,b,c,d,e,f))) {
	    improve=1;
	    dl_bits[a]=0;
	    dl_bits[b]=0;
	    dl_bits[c]=0;
	    dl_bits[d]=0;
	    dl_bits[e]=0;
	    dl_bits[f]=0;
	    break;
	  }
	}
      }
      if (improve==0) {
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
