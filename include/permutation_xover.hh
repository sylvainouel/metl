#ifndef PERMUTATION_XOVER
#define PERMUTATION_XOVER


/*
metl: A generic framework for sequential and parallel metaheuristics
Copyright (c) 2005-2015, Sylvain Ouellet


Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

*/

namespace metl {

// implementation of path relinking (path crossover)
// this is usefull for problems in which a solution is a
// permutation of elements.
// It is an heuristic crossover that is a least respectfull
template <class prob_t>
struct path_xover: public abstract_crossover<prob_t> {
  
  // alpha: probability of not selecting best child
  void set_alpha(float _alpha) { alpha = _alpha; }


protected:
  path_xover()
    : minimum(std::numeric_limits<unsigned>::max()),
      alpha(0)
  {}

  template <class op>
  void operator()(const typename prob_t::soleval_t& parent1, 
		  const typename prob_t::soleval_t& parent2, 
		  typename prob_t::soleval_t& child, 
		  op& my_op) const 
  {
    // only works with random access containers

    if (parent1==parent2) {
      child = parent1;
      return;
    }

    bool nc = true;

    // make copies of the two parents
    typename prob_t::soleval_t p1 = parent1;
    typename prob_t::soleval_t p2 = parent2;

    // choose a random start position
    const unsigned start = rng(parent1.first.size());     
    typename prob_t::sol_t::iterator a = p1.first.begin()+start;
    typename prob_t::sol_t::iterator b = p2.first.begin()+start;
    const typename prob_t::sol_t::const_iterator end = a;

    const prob_t& instance = prob_t::instance();

    do {
      if (*a!=*b) {
	my_op(p1,p2,a,b);

	if (p1.second <= p2.second) {
	  if (nc || (p1.second < child.second && rng()>alpha)) {
	    child = p1;
	    nc=false;
	  }
	} else {
	  if (nc || (p2.second < child.second && rng()>alpha)) {
	    child = p2;
	    nc=false;
	  }
	}
      } 

      ++a;
      ++b;
      if (a==p1.first.end()) {
	a=p1.first.begin();
	b=p2.first.begin();
      }
    } while (end!=a);

  
    return;
  }

private:
  unsigned minimum;

protected:
  float alpha;
};



template <class prob_t>
struct op_path_xover {
  op_path_xover(const typename prob_t::soleval_t& parent1, 
		const typename prob_t::soleval_t& parent2)
    : minimum(*std::min_element(parent1.first.begin(), parent1.first.end())),
      map1(parent1.first.size()),
      map2(parent2.first.size())
  {
    // compute the mapx (give position of element i in each parent
    for (unsigned i=0;i<parent1.first.size(); ++i) {
      map1[parent1.first[i]-minimum]=i; 
      map2[parent2.first[i]-minimum]=i; 
    }
  }

  virtual ~op_path_xover() {}

  
  virtual void operator()(typename prob_t::soleval_t& p1,
			  typename prob_t::soleval_t& p2,
			  const typename prob_t::sol_t::iterator& a,
			  const typename prob_t::sol_t::iterator& b)=0;

protected:
  unsigned minimum;
  std::vector<unsigned> map1;
  std::vector<unsigned> map2;
};



  // this is the internal operation for the swap version of path xover
template <class prob_t, class permut_move=permutation_move<prob_t> >
struct op_path_xover_swap: public op_path_xover<prob_t> {
  typedef op_path_xover<prob_t> base;

  op_path_xover_swap(const typename prob_t::soleval_t& parent1, 
		     const typename prob_t::soleval_t& parent2)
    : base(parent1,parent2)
  {}

  void operator()(typename prob_t::soleval_t& p1,
		  typename prob_t::soleval_t& p2,
		  const typename prob_t::sol_t::iterator& a,
		  const typename prob_t::sol_t::iterator& b)
  {
    const typename prob_t::sol_t::value_type ta = *a-this->minimum;
    const typename prob_t::sol_t::value_type tb = *b-this->minimum;
    const permut_move m1(this->map1[ta], this->map1[tb]);
    const permut_move m2(this->map2[ta], this->map2[tb]);

    p1.second +=  m1.internal_cost(p1.first);
    p2.second +=  m2.internal_cost(p2.first);    
    m1(p1.first);
    m2(p2.first);
    
    std::swap(this->map1[ta],this->map1[tb]);
    std::swap(this->map2[ta],this->map2[tb]);
  }
};


  // this is a partial specialization in case the user did not define
  // a specialized permutation move. It is just a bit faster.
template <class prob_t>
struct op_path_xover_swap<prob_t, permutation_move<prob_t> >: public op_path_xover<prob_t> {
  typedef op_path_xover<prob_t> base;

  op_path_xover_swap(const typename prob_t::soleval_t& parent1, 
		     const typename prob_t::soleval_t& parent2)
    : base(parent1,parent2)
  {}

  void operator()(typename prob_t::soleval_t& p1,
		  typename prob_t::soleval_t& p2,
		  const typename prob_t::sol_t::iterator& a,
		  const typename prob_t::sol_t::iterator& b)
  {
    const typename prob_t::sol_t::value_type ta = *a-this->minimum;
    const typename prob_t::sol_t::value_type tb = *b-this->minimum;

    std::swap(p1.first[ta], p1.first[tb]);
    std::swap(p2.first[ta], p2.first[tb]);
    std::swap(this->map1[ta],this->map1[tb]);
    std::swap(this->map2[ta],this->map2[tb]);

    const prob_t& instance = prob_t::instance();
    p1.second = instance.evaluation(p1.first);
    p2.second = instance.evaluation(p2.first);

  }
};



  // this is the internal operation for insert path xover
template <class prob_t>
struct op_path_xover_insert: public op_path_xover<prob_t> {
  typedef op_path_xover<prob_t> base;

  op_path_xover_insert(const typename prob_t::soleval_t& parent1, 
		       const typename prob_t::soleval_t& parent2)
    : base(parent1,parent2)
  {}

  void operator()(typename prob_t::soleval_t& p1,
		  typename prob_t::soleval_t& p2,
		  const typename prob_t::sol_t::iterator& a,
		  const typename prob_t::sol_t::iterator& b)
  {
    const typename prob_t::sol_t::value_type& ta = *a-this->minimum;
    const typename prob_t::sol_t::value_type& tb = *b-this->minimum;


    // perform move in p1

    typename prob_t::sol_t tp1 = p1.first;
    do_insert(tp1, p2.first, p1.first, 1, tb, a - p1.first.begin());
    do_insert(tp1, p2.first, p2.first, 2, ta, b - p2.first.begin());

    const prob_t& instance = prob_t::instance();
    p1.second = instance.evaluation(p1.first);
    p2.second = instance.evaluation(p2.first);
  }

  private:
    
    void do_insert(const typename prob_t::sol_t& tp1, 
		   const typename prob_t::sol_t& p2,
		   typename prob_t::sol_t& pa, 
		   unsigned which,
		   const typename prob_t::sol_t::value_type& t,
		   unsigned start)
  {
    std::vector<unsigned>& mapx = which==1 ? this->map1 : this->map2;
    
    unsigned hole = mapx[t];
    unsigned pos = hole;
    unsigned pos1 = hole;
    

    // on decale le contenue du vecteur pa pour avoir de la place pour
    // faire l'insertion.

    // on copie de pos1 vers pos

    do {
      // avance pos1 pour trouver la position de lecture
      // on saute chaque fois que les 2 elemets sont egaux
      do {
	if (pos1==0) pos1=pa.size();
	--pos1;
	if (pos1 == start) break;
      } while (tp1[pos1] == p2[pos1]);

      pa[pos] = pa[pos1];   // copie
      mapx[pa[pos]] = pos;  // update map
      pos = pos1;   // avance pos

    } while (pos1!=start);
    
    pa[start] = t;   // effectue l'insertion
    mapx[t] = start;
  }

};




// The user should (and probably does) provide a specialized move
// function (that is a permutation_move) for is problem. In such a
// case, path_xover_swap will be much faster.
template <class prob_t,class permut_move=permutation_move<prob_t> >
struct path_xover_swap: public path_xover<prob_t> { 
  typedef path_xover<prob_t> base;

  void operator()(const typename prob_t::soleval_t& parent1, 
		  const typename prob_t::soleval_t& parent2, 
		  typename prob_t::soleval_t& child) const {

    op_path_xover_swap<prob_t, permut_move> op(parent1, parent2);

    base::operator()(parent1, parent2, child, op);
  }
};


  // this is the insertion version of path relinking
template <class prob_t>
struct path_xover_insert: public path_xover<prob_t> { 
  typedef path_xover<prob_t> base;

  void operator()(const typename prob_t::soleval_t& parent1, 
		  const typename prob_t::soleval_t& parent2, 
		  typename prob_t::soleval_t& child) const {

    op_path_xover_insert<prob_t> op(parent1, parent2);

    base::operator()(parent1, parent2, child, op);
  }
};

}

#endif
