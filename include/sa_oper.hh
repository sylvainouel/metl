#ifndef SA_OPER_HH
#define SA_OPER_HH

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

template <class prob_t, class _move>
struct sa_accept_scheme {
  virtual ~sa_accept_scheme() {}
  virtual bool operator()(const _move m)=0;

  void reset_counters() {
    reject=0;
  }

  unsigned get_rejects() const { return reject; }

  // FIXME: should realy be a friend function because it does not make sense to have it in the public interface
  void set_sol_and_eval(typename prob_t::sol_t* so , typename prob_t::eval_t* e)
  {
    s = so;
    _sol_eval = e;
  }


protected:
  sa_accept_scheme(double &temp)
    : s(0),
      T(temp), 
      _sol_eval(0),
      reject(0)
  {}

  
  inline bool move_accept(const _move& m, const typename prob_t::eval_t& e) {
    m(*s);
    *_sol_eval += e;
    reject=0;
    return true;
  }

  inline bool move_reject() {
    ++reject;
    return false;
  }

  
  sa_accept_scheme& operator=(const sa_accept_scheme& r)
  {
    if (this == &r) return *this;
    s=r.s;
    // no assign T
    _sol_eval =r._sol_eval;
    reject = r.reject;

    return *this;
  }

protected:
  //  typename prob_t::sol_t& s; 
  typename prob_t::sol_t* s;
  double& T;

private:
  //  typename prob_t::eval_t& _sol_eval;
  typename prob_t::eval_t* _sol_eval;
  unsigned reject;

  sa_accept_scheme(const sa_accept_scheme<prob_t, _move>&);
};



template <class prob_t, class _move>
struct metropolis: public sa_accept_scheme<prob_t, _move> {
  metropolis(double &temp) 
    : sa_accept_scheme<prob_t,_move>(temp)
  {}


  inline bool operator()(const _move m)
  {
    const typename prob_t::eval_t e=m.internal_cost(*(this->s));

    if (e<=0 || rng() < exp(double(-e)/this->T)) {
      return move_accept(m,e);
    }
    return this->move_reject();
  }
};


template <class prob_t, class _move>
struct threshold_accept: public sa_accept_scheme<prob_t, _move> {
  threshold_accept(double &temp) 
    : sa_accept_scheme<prob_t,_move>(temp),
      _threshold(0)
  {}

  void set_threshold(typename prob_t::eval_t threshold)
  {
    _threshold = threshold;
  }

  inline bool operator()(const _move m)
  {
    const typename prob_t::eval_t e=m.internal_cost(*(this->s));

    if (e<_threshold) {
      move_accept(m,e);
      return true;
    }
    this->move_reject();
    return false;
  }
private:
  typename prob_t::eval_t _threshold;
};



template <class _neighborhood>
struct cooling_scheme_base {

  void set_neighborhood(_neighborhood& _nh)
  {
    // the cooling scheme is specialized for the neighborhood so that
    // it is possible to implement a specialized cooling scheme that
    // "act" on the neighborhood when changing the temperature.
    nh = &_nh;
  }

  cooling_scheme_base& operator=(const cooling_scheme_base& r)
  {
    if (this == &r) return *this;
    // NO assignment to T. (anyway it's const to protect it)
    it = r.it;
    return *this;
  }

protected:
  cooling_scheme_base(double* temp)
    : T(temp),
      it(0),
      nh(0)
  {}

  virtual ~cooling_scheme_base() {}
  // returns true if temperature change
  virtual bool operator()() {
    ++it;
    return false;
  }


  double* const T;  // cannot change where T points
  unsigned it;
  _neighborhood* nh;

private:
  cooling_scheme_base(const cooling_scheme_base&); // forbid copy construction
    
};



// this is a geometric cooling scheme with steps
template <class _neighborhood>
struct cooling_geometric_steps : public cooling_scheme_base<_neighborhood> {
  cooling_geometric_steps(double* temp)
    : cooling_scheme_base<_neighborhood>(temp),
      cf(0.05),   // useless default value
      sl(1),     // useless default value
      next_step(sl)
  {}

  void set_cooling_factor(double cooling_factor) { cf=cooling_factor; }
  void set_step_length(unsigned step_length) { 
    sl = step_length; 
    next_step = step_length;
  }
  
  
  bool operator()() {    // returns true if temperature was modified
    cooling_scheme_base<_neighborhood>::operator()();
    if (this->it >= next_step)
      {
	next_step+=sl;
	(*(this->T)) *= (1-cf);   // modify the temperature
	return true;
      }
    return false;
  }

  
protected:
  double cf;
  unsigned sl;
  unsigned next_step;
};



}

#endif
