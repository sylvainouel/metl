#ifndef EXCHANGE_OPER_HH
#define EXCHANGE_OPER_HH

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


#include <iostream>
#include <list>

// this is a workaround for a bug in BOOST 1_32. 
// If your compiler does not support exceptions or you want to compile without exceptions, then define BOOST_NO_EXCEPTIONS.
#ifdef BOOST_NO_EXCEPTIONS
#undef BOOST_NO_EXCEPTIONS
#include <boost/archive/archive_exception.hpp>
#define BOOST_NO_EXCEPTIONS
#endif

#include <boost/shared_array.hpp>

#include "mpi_reduce_op.hh"

#include "metl_def.hh"

#ifdef USE_MPI
#ifdef HAVE_MPIPP_H     
#include <mpi++.h>
#else
#include <mpi.h>
#endif  // HAVE_MPIPP_H
#endif  // USE_MPI

#ifdef _OPENMP
#ifndef INCLUDED_OMP_H   
#define INCLUDED_OMP_H   
#include <omp.h>         
#endif                   
#else                    
#include "omp_stub.h"    
#endif


#ifdef USE_MPI
#include <sstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include "metl_config.hh"

namespace metl {


const int tag_done = 0;
const int tag_xchange = 1;
const int tag_xchange_reply = 2;
}

#endif  // USE_MPI

namespace metl {

enum AsyncExchangePolicy {BEST, RANDOM, RANDOM_BETTER};


template <class individu>
bool individus_compare(const individu& a, const individu& b) {
  // only compare on the eval
  return a.second<b.second;
}



////////////// structure qui permet the conserver les x meilleurs solutions
template <class prob_t>
class best_pool {
public:  

  best_pool(unsigned pool_size) 
    : max_size(pool_size),
      pool()
  {
    //    pool.reserve(max_size+1);
  }
    
  void insert(const typename prob_t::soleval_t& se) {
#pragma omp critical (OMP_BEST_POOL_ACCESS)
    {
      if (se.second < pool.back().second || pool.size()<max_size) {
	pool.insert(lower_bound(pool.begin(), pool.end(), se, individus_compare<typename prob_t::soleval_t>), se);
	if (pool.size() > max_size)
	  pool.pop_back();
      }
    }
  }

  void get_best(typename prob_t::soleval_t& se) const {
#pragma omp critical (OMP_BEST_POOL_ACCESS)
      {
	se = pool.front();
      }
  }

  void get_random(typename prob_t::soleval_t& se) const
  {
#pragma omp critical (OMP_BEST_POOL_ACCESS)
    {
      typename std::list<typename prob_t::soleval_t>::const_iterator x=pool.begin();
      advance(x,rng(pool.size()));
      se = *x;
    }
  }

  void get_random_better(typename prob_t::soleval_t& se) const
  {
#pragma omp critical (OMP_BEST_POOL_ACCESS)
    {
      typename std::list<typename prob_t::soleval_t>::const_iterator x=pool.begin();
      typename std::list<typename prob_t::soleval_t>::const_iterator ep=lower_bound(pool.begin(), pool.end(), se, individus_compare<typename prob_t::soleval_t>);
      advance(x,rng(distance(pool.begin(),ep)));
      se = *x;
    }
  }

private:
  unsigned max_size;
  std::list<typename prob_t::soleval_t> pool;  // FIXME: est-ce la meilleure structure
			     // ??  Une liste est quand meme pas
			     // si mal tant que la structure reste
			     // petite
};


template <class prob_t>
class exchange_op {
public:
  
  // does both send and receive
  virtual bool operator()(typename prob_t::soleval_t& se)
  {
    send(se);
    return recv(se);
  }

  virtual void send(const typename prob_t::soleval_t& se)=0;
  // returns true if sol is modified false otherwise
  virtual bool recv(typename prob_t::soleval_t& se)=0;

protected:
  exchange_op(unsigned exchange_period, bool _eval_recv=false)
    : period(exchange_period), cycle(0), eval_recv(_eval_recv)
  {  }

  virtual ~exchange_op() {}

  const unsigned period;
  unsigned cycle;
  bool eval_recv;
};



template <class prob_t>
struct blackboard_coop_op: public exchange_op<prob_t> {

protected:
  blackboard_coop_op(unsigned exchange_period, AsyncExchangePolicy exchange_policy, bool _eval_recv=false)
    : exchange_op<prob_t>(exchange_period, _eval_recv),
      policy(exchange_policy)
  { }

  void do_exchange(best_pool<prob_t>& pool, typename prob_t::soleval_t& se) const
  {
    switch (policy) {
    case RANDOM:
      pool.get_random(se);
      break;
    case BEST:
      pool.get_best(se);
      break;
    case RANDOM_BETTER:
      // retourne une solution aleatoirement parmis celles qui sont meilleures que la solution courante.
      pool.get_random_better(se);
      break;
    }
  }

private:
  const AsyncExchangePolicy policy;
};


// Une operation d'échange asynchrone avec OpenMP
template <class prob_t>
class omp_blackboard_coop_op: public blackboard_coop_op<prob_t> {
  typedef blackboard_coop_op<prob_t> base;
public:
  omp_blackboard_coop_op(best_pool<prob_t>& shared_pool, unsigned exchange_period, AsyncExchangePolicy exchange_policy=RANDOM, bool _eval_recv=false)
    : base(exchange_period, exchange_policy, _eval_recv),
      pool(shared_pool)
  {}

  void send(const typename prob_t::soleval_t& se)
  {
    if (++(this->cycle) == this->period) {
      pool.insert(se);
    }
  }

  bool recv(typename prob_t::soleval_t& se)
    {
      if (this->cycle == this->period) {
	this->cycle=0;

	base::do_exchange(pool,se);
	if (this->eval_recv)
	  se.second = prob_t::instance().evaluation(se.first);
	return true;
      }
      return false;
    }

private:
  best_pool<prob_t>& pool;
};




#ifdef USE_MPI

template <class prob_t>
class mpi_blackboard_coop_op : public blackboard_coop_op<prob_t> {
  typedef blackboard_coop_op<prob_t> base;
public:
  mpi_blackboard_coop_op(unsigned exchange_period, AsyncExchangePolicy exchange_policy=RANDOM, bool _eval_recv=false)
    : base(exchange_period, exchange_policy, _eval_recv)
  {}


  void tell_master(const typename prob_t::soleval_t& se) const {
    // workers must tell the master they are done
    // Once they finish, they send their best solution with it's
    // evaluation to the master
    if (MPI::COMM_WORLD.Get_rank()!=1) {

      std::ostringstream send_buf(std::ios::binary);
      // serialize solution into buffer
      boost::archive::binary_oarchive oa(send_buf);  // create output archive
      oa << se.first << se.second;
      
      // send it
      const std::string sbuf = send_buf.str();
      assert(sbuf.size() < max_buf_size);
      MPI::COMM_WORLD.Send(sbuf.data(), sbuf.size(), MPI::CHAR, 1, tag_done);
    }
  }

  void master(best_pool<prob_t>& pool) {
    unsigned workers = MPI::COMM_WORLD.Get_size()-1;
    char rbuf[max_buf_size];
    MPI::Status rstatus;
    typename prob_t::soleval_t se;

    while(workers>0) {
      MPI::COMM_WORLD.Recv(rbuf, max_buf_size, MPI::CHAR, MPI::ANY_SOURCE, MPI::ANY_TAG, rstatus);

      // un-serialize buffer
      std::istringstream receive_buf(std::string(rbuf,rstatus.Get_count(MPI::CHAR)));
      boost::archive::binary_iarchive ia(receive_buf);
      ia >> se.first >> se.second;
      // add item into best_pool
      pool.insert(se);


      switch (rstatus.Get_tag()) {
      case tag_done:
	--workers;
	break;
      case tag_xchange:
	//	std::cout << "worker : " << rstatus.Get_source() << " insert: " << se.second;
	// retreive best sol
	base::do_exchange(pool,se);
	//	std::cout << " retreive: " << se.second << std::endl;

	// serialize new solution
	std::ostringstream send_buf(std::ios::binary);
	// serialize solution into buffer
	boost::archive::binary_oarchive oa(send_buf);  // create output archive
	const typename prob_t::soleval_t& cse = se;
	oa << cse.first << cse.second;

	// send it
	const std::string sbuf = send_buf.str();
	assert(sbuf.size() < max_buf_size);
	MPI::COMM_WORLD.Send(sbuf.data(), sbuf.size(), MPI::CHAR, rstatus.Get_source(), tag_xchange_reply);
	break;
      }
    }
  }
 

  void send(const typename prob_t::soleval_t& se)
  {
    if (++(this->cycle) == this->period) {
      std::ostringstream send_buf(std::ios::binary);
      // serialize solution into buffer
      boost::archive::binary_oarchive oa(send_buf);
      oa << se.first << se.second;
      
      const std::string sbuf = send_buf.str(); // string is reference-counted, so copy is cheap.
      assert(sbuf.size() < max_buf_size);
      
      MPI::COMM_WORLD.Isend(sbuf.data(), sbuf.size(), MPI::CHAR, 1, tag_xchange);
    }
  }
  
  bool recv(typename prob_t::soleval_t& se)
  {
    if (this->cycle == this->period) {
      this->cycle=0;

      MPI::Status rstatus;
      MPI::COMM_WORLD.Recv(rbuf, max_buf_size, MPI::CHAR, 1, tag_xchange_reply, rstatus);
      // unserialize solution
      std::istringstream receive_buf(std::string(rbuf,rstatus.Get_count(MPI::CHAR)));
      boost::archive::binary_iarchive ia(receive_buf);
      ia >> se.first >> se.second;

      if (this->eval_recv)
	se.second = prob_t::instance().evaluation(se.first);

      CHECK_EVAL(se);

      return true;
    }
    return false;
  }

private:
  char rbuf[max_buf_size];
};


template <class prob_t>
struct mpi_ring_coop_op : public exchange_op<prob_t> {
  mpi_ring_coop_op(unsigned exchange_period, bool _eval_recv=false) 
    : exchange_op<prob_t>(exchange_period, _eval_recv),
      rank(MPI::COMM_WORLD.Get_rank()),
      size(MPI::COMM_WORLD.Get_size())
  {  }

  void send(const typename prob_t::soleval_t& se)
  {
    if (++(this->cycle) == this->period) {
      const unsigned send_to = rank==size-1 ? 0 : rank+1;

      std::ostringstream send_buf(std::ios::binary);
      // serialize solution into buffer
      boost::archive::binary_oarchive oa(send_buf);
      oa << se.first << se.second;
      
      const std::string sbuf = send_buf.str();
      assert(sbuf.size() < max_buf_size);
      
      MPI::COMM_WORLD.Isend(sbuf.data(), sbuf.size(), MPI::CHAR, send_to, tag_xchange);
    }
  }


  bool recv(typename prob_t::soleval_t& se)
  {
    if (this->cycle == this->period) {
      this->cycle=0;
      
      const unsigned recv_from = rank==0 ? size-1 : rank-1;
      
      MPI::Status rstatus;
      MPI::COMM_WORLD.Recv(rbuf, max_buf_size, MPI::CHAR, recv_from, tag_xchange,rstatus);
      // unserialize solution
      std::istringstream receive_buf(std::string(rbuf,rstatus.Get_count(MPI::CHAR)));
      boost::archive::binary_iarchive ia(receive_buf);
      ia >> se.first >> se.second;
      
      if (this->eval_recv)
	se.second = prob_t::instance().evaluation(se.first);

      return true;  // we modified the solution
    }
    return false;
  }

private:
  char rbuf[max_buf_size];
  const unsigned rank;
  const unsigned size;
  
};



template <class prob_t>
struct mpi_reduce_coop_op : public exchange_op<prob_t> {
  mpi_reduce_coop_op(unsigned exchange_period, bool _eval_recv=false) 
    : exchange_op<prob_t>(exchange_period, _eval_recv)
  {}

  void send(const typename prob_t::soleval_t& se)
  {}


  // recv actual does both send and receive
  bool recv(typename prob_t::soleval_t& se)
  {
    if (++(this->cycle) == this->period) {
      this->cycle=0;
      
      sol_reducter<prob_t>::Allreduce(se);
      
      if (this->eval_recv)
	se.second = prob_t::instance().evaluation(se.first);

      return true;  // we modified the solution
    }
    return false;
  }

};


#endif   // USE_MPI



template <class prob_t>
struct omp_ring_coop_op : public exchange_op<prob_t> {
  omp_ring_coop_op(unsigned exchange_period, boost::shared_array<typename prob_t::soleval_t>& shared_vect, bool _eval_recv=false) 
    : exchange_op<prob_t>(exchange_period, _eval_recv),
      rank(omp_get_thread_num()),
      size(omp_get_max_threads()),
      _shared(shared_vect)
  {  }

  void send(const typename prob_t::soleval_t& se)
  {
    if (++(this->cycle) == this->period) {
      const unsigned send_to = rank==size-1 ? 0 : rank+1;

      _shared[send_to] = se;
    }
  }


  bool recv(typename prob_t::soleval_t& se)
  {
    if (this->cycle == this->period) {
      this->cycle=0;
      const unsigned recv_from = rank==0 ? size-1 : rank-1;
#pragma omp barrier
      se = _shared[recv_from];
#pragma omp barrier
      if (this->eval_recv)
	se.second = prob_t::instance().evaluation(se.first);

      return true;  // we modified the solution
    }
    return false;
  }

private:
  const unsigned rank;
  const unsigned size;
  boost::shared_array<typename prob_t::soleval_t>& _shared;
};



template <class prob_t>
struct omp_reduce_coop_op : public exchange_op<prob_t> {
  omp_reduce_coop_op(unsigned exchange_period, typename prob_t::soleval_t* shared_individu, bool _eval_recv=false) 
    : exchange_op<prob_t>(exchange_period, _eval_recv),
      _shared(shared_individu)
  {  }

  void send(const typename prob_t::soleval_t& se)
  {
    if (++(this->cycle) == this->period) {
#pragma omp critical (OMP_ACCESS_SHARED_SOL)
      {
	if (se.second < _shared->second) {
	  *_shared = se;
	}
      }
    }
  }


  bool recv(typename prob_t::soleval_t& se)
  {
    if (this->cycle == this->period) {
      this->cycle=0;
      
      // wait to make sure everyone has done the send()
#pragma omp barrier

      se = *_shared;

      // wait to make sure everyone has finished doing the recv()
#pragma omp barrier
#pragma omp single nowait
      {
	_shared->second = std::numeric_limits<typename prob_t::eval_t>::max();
      }

      if (this->eval_recv)
	se.second = prob_t::instance().evaluation(se.first);

      return true;  // we modified the solution
    }
    return false;
  }

private:
  typename prob_t::soleval_t* _shared;
};


}

#endif
