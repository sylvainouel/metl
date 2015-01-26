#ifndef MPI_REDUCE_OP_HH
#define MPI_REDUCE_OP_HH
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


#include "metl_def.hh"

#ifdef _OPENMP
#ifndef INCLUDED_OMP_H   
#define INCLUDED_OMP_H   
#include <omp.h>         
#endif                   
#else                    
#include "omp_stub.h"    
#endif


#ifdef USE_MPI

// this is a workaround for a bug in BOOST 1_32. 
// If your compiler does not support exceptions or you want to compile without exceptions, then define BOOST_NO_EXCEPTIONS.
#ifdef BOOST_NO_EXCEPTIONS
#undef BOOST_NO_EXCEPTIONS
#include <boost/archive/archive_exception.hpp>
#define BOOST_NO_EXCEPTIONS
#endif

#include <sstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include "metl_config.hh"

namespace metl {


  template <class T>
  void pair_reduce(const void *invec, void *inoutvec, int len, const MPI::Datatype& datatype) {
    const T* in = static_cast<const T*>(invec);
    T* inout = static_cast<T*>(inoutvec);
  
    //  cerr << "len: " << len << endl;
    for (int i=0; i<len; ++i) {
      //    cerr << "reduction: " << in[i].second << " " << inout[i].second << endl;
      if (in[i].second < inout[i].second) 
	inout[i] = in[i];
    }
  }


  template <class sol_t, class eval_t>
  void pair_reduce_serialized(const void *invec, void *inoutvec, int len, const MPI::Datatype& datatype) {
    const unsigned size = datatype.Get_size();

    const char* in = static_cast<const char*>(invec);
    char* inout = static_cast<char*>(inoutvec);
  
    sol_t sol1, sol2;
    eval_t eval1, eval2;

    //  cerr << "len: " << len << endl;
    for (int i=0; i<len; ++i) {
      std::istringstream in_buf(std::string(in,size));
      boost::archive::binary_iarchive ia(in_buf);
      ia >> sol1 >> eval1;

      std::istringstream inout_buf(std::string(inout,size));
      boost::archive::binary_iarchive ioa(inout_buf);
      ioa >> sol2 >> eval2;

      //    cerr << "reduction: " << in[i].second << " " << inout[i].second << endl;
      if (eval1 < eval2) {
	memcpy(inout, in, size);
      }
      in+=size;
      inout+=size;
    }
  }


  template <class prob_t>
  struct sol_reducter {    // this is a singleton

    static void Allreduce(typename prob_t::soleval_t& se) {
      static sol_reducter<prob_t> instance(se);

      std::ostringstream send_buf(std::ios::binary);
      // serialize solution into buffer
      boost::archive::binary_oarchive oa(send_buf);
      const typename prob_t::soleval_t& cse=se;
      oa << cse.first << cse.second;
    
      const std::string sbuf = send_buf.str();
    
      if (sbuf.size() > max_buf_size) {
	std::cerr << "max_buf_size overflow" << std::endl;
	abort();	
      }

      if (sbuf.size() != instance.serialized_size) {
	// REQUIS: serialized size must be the same !
	std::cerr << "Serialized size must be constant" << std::endl;
	abort();
      }
      
      MPI::COMM_WORLD.Allreduce(sbuf.data(), instance.buf, 1, instance.type_pair, instance.pair_reduce_op);
      // unserialize solution
      std::istringstream receive_buf(std::string(instance.buf,instance.serialized_size));
      boost::archive::binary_iarchive ia(receive_buf);
      ia >> se.first >> se.second;
    }


  private:

    sol_reducter(const typename prob_t::soleval_t& se) 
      : pair_reduce_op(),
	type_pair(),
	serialized_size(0)
    {
    
      pair_reduce_op.Init(pair_reduce_serialized<typename prob_t::sol_t, typename prob_t::eval_t>,0);
      
      std::ostringstream send_buf(std::ios::binary);
      // serialize solution into buffer
      boost::archive::binary_oarchive oa(send_buf);
      
      oa << se.first << se.second;
      
      serialized_size = send_buf.str().size();
      //      std::cout << "serialized solution size: " << serialized_size << std::endl;
      type_pair = MPI::CHAR.Create_contiguous(serialized_size);
      type_pair.Commit();
    }

    ~sol_reducter() {
      pair_reduce_op.Free();
      type_pair.Free();
    }

    sol_reducter(const sol_reducter&);    // forbid copy construction
    sol_reducter& operator=(sol_reducter); // forbid assignment

    char buf[max_buf_size];
    
    MPI::Op pair_reduce_op;  
    MPI::Datatype type_pair;
    unsigned serialized_size;
  };


}

#endif

#endif
