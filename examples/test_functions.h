#ifndef TEST_FUNCTIONS_H
#define TEST_FUNCTIONS_H

#include <sys/time.h>
#include <cmath>


double dt(const timeval& t1, const timeval& t0) {
  double ret=(t1.tv_usec - t0.tv_usec);
  ret/=1000000;
  return ret+t1.tv_sec-t0.tv_sec;
}

template <class _prob>
void test_generator(metl::metaheuristic<_prob>* mh) {
  timeval timestart;
  timeval timeend;
  
  gettimeofday(&timestart,0);
  typename _prob::soleval_t se = (*mh)();
  gettimeofday(&timeend,0);
#ifdef USE_MPI
  if (MPI::COMM_WORLD.Get_rank()==0)
#endif
    {
      typename _prob::eval_t x_eval = _prob::instance().evaluation(se.first);
      
      if (fabs(double(x_eval - se.second)) > 0.1) {
	std::cout << "WARNING: metaheuristic returned : "<< se.second<< "  evaluation is: " << x_eval << std::endl;
      }

      std::cout << mh->name() << ": "<<  
	se.second << " " << 
	_prob::instance().is_valid(se.first) 
		<< " time: " << dt(timeend, timestart)<<std::endl; 
    }
}



template <class _prob>
void test_metaheuristic(metl::metaheuristic<_prob>* mh, const typename _prob::soleval_t& solution) {
  timeval timestart;
  timeval timeend;
  
  gettimeofday(&timestart,0);
  typename _prob::soleval_t se = (*mh)(solution);
  gettimeofday(&timeend,0);
#ifdef USE_MPI
  if (MPI::COMM_WORLD.Get_rank()==0)
#endif
    {
      typename _prob::eval_t x_eval = _prob::instance().evaluation(se.first);
      
      if (fabs(double(x_eval - se.second)) > 0.1) {
	std::cout << "WARNING: metaheuristic returned : "<< se.second<< "  evaluation is: " << x_eval << std::endl;
      }

      std::cout << mh->name() << ": "<<  
	se.second << " " << _prob::instance().is_valid(se.first) << " time: " << dt(timeend, timestart)<<std::endl; 
    }
}


#endif
