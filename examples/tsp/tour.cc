#include "tour.hh"
#include "tsp_prob.hh"
#include "meta_utility.hh"
#include <iostream>

using namespace std;

tour::tour()
  : A(), B()
{
  // create a cannonical tour
  const unsigned size=tsp_prob::instance().size();
  A.reserve(size);
  B.reserve(size);
  for (unsigned i=0; i<size; ++i) {
    A.push_back(i);
  }
  updateB();
}

tour::tour(const std::vector<unsigned>& x) 
  : A(x), B()
{
  updateB();
}


void tour::updateB() 
{
  if(B.size()!=A.size()) {
    B = vector<unsigned>(A.size(),0);
  }

  for (unsigned i=0; i!=A.size(); ++i)
    B[A[i]] = i;
}

bool tour::checkB() const {
  if (B.size()!=A.size()) return false;
  for (unsigned i=0; i!=A.size(); ++i)
    if (B[A[i]]!=i) return false;
  return true;
}

void tour::flip(unsigned a, unsigned b, unsigned c, unsigned d)
{
  // remove edge ab, add edge ad
  // remove edge cd, add edge bc
  //  std::cout << "a: " << a << " b: " << b << " c: " << c <<" d: " <<d<< std::endl;
  //  print();
  //      std::cout << next(b) << " " << next(c) << std::endl;
  //      std::cout << next(a) << " " << next(d) << std::endl;
  if (next(b)==a) {
    std::swap(a,c);
    std::swap(b,d);
  }
  assert(next(a)==b);
  assert(next(d)==c);

  unsigned i=B[b];
  unsigned j=B[d];

  if (j<i) {
    j = B[a];
    i = B[c];
  }
  
  bool flipdir=0;

  

  if(flipdir==0) {
  //        std::cout << i << "   " << j << std::endl;
  //    if (j<i) std::swap(i,j);
    //    assert(i<j);

    while(i<j) {
      B[A[i]] = j;
      B[A[j]] = i;
      std::swap(A[i], A[j]);
	
      i++; j--;
    }
  } else {
    if (i==0) i=size();
    i--; 
    j++;
    if (j==size()) j=0;
     
    bool done=0;
    while (!done || i>j) {
      B[A[i]] = j;
      B[A[j]] = i;
      std::swap(A[i], A[j]);
      if (j==size()-1 && i==0) {
	break;  // we're done
      }
	
      if (j==size()-1) {
	i--;
	j=0;
	done=true;
	continue;
      }
      if (i==0) {
	i=size()-1;
	j++;
	done=true;
	continue;
      }
      i--; 
      j++;
    }
  }

// //   if (j<i) {
// //     j = B[a];
// //     i = B[c];
// //   }
//   //      std::cout << i << "   " << j << std::endl;
//   //  if (j-i < size()/2) {
//   if(1) {   /// FIXME !!!!!!!!!!
//     //    std::cout << "i " << i << "  j " << j << std::endl;
//     while(i<j) {
//       B[A[i]] = j;
//       B[A[j]] = i;
//       std::swap(A[i], A[j]);
	
//       i++; j--;
//     }
//   } else {
//     //      print();
//     // swap shortest subtour
//     //    std::cout << "save: " << i << " " << j << std::endl;
//     if (i==0) i=size();
//     i--; 
//     j++;
//     if (j==size()) j=0;
     
//     bool done=0;
//     while (!done || i>j) {
//       B[A[i]] = j;
//       B[A[j]] = i;
//       std::swap(A[i], A[j]);
//       if (j==size()-1 && i==0) {
// 	break;  // we're done
//       }
	
//       if (j==size()-1) {
// 	i--;
// 	j=0;
// 	done=true;
// 	continue;
//       }
//       if (i==0) {
// 	i=size()-1;
// 	j++;
// 	done=true;
// 	continue;
//       }
//       i--; 
//       j++;
//     }
//     //    print();
//   }
  assert(checkB());
}


// using four random points, do a "double bridge" move on the current
// solution.
// Cut the initial tour into 4 sub-tours: A B C D
// reconnect the subtours in the following order: D C B A
//  void tour::double_bridge()
//  {
//    int cut_points[5];
//    vector<unsigned> new_sol;
//    new_sol.reserve(size());
  
//     for (int i=0; i<4; i++)
//      cut_points[i]=metl::rng(size());

//    sort(cut_points, cut_points+4);
//    cut_points[4]=cut_points[0];

//    for (int i=3; i>=0; i--) 
//      for (int j=cut_points[i]; j!=cut_points[i+1]; j=(j+1)%(size())) {
//        new_sol.push_back(A[j]);
//      }

//    //  cout <<new_sol.size() << endl;
//    if (A.size()==new_sol.size()) {
//      A=new_sol;  // update the tour
//      updateB();
//    } else {
//      cout << "doh !"<< endl;
//    }
//  }



void tour::double_bridge()
{
  int cut_points[5];
  vector<unsigned> new_sol;
  new_sol.reserve(size());
  
  
  for (int i=0; i<4; i++)
    cut_points[i]=metl::rng(size()-1);
  
  sort(cut_points, cut_points+4);
  cut_points[4]=cut_points[0];


  for (int i=3; i>=0; i--) 
    for (int j=cut_points[i]; j!=cut_points[i+1]; j=(j+1)%size()) {
      new_sol.push_back(A[j]);
    }

  if (new_sol.size()==A.size()) {
    A=new_sol;  // update the tour
    updateB();
    //    print();
  }
}



void tour::print() {
  std::cout << "A: ";
  for (std::vector<unsigned>::iterator i=A.begin(); i!=A.end(); ++i)
    std::cout << *i << " ";
  std::cout << std::endl;
  std::cout << "B: ";
  for (std::vector<unsigned>::iterator i=B.begin(); i!=B.end(); ++i)
    std::cout << *i << " ";
  std::cout << std::endl;

}
