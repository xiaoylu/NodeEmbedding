#ifndef MATRIX_HEADER
#define MATRIX_HEADER 

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <set>
#include <atomic>

#include "cascade.h"
#include "vec.h"
#include "math.h"

using namespace std;


/* The global state of all node vectors
   A[i] is the influence vector of node i
   B[i] is the selectivity vector of node i
   Both matrices sit in the shared memory
  */
class Matrix2
{
private:
  atomic<double> **A; // infection matrix  
  atomic<double> **B; // suspection matrix  

  //vector<vector<atomic<double>>> A; // lockless
  //vector<vector<atomic<double>>> B; // lockless

  int n; // number of nodes
  int m; // number of topics
  Vec * Bsum; // the sum of B_v for all nodes v in the network  

  // private init function
  Matrix2() { };

  // private destory function, release memory until program exits
  ~Matrix2(); 

  // private copy function 
  //Matrix2(Matrix2 const&);
  //void operator=(Matrix2 const&);

public:
  // C++11, disallow copy
  Matrix2(Matrix2 const&) = delete;
  void operator=(Matrix2 const&) = delete;

  friend class Cascade; // so that Cascade can access the A,B matrices
                        // during training process

  // Singleton pattern, return a reference
  //   to prevent the delete operation of pointers
  static Matrix2 &instance()
  {
      static Matrix2 s_instance;
      return s_instance; 
  }

  // Allocate memory and initlization A, B
  void init(int n, int m);
  
  // compute the summation of all B rows
  void sumB();

  int getN() { return n; };

  // output the A,B matrices to txt files
  void writeMat(string path);

  //==================utility function=================
  void getMaxNorm(vector<pair<int,double> >& ret);
  void getMostInfluencialSites(vector<pair<int,double> >& ret, int k);
  void printMaxComponent()
  {
    for (int i = 0;i < n;++i)
    {
      if (i > 10 && i < n-30)
      {
         cout << "..." << endl;
         i = n-30;
      }
      cout << i << " : " << "A" << Vec(m, A[i]).amax(); 
      cout << " B" << Vec(m, B[i]).amax() << endl;
    }
  }
  
};

#endif
