#ifndef cascade_header
#define cascade_header

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <float.h>
#include <tuple>

#include "vec.h"
#include "matrix.h"

using namespace std;

class Cascade
{
private:
  // hash : cascade ID --->  node id_1, infection time
  //                   --->  node id_2, infection time
  //                   --->  ...
  //                   --->  node id_n, infection time
  vector<pair<int, double> > * casc;

  // the string ID of this cascade
  double dimension;
  double OW;
  string strID;
  double stepsize;
  double epsilon;  
  int negative_sampling_size;

  Vec **H, **G, **P, **Q;
  // H[i] = \sum_i from 0 to i-1 A[node_i]
  // P[i] = \sum_i from i+1 to size-1 B[node_i]
  // Note: H[size] = \sum_i for all i A[node_i]
  //       P[size] = \sum_i for all i B[node_i]

public:
  Cascade(string strID_)
       : strID(strID_), H(nullptr), G(nullptr), P(nullptr), Q(nullptr) 
  {  
    casc = new vector<pair<int, double> >();
  };
  ~Cascade();

  // Create a cascade
  void addInfection(int node, double time);
  // Allocate memory for each seperate cascade
  void createBuffer();
  // Set parameters for the inference algorithm
  void setParams(int dimension_, double OW_, double stepsize_, double epsilon_, int negative_sampling_size_);

  string getID() const { return strID; }; 
  int getLength() const { return casc->size(); };

  void getInfectionAt(int i, int& node, double& time) const;
  double getDuration(double& min_t, double& max_t);
  void normalizeTime(double scale);
  void alignTime(double start_t);

  void update();
  // the first sweep to compute H,G,P,Q
  void firstsweep();
  double logLikelihood();
  void print();

  void descent(atomic<double> ** matrix, int node, const Vec& vc);

  inline void release(Vec ** v);
};

#endif
