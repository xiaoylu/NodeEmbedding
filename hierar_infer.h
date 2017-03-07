#ifndef HIERAR_INFER_HEADER
#define HIERAR_INFER_HEADER

#include <omp.h>
#include <iostream>
#include <sstream>
#include "matrix.h"

/* Algorithm specific parameters */



class HierarInfer
{
private:
  // mapping the compact index of nodes to their raw IDs 
  map<int, string> node2str; 
  // mapping a node to the cascades it belongs to
  map<int, set<string> *> node2casc;
  // mapping a cascade ID to the corresponding C++ instance 
  map<string, Cascade *> casc_list; 

  // the keys of the casc_list;
  vector<string> casc_IDs; 

  
  // the total number of nodes involved
  int num_nodes;

  // the total number of cascades involved 
  int num_cascades;

  double OW;

  int max_iterations; 
  int negative_sampling_size; 

  void constructCascIDs();
  void loadData(string txtfile);
  double averageDuration();

public:
  HierarInfer(string txtfile);
  ~HierarInfer();

  void filter(double min_duration, double max_duration);
  void normalize();

  void infer(int dimension, double OW,
             double stepsize, double epsilon,
             int negative_sampling, int max_iterations,
             int early_stop);

  double loglikelihood();

  void output(string output_path);

  //-------------------
  void split(const string &s, char delim, vector<string> &elems);
};

#endif
