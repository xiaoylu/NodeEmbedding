#include "matrix.h"

// initlization of the A (influence), B (selectivity) matrices
void Matrix2::init(int n, int m)
{
  this->n = n;
  this->m = m;
  cout << "Creating two " << n << "*" << m << " matrix" << endl;
  A = new atomic<double> *[n];  B = new atomic<double> *[n];
  for (int i = 0;i < n; ++i)
  {
    A[i] = new atomic<double>[m];  B[i] = new atomic<double>[m];
    for (int j = 0; j < m - 1; ++j)
    {
      double r = ((double) rand() / (RAND_MAX));
      A[i][j] = 0.002 * r;
      B[i][j] = 0.002 * r; 
    }
    A[i][m - 1] = B[i][m - 1] = 0.01;
  }

  Bsum = new Vec(m);
  sumB();
}

void Matrix2::sumB()
{
  // sum up all the "B_v"s for later usage
  Bsum->clear();
  for (int i = 0; i < n; ++i)
    *Bsum += B[i];
}

Matrix2::~Matrix2() 
{
  for (int i = 0;i < n;++i)
  {
    delete A[i];
    delete B[i];
  }
  delete A;
  delete B;
  delete Bsum;
}

//======================UTILITY FUNCTIONS======================
bool pairCompare(const pair<int, double>& firstElem, const pair<int, double>& secondElem) 
{
  return firstElem.second > secondElem.second;
}

// return the top node with highest norms
void Matrix2::getMaxNorm(vector<pair<int,double> >& ret)
{
  for (int i = 0;i < n; ++i)
  {
    ret.push_back( make_pair(i, Vec(m, A[i]).norm()) );
  }
  sort(ret.begin(), ret.end(), pairCompare); 
}

// return the top node with highest norms
void Matrix2::getMostInfluencialSites(vector<pair<int,double> >& ret, int k)
{
  for (int i = 0;i < n; ++i)
  {
    ret.push_back( make_pair(i, Vec(m, A[i]).getCompAt(k) ) );
  }
  sort(ret.begin(), ret.end(), pairCompare); 
}

// output the matrix to a txt file
void Matrix2::writeMat(string path)
{
  ofstream outfile_A;
  outfile_A.open((path + "TopProp_A.csv").c_str(), ofstream::out);
  ofstream outfile_B;
  outfile_B.open((path + "TopProp_B.csv").c_str(), ofstream::out);

  cout << "Writing matrix A,B to";
  //=========header=======
  for (int j = 0; j < m; ++j)
  {
    if (j != 0) outfile_A << ",";
    outfile_A << j; 
    if (j != 0) outfile_B << ",";
    outfile_B << j;
  }
  outfile_A << "\n";
  outfile_B << "\n";
  //=========content=======
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      if (j != 0) outfile_A << ",";
      outfile_A << A[i][j]; 

      if (j != 0) outfile_B << ",";
      outfile_B << B[i][j];
    }
    outfile_A << "\n";
    outfile_B << "\n";
  }
  cout << "  to " << path << "TopProp_A.csv" << endl;
  cout << "  and " << path << "TopProp_B.csv" << endl;
}
