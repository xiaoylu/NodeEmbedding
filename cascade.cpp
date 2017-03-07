#include "cascade.h"

void Cascade::getInfectionAt(int i, int& node, double& time) const
{
  pair<int, double> mypair = casc->at(i);
  node = mypair.first;
  time = mypair.second;
}

double Cascade::getDuration(double& min_t, double& max_t)
{
  double mintime = DBL_MAX, maxtime = -1; 
  for (pair<int, double>& mypair : *casc)
  {
    mintime = mintime < mypair.second? mintime : mypair.second;
    maxtime = maxtime > mypair.second? maxtime : mypair.second;
  }
  max_t = maxtime;
  min_t = mintime;
  return maxtime - mintime;
}

void Cascade::normalizeTime(double scale)
{
  for (pair<int, double>& mypair : (*casc))
    get<1>(mypair) /= scale; // normalization
}

void Cascade::alignTime(double start_t)
{
  for (pair<int, double>& mypair : (*casc))
  {
    get<1>(mypair) -= start_t; // alignment 
  }
}

void Cascade::update()
{
  int node;
  double time;
  int n = Matrix2::instance().getN();
  double rate_sample = (double) (n) / negative_sampling_size; 

  // the first sweep to compute H,G,P,Q
  firstsweep();

  // compute the derivative: dL(A,B)/dBv
  int size = getLength();
  for (int v = 1; v < size; ++v) // start from 1
                                     // because node 0 is influenced by nobody.   
  {
    getInfectionAt(v, node, time);

    Vec dBv(
      *G[v] - (*H[v] * time) + ( *H[v] /  (*H[v] ^ Matrix2::instance().B[node]) )
    );

    descent(Matrix2::instance().B, node, dBv);
  }

  // Now, update the negative links
  for (int i = 1; i < negative_sampling_size; i++)
  {
    int node_a = rand() % n;
    Vec dA(dimension, Matrix2::instance().A[node_a]);
    for (int v = 0; v < size; ++v) 
    {
      getInfectionAt(v, node, time);
      descent(Matrix2::instance().B, node, dA * (time - OW) * rate_sample );
    }
  }

  firstsweep();

  // compute the derivative: dL(A,B)/dAu
  Vec lastterm(dimension);
  for (int u = size - 1; u >= 0 ; --u)  // reverse order, start from size-1
                                        // because the last node influences nobody 
  {
    getInfectionAt(u, node, time);
   
    Vec dAu(
            ((*P[u]) * time - *Q[u]) + lastterm
    );

    descent(Matrix2::instance().A, node, dAu);

    // update the accumulative term
    // summation from u to size-1
    lastterm += Vec(dimension, Matrix2::instance().B[node])
                / (*H[u] ^ Matrix2::instance().B[node] ) ;
  }

  // Now, update the negative links
  for (int i = 0; i < negative_sampling_size; i++)
  {
    int node_b = rand() % n;
    Vec dB(dimension, Matrix2::instance().B[node_b]);
    for (int v = 0; v < size; ++v) 
    {
      getInfectionAt(v, node, time);
      descent(Matrix2::instance().A, node, dB * (time - OW) * rate_sample );
    }
  }
}

void Cascade::firstsweep()
{
  Vec vh(dimension); Vec vg(dimension); 
  int node; double time;
  int size = casc->size();

  int i = 0;
  for (auto& mypair : *casc)
  {
    node = mypair.first;
    time = mypair.second;

    vh += Matrix2::instance().A[node];
    vg.accu(Matrix2::instance().A[node], time);

    *(H[i+1]) = vh;  *(G[i+1]) = vg;
    ++i;
  }
  // H 0 1   2    ...  size-1     size
  //   * A 0 A 0,1a    A 0,size-2 A all 

  Vec vp(dimension); Vec vq(dimension); 
  for (i = size - 1;i >= 0;--i)
  {
    node = casc->at(i).first;
    time = casc->at(i).second;

    vp += Matrix2::instance().B[node];
    vq.accu(Matrix2::instance().B[node], time);

    if (i > 0) // the first node is not influenced by anyone.
    {
      *(P[i-1]) = vp;  *(Q[i-1]) = vq;
    }
    else if(i == 0)
    {
      *(P[size]) = vp;  *(Q[size]) = vq;
    }
  }
  // P 0             1            ...   size-2   size-1     size
  //   B 1,size-1    B 2,size-1         B size-1 None       B all
}

// projected gradient descent
void Cascade::descent(atomic<double> ** matrix, int node, const Vec& vc)
{
  for (int i = 0;i < dimension - 1; ++i) // the last component is always fixed
  {
    //read
    double tmp_read = matrix[node][i].load(memory_order_relaxed);
    tmp_read += (stepsize * vc.vec[i]);
    if (tmp_read <= epsilon) tmp_read = 0;
    // write
    matrix[node][i].store(tmp_read, memory_order_relaxed);
  }
}

double Cascade::logLikelihood()
{
  double ret = 0;

  firstsweep();
  
  int size = getLength();
  int node; double time;
  for (int v = 1; v < size ; ++v) // start from 1
  {
    getInfectionAt(v, node, time);

#ifdef DEBUG_MODE
    cout << node << "\t" << time << "*" << endl;
    cout << "H" << v;  H[v]->print();
    cout << "G" << v;  G[v]->print();
#endif

    double hold = *H[v] ^ Matrix2::instance().B[node];
    double hold2 = *G[v] ^ Matrix2::instance().B[node];

#ifdef DEBUG_MODE
    cout << hold << " " << hold2 << endl;
#endif

    ret += ( hold2 - hold * time + log( hold ) ) ; 
  }

  // the remaining nodes are uninfected by the time OW
  //if (ret > 0)
  //{
  //  cout << "**************************prev:" << ret << endl; 
  //  //print();
  //}
  //cout << casc->size() << endl;

  ret = ret + (
      (*(Matrix2::instance().Bsum) - *(P[size])) 
      ^ ( *(G[size]) - *(H[size])  * OW )
              );
  return ret;
}

void Cascade::setParams(int dimension_, double OW_, double stepsize_, double epsilon_, int negative_sampling_size_)
{ 
  dimension = dimension_;
  OW = OW_;
  stepsize = stepsize_;
  epsilon = epsilon_;
  negative_sampling_size = negative_sampling_size_;
}

Cascade::~Cascade()
{ 
  release(H); release(G);
  release(P); release(Q);
  delete casc; 
} 

void Cascade::addInfection(int node, double time )
{
  casc->push_back( make_pair(node, time) );
}

void Cascade::createBuffer()
{
  int size = getLength();
  H = new Vec*[size + 1]; 
  G = new Vec*[size + 1];
  P = new Vec*[size + 1];
  Q = new Vec*[size + 1];
  for (int i = 0;i < size + 1; ++i)
  {
    H[i] = new Vec(dimension); G[i] = new Vec(dimension);
    P[i] = new Vec(dimension); Q[i] = new Vec(dimension);
  }
} 

void Cascade::print()
{
  cout << "===" << strID << "===" << endl;
  for (auto& it : *casc)
  {
    cout << it.first << ":" << it.second << " ";
  }
  cout << endl;
}

// release memory 
inline void Cascade::release(Vec ** v)
{
  int size = getLength();
  if (!v || size <= 0) return;
  for (int i = 0; i < size; ++i)
  {
    delete v[i];
  }
  delete[] v;
  v = NULL;
}
