#ifndef vec_header
#define vec_header

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <atomic>

using namespace std;

class Vec 
{
public:
  double * vec;
  int m;

  Vec(int m) { this->m = m; vec = new double[m]; clear();};

  ~Vec() { delete[] vec; };

  Vec(const Vec& v)
  { 
    this->m = v.m;
    vec = new double[m];
    for (int i = 0; i < m; ++i)
    {
      this->vec[i] = v.vec[i];
    };
  };

  Vec(int m, const atomic<double> *array)
  { 
    this->m = m;
    vec = new double[m];
    for (int i = 0; i < m; ++i)
    {
      this->vec[i] = array[i].load(memory_order_relaxed);
    };
    //note if array is dynamically constructed
    //then it should be released somewhere outside 
  };

  void clear()
  {
     for (int i = 0;i < m; ++i)
     {
       this->vec[i] = 0;
     }
  }


  int amax() const
  {
     int ret = 0;
     for (int i = 1;i < m; ++i)
     {
       if (this->vec[i] > this->vec[ret])
         ret = i;
     }
     return ret;
  }

  // for compute the accumulative vector H,G, P,Q..
  void accu(const atomic<double> * array,const double mul)
  {
     for (int i = 0;i < m; ++i)
     {
       double tmp = array[i].load(memory_order_relaxed);
       this->vec[i] += ( tmp * mul );
     }
  };

  Vec& operator+=(const atomic<double> *array)
  {
    for (int i = 0;i < m; ++i)
    {
      this->vec[i] += array[i].load(memory_order_relaxed);
    }
    return *this;
  };

  Vec& operator+=(const Vec& other) {
    if (this != &other) { 
       for (int i = 0;i < m; ++i)
       {
         this->vec[i] += other.vec[i];
       }
    }
    return *this;
  };

  // assignment
  Vec& operator=(const Vec& other) {
    if (this != &other) { 
       for (int i = 0;i < m; ++i)
       {
         this->vec[i] = other.vec[i];
       }
    }
    return *this;
  };

  // inner product
  double operator^(atomic<double> * array)
  {
    double ret = 0;
    for (int i = 0;i < m; ++i)
    {
      double tmp = array[i].load(memory_order_relaxed);
      ret += this->vec[i] * tmp;
    }
    return ret;
  };

  // inner product
  double operator^(const Vec& other)
  {
    double ret = 0;
    for (int i = 0;i < m; ++i)
    {
      ret += this->vec[i] * other.vec[i];
    }
    return ret;
  };

  Vec& operator-=(const Vec& other) {
    for (int i = 0;i < m; ++i)
    {
      this->vec[i] -= other.vec[i];
    }
    return *this;
  };

  Vec& operator*=(const double val)
  {
    for (int i = 0;i < m; ++i)
    {
      this->vec[i] *= val;
    }
    return *this;
  }; 

  // point-wise product
  Vec& operator%=(const Vec& other)
  {
    for (int i = 0;i < m; ++i)
    {
      this->vec[i] *= other.vec[i];
    }
    return *this;
  }; 

  void print()
  {
    cout << "<<";
    // for debug
    for (int i = 0; i < m; ++i)
    {
      if (i < m - 3 && i > 5) 
      { 
        cout << "...";
        i = m - 3;
      }
      cout << this->vec[i] << "\t" ;
    }
    cout << ">>" << m << endl; 
  };

  // the L2 norm of a vector
  double norm()
  {
    double sum_square = 0;
    for (int i = 0; i < m; ++i)
    {
      sum_square += vec[i] * vec[i];
    }
    return sqrt(sum_square);
  }

  double getCompAt(int i) { return vec[i]; }
};
 
// a vector addes another vector 
inline Vec operator+(Vec l, const Vec& r)
{
  l += r;
  return l;
}

// vector minus a vector 
inline Vec operator-(Vec l, const Vec& r)
{
  l -= r;
  return l;
}

// vector multiples a number 
inline Vec operator*(Vec l, const double r)
{
  l *= r;
  return l;
}

// vector divides a number 
inline Vec operator/(Vec l, const double r)
{
  l *= (1.0 / r);
  return l;
}

// point-wise product
inline Vec operator%(Vec l, const Vec& r)
{
  l %= r;
  return l;
}
#endif
