#ifndef HELP_H
#define HELP_H

#include <random>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include "headers.h"

class sample{

  int low;
  Index seed;
  std::mt19937 re;

 public:
  sample (Index seed, int l=0) : low(l), re(seed) {}
  int operator() (const std::initializer_list<double>& il)
  {
    std::discrete_distribution<> dist{il};
    return dist(re) + low;
  }
};

template<class F, size_t n> array<double, n> apply(array<double,n>& a, F f)
{
  for (auto& x : a)
    x = f(x);  
  return a;
}


template<size_t n> void operator += (array<double,n>& a, double val)
{
  for (auto& x : a)
    x += val;
}

template<size_t n> array<double,n> operator += (array<double,n>& a, const array<double,n>& b)
{
  for (int i=0; i<n; ++i)
    a[i] += b[i];
  return a;
}
  
template<size_t n> array<double,n> operator + (array<double,n> a, array<double,n> b)
{
  array<double,n> res;

  for (int i=0; i<n; ++i)
    res[i] = a[i] + b[i];
    
  return res;
}

template<class Cont, class Pred> std :: vector<int> find_idx_if (Cont& C, Pred P){

  std :: vector<int> idx;
  auto sz = C.size();
  
  for (int i=0; i<sz; ++i)
    if (P(C[i])) 
      idx.push_back(i);

  return idx;
}

template<class V, class F> Index count_if(const V& v, F f)
{
  return std::count_if(v.begin(), v.end(), f);
}

template<class V, class W> std :: vector<int> find_difference(V a, W b)
{ 
  std :: vector<int> res(a.size());
  auto p = std::set_difference(a.begin(), a.end(), b.begin(), b.end(), res.begin());
  res.erase(p, res.end());
  
  return res;
}


template<class Cont1, class Cont2>
std :: vector<int> match(Cont1& v1, Cont2& v2){

  std :: vector<int> res;

  for (auto& x : v2)
    {
      // find the index of x in v1
      for (int i=0; i<v1.size(); ++i)
	{
	  if (v1[i] == x)
	    {
	      res.push_back(i);
	      break;
	    }
	}
    }
  return res;
}

template<class F, class V> V transform(const V& v, F f)
{
  V res(v.size());
  std::transform(v.begin(), v.end(), res.begin(), f);
  return res;
}

template<class F, class V1, class V2> void transform (const V1& v1, V2& v2, F f)
{
  std::transform(v1.begin(), v1.end(), v2.begin(), f);
}

template<typename V, typename T> bool is_in (const T& t, const V& v)
{
  return std::find(v.begin(), v.end(), t) != v.end(); 
}


class equals{
  int x;
public:
  equals(int val) : x(val) {}
  bool operator () (int t) { return x == t; }
};

class not_equals{
  int x;
  public:
  not_equals(int x) {this->x = x;}
  bool operator () (int t) { return x != t; }
};


struct x_to_logx{ void operator () (double& x){ x = log(x);} };

#endif
