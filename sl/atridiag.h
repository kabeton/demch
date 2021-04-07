#include <cassert>
#include <cstdio>
#include <vector>
#include <functional>

#ifndef ATRIDIAG_H
#define ATRIDIAG_H 
template<typename T>
void update_container(T a, double lambda);

template<typename T>
class atridiag
{
private:
  T a, b, c;
  int L;
  std::vector<double> diag{0};
public:
  atridiag(int _L, T _a, T _b, T _c) 
  {
    L = _L;
    a = _a;
    b = _b;
    c = _c;
  }

  ~atridiag(){};

  void update(double lambda)
  {
    update_container(b, lambda);
  }

  void get_diag()
  {
    diag.resize(L + 1);
    diag[0] = 1;
    diag[1] = b[1];
    diag[2] = b[2] - (c[1] + a[1])*a[2]/b[1];
    for(int i = 3; i < L; i++) 
    {
      diag[i] = b[i] - c[i - 1]*a[i]/diag[i - 1];
    }
    diag[L] = 1;
  }

  double det()
  {
    get_diag();
    double ans = 1;
    for(int i = 1; i <= L; i++) 
    {
      ans *= diag[i];   
    }
    return ans;
  }
};

class model_container
{
public:
  int size;
  double h, lambda;
  char type;
  std::function<double(double)> p;
  std::function<double(double)> q;
public:
  model_container() 
  {
    size = 0;
  }

  model_container(char _type, double _h, double _lambda, int _size,
                  std::function<double(double)> _q,
                  std::function<double(double)> _p) 
  {
    size = _size;
    type =_type;
    h = _h;
    lambda =_lambda;
    q = _q;
    p = _p;
  }

  ~model_container(){};

  void set_lambda(double l) 
  {
    lambda = l;
  }

  virtual double operator[] (const int index) 
  {
    if(type == 'b') {
      if(!(index == 0 || index == size - 1)) return lambda/2*h*h - 2;
    }

    if(type == 'c') {
      if(index == 0) return 0;
    }

    if(type == 'a') {
      if(index == size - 1) return 0;
    }
    return 1;
  }
};

class function_container : public model_container
{
private:
public:
  function_container() 
  {
    size = 0;
  }

  function_container(char _type, double _h, double _lambda, int _size,
                    std::function<double(double)> _q,
                    std::function<double(double)> _p) 
  : model_container(_type, _h, _lambda, _size, _q, _p){};

  ~function_container() {};

  double operator[] (const int index) {
    if(type == 'a') 
    {
      if(index == size - 1) return 0;
      return 3*q(h*(index - 1)) + q(h*(index + 1));
    }

    if(type == 'b')
    {
      if(!(index == 0 || index == size - 1)) 
      {
        return 4*(h*h*lambda*p(h*index) - (q(h*(index + 1)) + q(h*(index - 1))));
      }
    }

    if(type == 'c') 
    {
      if(index == 0) return 0;
      return 3*q(h*(index + 1)) + q(h*(index - 1));
    }
    return 1;
  }
};

void update_container(model_container &a, double lambda) 
{
  a.set_lambda(lambda);
}

void update_container(function_container &a, double lambda) 
{
  a.set_lambda(lambda);
}
#endif