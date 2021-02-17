#ifndef SPLINE_H
#define SPLINE_H 1

#include "poly.h"

class Spline {
  std::vector<std::vector<double>> a;
  std::vector<double> grid;
public:
  Spline();

  Spline(Poly p, std::vector<double> x) {
    grid = x;
    Poly pd = p.derivative();
    for(int i = 0; i <= p.getdeg() - 1; i++) {
      double h = x[i+1] - x[i];
      std::vector<double> local(4);
      local[0] = (-1*pd(x[i+1])*x[i]*x[i]*x[i+1]*h + p(x[i+1])*x[i]*x[i]*(3*x[i+1] - x[i]))/(h*h*h) +
                 (p(x[i])*x[i+1]*x[i+1]*(x[i+1] - 3*x[i]) - pd(x[i])*x[i]*x[i+1]*x[i+1]*h)/(h*h*h);
      local[1] = (pd(x[i+1])*x[i]*(2*x[i+1] + x[i])*h - 6*(p(x[i+1]) - p(x[i]))*x[i]*x[i+1])/(h*h*h) + 
                 (pd(x[i])*x[i+1]*(x[i+1] + 2*x[i])*h)/(h*h*h);
      local[2] = (-1*pd(x[i+1])*h*(x[i+1] + 2*x[i]) + 3*(p(x[i+1]) - p(x[i]))*(x[i+1] + x[i]))/(h*h*h) -
                 (pd(x[i])*h*(x[i] + 2*x[i+1]))/(h*h*h);
      local[3] = (pd(x[i+1])*h - 2*(p(x[i+1]) - p(x[i])) + pd(x[i])*h)/(h*h*h);
      a.push_back(local);
    }
  }
  
  std::vector<double> get_ith_coef(int i) const {
    if(i < a.size()) return a[i];
    return {0};
  }

  int size() const {
    return a.size();
  }

  double value(double x) const {
    int i = 0;
    double start = grid[i];
    while(start < x) {
      start = grid[++i];
    }
    Poly c = Poly{get_ith_coef(i-1)};
    return c.value(x);
  }

  ~Spline() {};
};

std::ostream& operator<<(std::ostream &os, const Spline x) {
  for(int i = 0; i < x.size(); i++) {
    std::vector<double> c = x.get_ith_coef(i);
    os << "i=" << i << " " << c[0] << " " << c[1] << " x " << c[2] << " x^2 " << c[3] << " x^3" << std::endl;
  }
  return os;
}

#endif
