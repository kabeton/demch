#include <utility>
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#pragma once

typedef std::pair<double, double> point;

struct grid {
  double h;
  std::vector<point> x;
  grid() {
    this->h = 0;
    this->x = {std::make_pair(0.0, 0.0)};
  };
  grid(double _h, double x0, double xn, double y0) {
    this->h = _h;
    int n = (int)((xn - x0)/_h);
    this->x.push_back(std::make_pair(x0, y0));
    for(int i = 1; i <= n; i++) {
      this->x.push_back(std::make_pair(x0 + i*h, 0));
    }
  }
};

double dist(grid x, grid y) {
  grid grid_max, grid_min;
  if(x.h > y.h) {
    grid_max = x; grid_min = y;
  } else {
    grid_max = y; grid_min = x;
  }

  int n = (int)((grid_max.x.back().first - grid_max.x[0].first)/grid_max.h);
  int n2 = (int)((grid_min.x.back().first - grid_min.x[0].first)/grid_min.h);
  int ratio = (int)n2/n;
  std::vector<double> d(n);
  for(int i = 0; i <= n; i++) {
    d.push_back(fabs(grid_max.x[i].second - grid_min.x[i*ratio].second));
  }
  return *std::max_element(d.begin(), d.end());
}

class solver 
{
private:
  double x0, xn, y0;
protected:
  std::function<double(double, double)> f;
public:
  solver(double _x0, double _xn, double _y0,
      std::function<double(double, double)> _f) {
    x0 = _x0; xn = _xn; y0 = _y0; f = _f;
  }

  virtual double step(int i, double h, point xpr) = 0;

  std::pair<grid, grid> h_sol(double h){
    grid xh(h, x0, xn, y0);
    compute(xh);
    grid xh2(h/2, x0, xn, y0);
    compute(xh2);
    return std::make_pair(xh2, xh);
  }

  grid eps_sol(double error) {
    double h = (xn - x0)/10;
    std::pair<grid, grid> p = h_sol(h);
    while (dist(p.first, p.second)/15 > error) {
      h = h / 2;
      p = h_sol(h);
    }
    return p.first;
  }

  std::vector<grid> eps_sol_saved(double error) {
    double h = (xn - x0)/10;
    std::vector<grid> results;
    std::pair<grid, grid> p = h_sol(h);
    results.push_back(p.second);
    results.push_back(p.first);
    while (dist(p.first, p.second) > error) {
      h = h / 2;
      p = h_sol(h);
      results.push_back(p.first);
    }
    return results;
  }

  void compute(grid &x) {
    int n = (int)((xn - x0)/x.h);
    for(int i = 0; i < n; i++) {
      double y1 = step(i, x.h, x.x[i]);
      x.x[i+1].second = y1;
    }
  }

  void demo_print(double h) {
    std::pair<grid, grid> p = h_sol(h);
    int nf = (int)((xn - x0)/p.first.h);
    int ns = (int)((xn - x0)/p.second.h);
    int n = std::min(nf, ns);
    std::cout << "Grid with step size h/2: " << std::endl;
    for(auto i = 0; i <= 2*n; i+=2) {
      std::cout << "x: " << p.first.x[i].first << "  y: "
        << p.first.x[i].second << std::endl;
    }
    std::cout << "Grid with step size h: " << std::endl;
    for(auto i = 0; i <= n; i++) {
      std::cout << "x: " << p.second.x[i].first << "  y: " 
        << p.second.x[i].second << std::endl;
    }
    std::cout << "Difference: " << std::endl;
    for(auto i = 0; i <= n; i++) {
      std::cout << "x: " << p.second.x[i].first << "  |y_h/2 - y_h|: " 
        << fabs(p.second.x[i].second - p.first.x[2*i].second) << std::endl;
    }
  }
};

class rk4: public solver 
{
public:
  rk4(double _x0, double _xn, double _y0, 
      std::function<double(double, double)> _f) : solver(_x0, _xn, 
        _y0, _f){};

  double step(int i, double h, point xpr) {
    double f1, f2, f3, f4;
    double xi = xpr.first, yi = xpr.second;
    f1 = f(xi, yi);
    f2 = f(xi + h/2, yi + h*f1/2);
    f3 = f(xi + h/2, yi + h*f2/2);
    f4 = f(xi + h, yi + h*f3);
    return yi + h/6*(f1 + 2*f2 + 2*f3 + f4);
  }
};
